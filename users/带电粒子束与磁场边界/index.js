/**
 * PhysLab 插件 · 带电粒子束与磁场边界
 * ------------------------------------------------------------------
 * 在官方「电磁·带电粒子」之外，新增一个模型：
 *   beam-particle  带电粒子束（n 个）+ 圆形 / 方形磁场边界
 *
 * 特点：
 *   1. 粒子个数 1–12，按间距 d 垂直入射方向排成一束（也可角度散开）；
 *   2. 磁场只存在于边界内：圆（半径 R）或矩形（宽 w、高 h，中心 bx,by）；
 *   3. 场外不受力 → 匀速直线；场内 → 圆周运动（有 E 时叠加 E×B 漂移）；
 *   4. 出射点、偏转角、边界内粒子数实时给出，画面上标出射点。
 *
 * 物理：场内为解析解（圆周/漂移），出射时刻用「二分法」在边界上求根，
 *       场外为匀速直线 —— 不是数值积分，出射点位置可以拿去对答案。
 */

PhysLab.definePlugin({
  id: 'beam-field',
  name: '带电粒子束与磁场边界',
  version: '1.0.1',
  author: 'PhysLab 用户插件',
  category: 'users',
  kinds: ['module', 'chart'],
  apiVersion: '1.0',

  activate(api) {
    const H = api.helpers, kit = api.kit;
    const { fmt, nfix, clamp, sArrow, sText, sHUD, FONT } = H;
    const TWO_PI = 2 * Math.PI;
    const p = (key, sym, name, unit, def, x) => Object.assign({ key, sym, name, unit, def }, x || {});

    /* ---------- 边界 ---------- */
    const hasBoundary = P => (P.R > 0) || (P.w > 0 && P.h > 0);
    const insideOf = (P, x, y) => {
      if (P.R > 0) return (x - P.bx) * (x - P.bx) + (y - P.by) * (y - P.by) <= P.R * P.R;
      if (P.w > 0 && P.h > 0) return Math.abs(x - P.bx) <= P.w / 2 && Math.abs(y - P.by) <= P.h / 2;
      return true;
    };
    const bboxOf = P => {
      if (P.R > 0) return [P.bx - P.R, P.bx + P.R, P.by - P.R, P.by + P.R];
      if (P.w > 0 && P.h > 0) return [P.bx - P.w / 2, P.bx + P.w / 2, P.by - P.h / 2, P.by + P.h / 2];
      return null;
    };

    /* ---------- 场内 / 场外解析推进 ---------- */
    function prop(P, s, t0, t1, inField) {
      const dt = t1 - t0;
      const om = P.q * P.B / P.m;
      if (!inField || Math.abs(P.B) < 1e-12) {          // 场外，或只有电场
        if (!inField) return { x: s.x + s.vx * dt, y: s.y + s.vy * dt, vx: s.vx, vy: s.vy };
        const a = P.q * P.E / P.m;
        return { x: s.x + s.vx * dt, y: s.y + s.vy * dt + 0.5 * a * dt * dt, vx: s.vx, vy: s.vy + a * dt };
      }
      if (Math.abs(om) < 1e-9) {                        // 只有电场
        const a = P.q * P.E / P.m;
        return { x: s.x + s.vx * dt, y: s.y + s.vy * dt + 0.5 * a * dt * dt, vx: s.vx, vy: s.vy + a * dt };
      }
      const vd = P.E / P.B, u0x = s.vx - vd, u0y = s.vy;
      const c = Math.cos(om * dt), si = Math.sin(om * dt);
      return {
        x: s.x + vd * dt + (u0x * si + u0y * (1 - c)) / om,
        y: s.y + (u0y * si - u0x * (1 - c)) / om,
        vx: vd + u0x * c + u0y * si,
        vy: u0y * c - u0x * si,
      };
    }
    const speed = s => Math.hypot(s.vx, s.vy);

    /* ---------- 求下一次穿过边界的时刻（二分） ---------- */
    function nextCross(P, s, t0, inField) {
      if (!hasBoundary(P)) return P.tMax;
      const K = 160;
      for (let i = 1; i <= K; i++) {
        const t = t0 + (P.tMax - t0) * i / K;
        const e = prop(P, s, t0, t, inField);
        if (insideOf(P, e.x, e.y) !== inField) {
          let a = t0 + (P.tMax - t0) * (i - 1) / K, b = t;
          for (let k = 0; k < 46; k++) {
            const m = (a + b) / 2, em = prop(P, s, t0, m, inField);
            if (insideOf(P, em.x, em.y) === inField) a = m; else b = m;
          }
          return (a + b) / 2;
        }
      }
      return P.tMax;
    }

    /* ---------- 单个粒子的分段表 ---------- */
    function segsOf(P, ini) {
      const segs = [];
      let s = Object.assign({}, ini), t = 0, inField = insideOf(P, s.x, s.y);
      for (let g = 0; g < 12 && t < P.tMax - 1e-12; g++) {
        let tc = nextCross(P, s, t, inField);
        if (tc <= t + 1e-9) tc = Math.min(P.tMax, t + 1e-6);   // 擦边：判据退化时给一个极小步长
        const t1 = Math.min(tc, P.tMax);
        segs.push({ t0: t, t1, inField, s0: Object.assign({}, s) });
        const e = prop(P, s, t, t1, inField);
        s = Object.assign({}, s, e); t = t1;
        if (tc >= P.tMax - 1e-12) break;
        inField = !inField;
      }
      if (t < P.tMax - 1e-12) segs.push({ t0: t, t1: P.tMax, inField, s0: Object.assign({}, s) });
      if (!segs.length) segs.push({ t0: 0, t1: P.tMax, inField, s0: Object.assign({}, s) });
      return segs;
    }

    /* ---------- 粒子初值 ---------- */
    function initParticles(P) {
      const n = clamp(Math.round(P.n), 1, 12), ph = P.phi * Math.PI / 180;
      const out = [];
      for (let i = 0; i < n; i++) {
        const off = (i - (n - 1) / 2) * P.d;
        out.push({
          i, off,
          x: P.x0 - Math.sin(ph) * off,
          y: P.y0 + Math.cos(ph) * off,
          vx: P.v0 * Math.cos(ph),
          vy: P.v0 * Math.sin(ph),
        });
      }
      return out;
    }

    function particleAt(P, segs, t) {
      let g = segs[segs.length - 1];
      for (const s of segs) { if (t <= s.t1 + 1e-12) { g = s; break; } }
      const e = prop(P, g.s0, g.t0, Math.min(Math.max(t, g.t0), g.t1), g.inField);
      e.v = speed(e); e.inside = g.inField ? 1 : 0;
      return e;
    }

    /* ---------- 主模拟 ---------- */
    function simulateBeam(P) {
      const inits = initParticles(P), n = inits.length;
      const parts = inits.map(q => segsOf(P, q));
      const T = Math.abs(P.q * P.B / P.m) > 1e-12 ? TWO_PI * P.m / Math.abs(P.q * P.B) : Infinity;
      const Rc = Math.abs(P.q * P.B / P.m) > 1e-12 ? P.m * P.v0 / Math.abs(P.q * P.B) : Infinity;

      const stateAt = t => {
        const ps = parts.map(segs => particleAt(P, segs, t));
        const p0 = ps[0];
        return {
          t, x: p0.x, y: p0.y, vx: p0.vx, vy: p0.vy, v: p0.v,
          Ek: 0.5 * P.m * p0.v * p0.v, inside: p0.inside,
          nkin: ps.reduce((a, q) => a + q.inside, 0), parts: ps,
        };
      };

      /* 出射点：每颗粒子第一次“场内 → 场外”的位置 */
      const exits = [];
      parts.forEach((segs, i) => {
        for (const g of segs) {
          if (g.inField && g.t1 < P.tMax - 1e-9) {
            const e = prop(P, g.s0, g.t0, g.t1, true);
            const ang = Math.atan2(e.vy, e.vx) * 180 / Math.PI;
            exits.push({ i, x: e.x, y: e.y, t: g.t1, ang, v: speed(e) });
            break;
          }
        }
      });

      const bb = bboxOf(P) || [0, 0, 0, 0];
      let xmin = bb[0], xmax = bb[1], ymin = bb[2], ymax = bb[3];
      for (let i = 0; i <= 60; i++) {
        const st = stateAt(P.tMax * i / 60);
        for (const q of st.parts) {
          xmin = Math.min(xmin, q.x); xmax = Math.max(xmax, q.x);
          ymin = Math.min(ymin, q.y); ymax = Math.max(ymax, q.y);
        }
      }
      if (!(xmax > xmin)) { xmin -= 1; xmax += 1; }
      if (!(ymax > ymin)) { ymin -= 1; ymax += 1; }
      const pad = 0.08 * Math.max(xmax - xmin, ymax - ymin);
      xmin -= pad; xmax += pad; ymin -= pad; ymax += pad;

      const phases = [];
      const g0 = parts[0];
      g0.forEach((g, k) => {
        const a = prop(P, g.s0, g.t0, g.t0, g.inField), b = prop(P, g.s0, g.t0, g.t1, g.inField);
        phases.push({
          t0: g.t0, t1: g.t1,
          title: g.inField ? (Math.abs(P.B) < 1e-12 ? '场区内：类平抛' : '场区内：圆周运动') : '场区外：匀速直线',
          lines: [
            '位置：(' + fmt(a.x) + ', ' + fmt(a.y) + ') → (' + fmt(b.x) + ', ' + fmt(b.y) + ') m',
            '速度：(' + fmt(a.vx) + ', ' + fmt(a.vy) + ') → (' + fmt(b.vx) + ', ' + fmt(b.vy) + ') m/s',
            '速率：' + fmt(speed(a)) + ' → ' + fmt(speed(b)) + ' m/s',
          ],
          end: k === g0.length - 1 ? '到达最长模拟时间' : '穿过边界',
        });
      });

      const st0 = stateAt(0), st1 = stateAt(P.tMax);
      const sim = kit.makeSim({
        P, tEnd: P.tMax,
        phases, samples: 60,
        state: t => stateAt(t),
        xmin, xmax,
        tiles: [
          ['粒子个数', String(n)],
          ['回旋半径 R', isFinite(Rc) ? fmt(Rc) + ' m' : '—'],
          ['回旋周期 T', isFinite(T) ? fmt(T) + ' s' : '—'],
          ['已出场', String(exits.length) + ' / ' + n],
        ],
        names: { v: 'v', vx: 'vx', vy: 'vy' },
        rows: () => {
          const R = [
            ['粒子个数 n', String(n)],
            ['回旋半径 R = mv / (|q|B)', isFinite(Rc) ? fmt(Rc) + ' m' : '—'],
            ['回旋周期 T = 2πm / (|q|B)', isFinite(T) ? fmt(T) + ' s' : '—'],
            ['场区形状', P.R > 0 ? '圆形（R = ' + fmt(P.R) + ' m）' : (P.w > 0 && P.h > 0 ? '方形（' + fmt(P.w) + ' × ' + fmt(P.h) + ' m）' : '无边界')],
            ['边界中心', '(' + fmt(P.bx) + ', ' + fmt(P.by) + ') m'],
            ['末时刻场内粒子数', String(st1.nkin) + ' / ' + n],
            ['出场粒子数', String(exits.length) + ' / ' + n],
          ];
          exits.forEach(e => R.push(['第 ' + (e.i + 1) + ' 个粒子出射', '(' + fmt(e.x) + ', ' + fmt(e.y) + ') m，t = ' + fmt(e.t) + ' s，速度方向 ' + fmt(e.ang) + '°']));

          return R;
        },
        checks: () => {
          const out = [];
          let dmax = 0;
          for (let i = 0; i <= 80; i++) {
            const st = stateAt(P.tMax * i / 80);
            if (Math.abs(P.E) < 1e-12 && Math.abs(P.B) > 1e-12) {
              for (const q of st.parts) dmax = Math.max(dmax, Math.abs(q.v - P.v0));
            }
          }
          if (Math.abs(P.E) < 1e-12) out.push({ status: dmax < 1e-6 * Math.max(1, P.v0) ? 'pass' : 'fail', item: '只有磁场时速率不变（洛伦兹力不做功）', detail: '全束最大速率偏差 ' + dmax.toExponential(1) + ' m/s' });
          if (hasBoundary(P)) {
            let bad = 0;
            parts.forEach(segs => segs.forEach(g => {
              if (g.t1 < P.tMax - 1e-9) {
                const e = prop(P, g.s0, g.t0, g.t1, g.inField);
                const onEdge = P.R > 0
                  ? Math.abs(Math.hypot(e.x - P.bx, e.y - P.by) - P.R)
                  : Math.min(Math.abs(Math.abs(e.x - P.bx) - P.w / 2), Math.abs(Math.abs(e.y - P.by) - P.h / 2));
                if (onEdge > 1e-6) bad++;
              }
            }));
            out.push({ status: bad ? 'fail' : 'pass', item: '每段端点都落在边界上（出射时刻求根）', detail: bad ? bad + ' 个端点偏离边界' : '全部端点误差 < 1e-6 m' });
          }
          out.push({ status: 'info', item: '建模假设', detail: '边界处磁场突变，忽略边缘效应；粒子视为质点、忽略重力与相对论效应。' });
          return out;
        },
        solution: () => {
          const st = [];
          st.push({
            title: '受力与运动判定', items: [
              { tag: '洛伦兹力', tex: '$f = |q|vB$，方向由左手定则，始终垂直于速度、不做功' },
              { tag: '场外', tex: '边界外无磁场 → 匀速直线运动' },
            ],
          });
          if (isFinite(Rc)) st.push({
            title: '场内圆周运动', items: [
              { tag: '向心力', tex: '$|q|vB = m v^2 / R$' },
              { tag: '半径', tex: '$R = m v / (|q|B) = #' + fmt(Rc) + '#$ m' },
              { tag: '周期', tex: '$T = 2πm / (|q|B) = #' + fmt(T) + '#$ s（与速率无关）' },
            ],
          });
          if (exits.length) {
            const e0 = exits[0];
            st.push({
              title: '出射点（几何法）', items: [
                { tag: '出射点', tex: '$(' + fmt(e0.x) + ', ' + fmt(e0.y) + ')$ m' },
                { tag: '出射速度方向', tex: '$' + fmt(e0.ang) + '°$（与 +x 轴夹角）' },
                { tag: '方法', tex: '先定轨迹圆心，再求轨迹圆与边界的交点，最后按几何求偏转角' },
              ],
            });
          }
          return st;
        },
        infoBox: st => '粒子个数 <b>' + st.parts.length + '</b>（场内 <b>' + st.nkin + '</b>）<br>第 1 个粒子：(' +
          nfix(st.x, 3) + ', ' + nfix(st.y, 3) + ') m<br>速率 <b>' + nfix(st.v, 3) + ' m/s</b>' +
          (isFinite(Rc) ? '<br>回旋半径 <b>' + nfix(Rc, 3) + ' m</b>' : ''),
        render2D: (ctx, W, Hh, s, st, o) => {
          const vw = (o && o.view) || {};
          const sx0 = s.xmin, sx1 = s.xmax, sy0 = s.ymin, sy1 = s.ymax;
          const sc = (vw.zoom || 1) * Math.min((W - 70) / (sx1 - sx0), (Hh - 90) / (sy1 - sy0));
          const cx = W / 2 + (vw.x || 0) - ((sx0 + sx1) / 2) * sc;
          const cy = Hh / 2 + 8 + (vw.y || 0) + ((sy0 + sy1) / 2) * sc;
          const X = x => cx + x * sc, Y = y => cy - y * sc;

          ctx.fillStyle = '#fffdf8'; ctx.fillRect(0, 0, W, Hh);

          /* 边界 */
          if (hasBoundary(P)) {
            ctx.save();
            ctx.beginPath();
            if (P.R > 0) ctx.arc(X(P.bx), Y(P.by), P.R * sc, 0, TWO_PI);
            else ctx.rect(X(P.bx - P.w / 2), Y(P.by + P.h / 2), P.w * sc, P.h * sc);
            ctx.fillStyle = 'rgba(28,126,214,.055)'; ctx.fill();
            ctx.setLineDash([7, 5]); ctx.lineWidth = 1.8; ctx.strokeStyle = 'rgba(28,126,214,.75)'; ctx.stroke();
            ctx.setLineDash([]);
            ctx.restore();
          }

          /* 磁场符号（只在边界内） */
          if (Math.abs(P.B) > 1e-12) {
            ctx.strokeStyle = ctx.fillStyle = 'rgba(28,126,214,.32)'; ctx.lineWidth = 1.2;
            const gap = 34;
            for (let gx = 24; gx < W; gx += gap) {
              for (let gy = 52; gy < Hh - 8; gy += gap) {
                const wx = (gx - cx) / sc, wy = (cy - gy) / sc;
                if (!insideOf(P, wx, wy)) continue;
                if (P.B > 0) { ctx.beginPath(); ctx.arc(gx, gy, 2.2, 0, TWO_PI); ctx.fill(); ctx.beginPath(); ctx.arc(gx, gy, 5.6, 0, TWO_PI); ctx.stroke(); }
                else { ctx.beginPath(); ctx.moveTo(gx - 3.6, gy - 3.6); ctx.lineTo(gx + 3.6, gy + 3.6); ctx.moveTo(gx + 3.6, gy - 3.6); ctx.lineTo(gx - 3.6, gy + 3.6); ctx.stroke(); }
              }
            }
          }

          /* 坐标轴 */
          ctx.strokeStyle = '#ded7c6'; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
          ctx.beginPath(); ctx.moveTo(0, Y(0)); ctx.lineTo(W, Y(0)); ctx.moveTo(X(0), 44); ctx.lineTo(X(0), Hh); ctx.stroke();
          ctx.setLineDash([]);

          /* 每颗粒子的轨迹与当前位置 */
          const N = 140, colors = ['#e8590c', '#1c7ed6', '#2b8a3e', '#7048e8', '#d6336c', '#0c8599', '#a13a0a', '#5c940d', '#6741d9', '#c2255c', '#0b7285', '#e67700'];
          const n = st.parts.length;
          for (let i = 0; i < n; i++) {
            const c = colors[i % colors.length];
            const segs = parts[i];
            ctx.strokeStyle = c; ctx.lineWidth = 1.9; ctx.beginPath();
            for (let k = 0; k <= N; k++) {
              const e = particleAt(P, segs, st.t * k / N);
              k ? ctx.lineTo(X(e.x), Y(e.y)) : ctx.moveTo(X(e.x), Y(e.y));
            }
            ctx.stroke();
          }
          for (let i = 0; i < n; i++) {
            const q = st.parts[i], c = colors[i % colors.length];
            ctx.fillStyle = c; ctx.beginPath(); ctx.arc(X(q.x), Y(q.y), 6.5, 0, TWO_PI); ctx.fill();
            ctx.strokeStyle = '#fffdf8'; ctx.lineWidth = 1.6; ctx.stroke();
          }

          /* 出射点 */
          ctx.lineWidth = 1.6;
          exits.forEach(e => {
            if (e.t > st.t + 1e-9) return;
            ctx.strokeStyle = 'rgba(224,49,49,.85)';
            ctx.beginPath(); ctx.arc(X(e.x), Y(e.y), 4, 0, TWO_PI); ctx.stroke();
          });

          /* 第 1 个粒子的速度矢量 */
          if (o && o.showVec !== false && st.v > 1e-9) {
            sArrow(ctx, X(st.x), Y(st.y), X(st.x) + st.vx / Math.max(s.vmax || st.v, 1e-9) * 58, Y(st.y) - st.vy / Math.max(s.vmax || st.v, 1e-9) * 58, '#2b8a3e', 2.4);
          }

          sText(ctx, (P.R > 0 ? '圆形边界 R = ' + nfix(P.R, 2) + ' m' : (P.w > 0 && P.h > 0 ? '方形边界 ' + nfix(P.w, 2) + '×' + nfix(P.h, 2) + ' m' : '无边界')) +
            '　B ' + (Math.abs(P.B) < 1e-12 ? '= 0' : (P.B > 0 ? '垂直纸面向外' : '垂直纸面向内')), W - 12, 22, '#566078', '12px ' + FONT, 'right');
          sHUD(ctx, W, Hh, s, st);
        },
      });
      sim.xmin = xmin; sim.xmax = xmax; sim.ymin = ymin; sim.ymax = ymax;
      sim.vmax = Math.max(st1.v, P.v0, 1e-6);
      sim.beam = { parts, exits, Rc, T, n };
      return sim;
    }

    /* ---------- 注册模型 ---------- */
    api.models.register({
      id: 'beam-particle',
      name: '电磁·带电粒子束与磁场边界',
      short: '粒子束',
      carrier: '粒子',
      noCarrier: true,
      segLabel: '粒子束',
      stats: [['t', 't', 's'], ['x', 'x', 'm'], ['y', 'y', 'm'], ['vx', 'vx', 'm/s'], ['vy', 'vy', 'm/s'], ['v', 'v', 'm/s'], ['Ek', 'Ek', 'J'], ['nkin', '场内个数', '']],
      defaultCharts: ['beam_xy', 'beam_vt'],
      aiNote: 'x 向右、y 向上；B>0 垂直纸面向外。R>0 为圆形边界，否则 w、h 均>0 时为方形边界；场外不受力做匀速直线运动。',
      params: [
        p('m', 'm', '粒子质量', 'kg', '1', { group: '粒子', pos: true }),
        p('q', 'q', '电荷量（带符号）', 'C', '1', { group: '粒子', ai: '负电荷填负数' }),
        p('v0', 'v<sub>0</sub>', '初速率', 'm/s', '2', { group: '粒子', pos: true }),
        p('phi', 'φ', '入射方向（与 +x 夹角）', '°', '0', { group: '粒子', hint: '0 = 沿 +x；90 = 沿 +y' }),

        p('n', 'n', '粒子个数', '个', '5', { group: '粒子束', pos: true, hint: '1–12，取整' }),
        p('d', 'd', '相邻粒子间距', 'm', '0.5', { group: '粒子束', nonneg: true, hint: '沿垂直于入射方向排开' }),
        p('x0', 'x<sub>0</sub>', '束起点 x', 'm', '-3', { group: '粒子束' }),
        p('y0', 'y<sub>0</sub>', '束起点 y', 'm', '0', { group: '粒子束' }),

        p('B', 'B', '磁感应强度（带符号）', 'T', '1', { group: '场', hint: '正 = 垂直纸面向外，0 = 无磁场', ai: '正 = 向外，负 = 向内' }),
        p('E', 'E', '电场强度（沿 y，可选）', 'V/m', '0', { group: '场', hint: '一般填 0；非 0 时叠加 E×B 漂移' }),

        p('R', 'R', '圆形边界半径', 'm', '3', { group: '磁场边界', nonneg: true, hint: 'R>0 用圆形边界；填 0 表示不用圆形' }),
        p('w', 'w', '方形边界宽', 'm', '0', { group: '磁场边界', nonneg: true, hint: 'R=0 时，w 与 h 都>0 才用方形边界' }),
        p('h', 'h', '方形边界高', 'm', '0', { group: '磁场边界', nonneg: true }),
        p('bx', 'x<sub>c</sub>', '边界中心 x', 'm', '0', { group: '磁场边界' }),
        p('by', 'y<sub>c</sub>', '边界中心 y', 'm', '0', { group: '磁场边界' }),

        p('tMax', 't<sub>max</sub>', '模拟时间', 's', '6', { group: '环境', pos: true }),
      ],
      build(cfg, P, errors) {
        if (!errors.q && Math.abs(P.q) < 1e-15) errors.q = '电荷量不能为 0';
        if (!errors.w && !errors.h && P.R <= 0 && (P.w > 0) !== (P.h > 0)) errors.h = '方形边界要同时填 w 和 h';
      },
      simulate: simulateBeam,
      presets: [
        {
          name: '圆形边界 · 平行束入射（会聚 / 散开）',
          params: { m: '1', q: '1', v0: '2', phi: '0', n: '5', d: '0.5', x0: '-3', y0: '0', B: '1', E: '0', R: '3', w: '0', h: '0', bx: '0', by: '0', tMax: '6' },
          problem: {
            text: '半径 R=3 m 的圆形区域内存在垂直纸面向外、B=1 T 的匀强磁场。质量 m=1 kg、电荷量 q=1 C 的粒子束以 v₀=2 m/s 沿 +x 方向从左侧射入，束宽 2 m（5 个粒子，间距 0.5 m）。求：\n(1) 单颗粒子的回旋半径；\n(2) 各粒子的出射点位置；\n(3) 说明为什么这束粒子出射后不再平行。',
            answer: '(1) R = mv/(qB) = 2 m  (2) 见「结果」页出射点表  (3) 不同入射弦长对应不同偏转角，出射后方向不同',
          },
        },
        {
          name: '圆形边界 · 增大 B（轨迹变紧）',
          params: { m: '1', q: '1', v0: '2', phi: '0', n: '3', d: '0.4', x0: '-3', y0: '0', B: '2', E: '0', R: '3', w: '0', h: '0', bx: '0', by: '0', tMax: '6' },
          problem: {
            text: '圆形磁场区域半径 R=3 m，B=2 T 垂直纸面向外。粒子束以 v₀=2 m/s 沿 +x 射入，m=1 kg、q=1 C。求轨迹半径，并与 B=1 T 时的出射点比较。',
            answer: 'R = mv/(qB) = 1 m；B 加倍、轨迹半径减半，出射点更靠近入射侧（见结果页）',
          },
        },
        {
          name: '方形边界 · 束流穿过（有场/无场对照）',
          params: { m: '1', q: '1', v0: '2', phi: '0', n: '4', d: '0.4', x0: '-3', y0: '0', B: '1', E: '0', R: '0', w: '4', h: '3', bx: '0', by: '0', tMax: '6' },
          problem: {
            text: '边长 4 m × 3 m 的方形区域内存在垂直纸面向外、B=1 T 的匀强磁场（区域外无磁场）。粒子束以 v₀=2 m/s 沿 +x 射入，m=1 kg、q=1 C。求粒子在区域内外的运动性质，并判断哪些粒子能穿出区域。',
            answer: '区域内做半径 2 m 的圆周运动，区域外匀速直线；出射情况见「结果」页',
          },
        },
        {
          name: '无边界对照 · 完整圆周运动',
          params: { m: '1', q: '1', v0: '2', phi: '0', n: '3', d: '1', x0: '0', y0: '0', B: '1', E: '0', R: '0', w: '0', h: '0', bx: '0', by: '0', tMax: '7' },
          problem: {
            text: '无边界匀强磁场 B=1 T，三颗粒子以 v₀=2 m/s 沿 +x 射入，m=1 kg、q=1 C，忽略重力。求轨迹半径、周期，并说明它们的轨迹关系。',
            answer: 'R = 2 m，T = 2π s；三颗粒子轨迹为等半径圆，圆心沿 y 方向平移',
          },
        },
      ],
    });

    /* ---------- 图表 ---------- */
    api.charts.register({
      id: 'beam_xy', name: 'y-x 轨迹（第 1 个粒子）', xl: 'x / m', yl: 'y / m', models: ['beam-particle'],
      series: [{ n: '第 1 个粒子', c: '#e8590c', x: 'x', y: 'y' }],
    });
    api.charts.register({
      id: 'beam_vt', name: 'v-t 速率（磁场中不变）', xl: 't / s', yl: 'v / (m·s⁻¹)', time: true, models: ['beam-particle'],
      series: [{ n: 'v', c: '#e8590c', x: 't', y: 'v' }, { n: 'vx', c: '#1c7ed6', x: 't', y: 'vx' }, { n: 'vy', c: '#2b8a3e', x: 't', y: 'vy' }],
    });
  },
});
