/* engine/inference.py 的 **JS 镜像**（纯函数、无 DOM）。
 *
 * 用途：界面按需调用「从部分作答预测某模块特质的最终得分 + 置信度」，
 * 不必把 Python 搬进浏览器。参数由 `python engine/inference.py --dump-params
 * <CODE> --out params.json` 导出（含 `pairB64`），两边用**同一份参数**。
 *
 * ⚠️ 两个实现必须逐值一致（本项目历史上「两套实现漂移」出过问题）。
 * 一致性由 `tools/check/test_inference_pyjs_consistency.py` 实测：
 * 同一批合成部分作答，逐项比对 `predicted_score` / `confidence_r` /
 * `r2_total`，并打印最大偏差。
 *
 * 模型（与 Python 逐行对应，别只改一边）
 * ======================================
 *   T_known = 已答的**本维度**题项之和（已知）
 *   T_U     = 未答的本维度题项之和（要预测）
 *   pred    = (T_known + Σ_{j∈U} μⱼ + σ_U · cᵀ Σ_AA⁻¹ z_A) / scale
 *   cᵢ      = (Σ_{j∈U} σⱼ ρᵢⱼ) / σ_U
 *   confidence_r（判据口径）= √(bᵀ Σ_AA⁻¹ b)，b = corrected 题↔维度相关，
 *              预测变量**只含本维度已答的题**（与 dim_ceiling 同口径）
 *   r2_total（预测准确度口径）= 1 − σ_U²(1 − R²_U)/σ_T²
 *
 * 输入/输出结构见 `predict()` 的注释；与 `engine/inference.py` 的 `predict()` 同形。
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.EmoHealthInference = factory();
}(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  var TOL = 1e-9;      // 与 Python 的 TOL 相同
  var RIDGE = 1e-6;    // 与 Python 的 RIDGE 相同（仅在解失败时加在对角线上）

  /* base64 → Float32 小端（不依赖 Buffer，浏览器里也能跑） */
  function decodePairB64(b64) {
    var bytes;
    if (typeof Buffer !== 'undefined' && typeof Buffer.from === 'function') {
      var buf = Buffer.from(b64, 'base64');
      bytes = new Uint8Array(buf.buffer, buf.byteOffset, buf.length);
    } else {
      var bin = (typeof atob === 'function' ? atob(b64) : '');
      bytes = new Uint8Array(bin.length);
      for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    }
    var n = bytes.length >> 2;
    var dv = new DataView(bytes.buffer, bytes.byteOffset, n * 4);
    var out = new Float64Array(n);
    for (var k = 0; k < n; k++) out[k] = dv.getFloat32(k * 4, true);
    return out;
  }

  /* 展平三角阵下标：a=max(i,j)、b=min(i,j)，**1-based 题序** */
  function triIndex(a, b) {
    if (a < b) { var t = a; a = b; b = t; }
    return (a - 1) * (a - 2) / 2 + (b - 1);
  }

  /* 从 pairB64 建 n×n 相关矩阵（行主序展平）。下标是**量表包题序**。 */
  function buildPairMatrix(pairB64, n) {
    var raw = decodePairB64(pairB64);
    var m = new Float64Array(n * n);
    for (var i = 0; i < n; i++) m[i * n + i] = 1;
    for (var a = 1; a <= n; a++) {
      for (var b = 1; b < a; b++) {
        var v = raw[triIndex(a, b)];
        m[(a - 1) * n + (b - 1)] = v;
        m[(b - 1) * n + (a - 1)] = v;
      }
    }
    return m;
  }

  /* 解 A x = b（行主序方阵）；近奇异时加岭重试。与 Python 的 `_solve` 同策略。 */
  function solve(a, b, n) {
    var m = new Float64Array(n * (n + 1));
    var i, j, k;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) m[i * (n + 1) + j] = a[i * n + j];
      m[i * (n + 1) + n] = b[i];
    }
    for (var c = 0; c < n; c++) {
      var piv = c;
      for (var r = c + 1; r < n; r++) {
        if (Math.abs(m[r * (n + 1) + c]) > Math.abs(m[piv * (n + 1) + c])) piv = r;
      }
      if (Math.abs(m[piv * (n + 1) + c]) < 1e-9) return null;
      if (piv !== c) {
        for (k = 0; k <= n; k++) {
          var tmp = m[c * (n + 1) + k];
          m[c * (n + 1) + k] = m[piv * (n + 1) + k];
          m[piv * (n + 1) + k] = tmp;
        }
      }
      for (r = 0; r < n; r++) {
        if (r === c) continue;
        var f = m[r * (n + 1) + c] / m[c * (n + 1) + c];
        if (!f) continue;
        for (k = c; k <= n; k++) m[r * (n + 1) + k] -= f * m[c * (n + 1) + k];
      }
    }
    var x = new Float64Array(n);
    for (i = 0; i < n; i++) {
      var d = m[i * (n + 1) + i];
      x[i] = Math.abs(d) < 1e-12 ? 0 : m[i * (n + 1) + n] / d;
    }
    return x;
  }

  function solveWithRidge(a, b, n) {
    var x = solve(a, b, n);
    if (x) return x;
    var d = new Float64Array(a);
    for (var i = 0; i < n; i++) d[i * n + i] += RIDGE;
    return solve(d, b, n);
  }

  function clamp01(v) { return Math.max(0, Math.min(1, v)); }

  /* 与 engine/scoring.py 的 `DimensionScore.to_dict()` 对齐：官方分四舍五入到 4 位。
     （不这么做，Python 侧 1.3571 与 JS 侧 1.3571428… 会被一致性测试判为偏差。） */
  function round4(v) { return Math.round((v + Number.EPSILON) * 1e4) / 1e4; }

  function subMatrix(m, n, idx) {
    var k = idx.length, out = new Float64Array(k * k);
    for (var i = 0; i < k; i++) {
      for (var j = 0; j < k; j++) out[i * k + j] = m[idx[i] * n + idx[j]];
    }
    return out;
  }

  /* ── Cholesky：前瞻要**一次分解、多次回代**，不能每个右端项都做一次消元 ──
     `Σ_AA` 是对称（半）正定矩阵，Cholesky 分解 O(m³/3)，之后每个右端项只要 O(m²)。
     若分解失败（近奇异）→ 加岭 1e-6 重试；再失败返回 null，调用方走慢路径
     （直接调 `dimConfidence` 精确重算）。与 `solveWithRidge` 同一策略。 */
  function cholesky(a, m, ridge) {
    var L = new Float64Array(m * m);
    var r = ridge || 0;
    for (var i = 0; i < m; i++) {
      for (var j = 0; j <= i; j++) {
        var sum = a[i * m + j] + (i === j ? r : 0);
        for (var k = 0; k < j; k++) sum -= L[i * m + k] * L[j * m + k];
        if (i === j) {
          if (!(sum > 0) || !isFinite(sum)) return null;
          L[i * m + i] = Math.sqrt(sum);
        } else {
          L[i * m + j] = sum / L[j * m + j];
        }
      }
    }
    return L;
  }

  function cholSolveVec(L, m, b) {
    var y = new Float64Array(m), x = new Float64Array(m);
    for (var i = 0; i < m; i++) {
      var s = b[i];
      for (var k = 0; k < i; k++) s -= L[i * m + k] * y[k];
      y[i] = s / L[i * m + i];
    }
    for (i = m - 1; i >= 0; i--) {
      var t = y[i];
      for (k = i + 1; k < m; k++) t -= L[k * m + i] * x[k];
      x[i] = t / L[i * m + i];
    }
    return x;
  }

  function dot(a, b) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  /* 原始作答 → 计分值（反向题 (min+max)−v）。与 engine/scoring.py 同规则。 */
  function codeResponses(params, responses) {
    var lo = params.response_min, hi = params.response_max;
    var out = {};
    for (var iid in responses) {
      if (!Object.prototype.hasOwnProperty.call(responses, iid)) continue;
      var key = params.key_of_item[iid];
      if (key === undefined) throw new Error('量表外题项：' + iid);
      var v = responses[iid];
      if (!(lo <= v && v <= hi)) throw new Error(iid + ' 的作答 ' + v + ' 超出 ' + lo + '..' + hi);
      out[iid] = key >= 0 ? v : (lo + hi) - v;
    }
    return out;
  }

  /* 官方分（缺失不猜）：只有答满才给值；聚合方式取量表声明。 */
  function officialScore(params, dim, coded) {
    var own = params.dim_items[dim];
    var vals = [], missing = 0;
    for (var i = 0; i < own.length; i++) {
      if (coded[own[i]] === undefined) { missing++; continue; }
      vals.push(coded[own[i]]);
    }
    var agg = params.dim_agg[dim];
    if (missing) {
      return { value: null, complete: false, n_answered: vals.length,
               n_items: own.length, aggregation: agg,
               blocker: '官方分缺失不猜：' + vals.length + '/' + own.length + ' 题已答' };
    }
    var v = null;
    if (agg === 'sum' || agg === 'weighted_sum') {
      v = round4(vals.reduce(function (s, x) { return s + x; }, 0));
    } else if (agg === 'mean') {
      v = round4(vals.reduce(function (s, x) { return s + x; }, 0) / vals.length);
    } else if (agg === 'max') {
      v = round4(Math.max.apply(null, vals));
    } else if (agg === 'count_true') {
      v = vals.filter(function (x) { return x > 0; }).length;
    }
    return { value: v, complete: true, n_answered: vals.length,
             n_items: own.length, aggregation: agg, blocker: null };
  }

  /* 本维度题的统计量（题级 SD、题号集合、完整分 SD）——`dimConfidence` 与
     `forecastGains` 共用同一套，保证两处口径不会漂移。 */
  function ownStats(params, pairs, n, dim, idxOf) {
    var own = params.dim_items[dim];
    var s = new Float64Array(own.length);
    for (var j = 0; j < own.length; j++) s[j] = params.item_sd[own[j]];
    var io = own.map(function (id) { return idxOf[id]; });
    var sigmaT = 0;
    for (var x1 = 0; x1 < own.length; x1++) {
      for (var x2 = 0; x2 < own.length; x2++) {
        sigmaT += s[x1] * s[x2] * pairs[io[x1] * n + io[x2]];
      }
    }
    var ownSet = {};
    for (j = 0; j < own.length; j++) ownSet[own[j]] = true;
    return { own: own, s: s, io: io, ownSet: ownSet,
             sigmaT: Math.sqrt(Math.max(sigmaT, TOL)) };
  }

  /* 一道题在某个维度判据里的系数 cᵢ（**与作答值无关**）：
       · 本题在该维度内 → corrected（去本题）
       · 本题在别的维度 → 与完整分的普通相关                             */
  function cForItem(params, pairs, n, dim, itemId, idxOf, st) {
    if (st.ownSet[itemId]) {
      var row = params.item_dim_corr[itemId] || {};
      return typeof row[dim] === 'number' ? row[dim] : 0;
    }
    var num = 0, base = idxOf[itemId] * n;
    for (var j = 0; j < st.own.length; j++) num += st.s[j] * pairs[base + idxOf[st.own[j]]];
    return num / st.sigmaT;
  }

  /* 判据口径 confidence_r：**可借跨维度信息**（与 engine/inference.py 的 dim_confidence 逐行对应）。
   *
   *   c_i = corr(x_i, T_f)          i 是别的维度的题（x_i 不在 T_f 里 → 普通相关）
   *   c_i = corr(x_i, T_f − x_i)    i 是本维度已答的题（corrected，去本题）
   *   confidence_r = √(cᵀΣ_AA⁻¹c)
   *
   * 硬规则：本维度一道题都没答 → status=NO_OWN_ITEM_ANSWERED、confidence_r=null。 */
  function dimConfidence(params, pairs, n, dim, coded, idxOf) {
    var own = params.dim_items[dim];
    var aOwn = own.filter(function (id) { return coded[id] !== undefined; });
    if (!aOwn.length) {
      return { dimension: dim, status: 'NO_OWN_ITEM_ANSWERED',
               confidence_r: null, confidence_r_own_only: null,
               cross_dim_borrow_gain: null, n_own_answered: 0 };
    }
    var answered = params.ids.filter(function (id) { return coded[id] !== undefined; });
    var st = ownStats(params, pairs, n, dim, idxOf);
    var c = new Float64Array(answered.length);
    for (var k = 0; k < answered.length; k++) {
      c[k] = cForItem(params, pairs, n, dim, answered[k], idxOf, st);
    }
    var ia = answered.map(function (id) { return idxOf[id]; });
    var xb = solveWithRidge(subMatrix(pairs, n, ia), c, answered.length);
    var rAll = xb ? Math.sqrt(clamp01((function () {
      var t = 0; for (var q = 0; q < c.length; q++) t += c[q] * xb[q]; return t;
    })())) : 0;
    // 旧口径（只用本维度已答、corrected）——只作对照，不用于门控
    var iaOwn = aOwn.map(function (id) { return idxOf[id]; });
    var b = new Float64Array(aOwn.length);
    for (var j = 0; j < aOwn.length; j++) {
      var r2 = params.item_dim_corr[aOwn[j]] || {};
      b[j] = typeof r2[dim] === 'number' ? r2[dim] : 0;
    }
    var xbOwn = solveWithRidge(subMatrix(pairs, n, iaOwn), b, aOwn.length);
    var rOwn = xbOwn ? Math.sqrt(clamp01((function () {
      var t = 0; for (var q = 0; q < b.length; q++) t += b[q] * xbOwn[q]; return t;
    })())) : 0;
    return { dimension: dim, status: 'OK', confidence_r: rAll,
             confidence_r_own_only: rOwn, cross_dim_borrow_gain: rAll - rOwn,
             n_own_answered: aOwn.length, n_pred: answered.length };
  }

  /* ══ 信息增益前瞻：下一题推荐的核心 ═══════════════════════════════════════

     给一个**候选题** x，预测「答了 x 之后各维度的判据置信度 conf_r 变成多少」，
     取收益最大者推荐。**只依赖题间相关结构（pairB64）与「当前答了哪些题」**：
     签名里**根本没有作答值**（`answeredIds` / `candidateIds` 都是题号数组），
     所以不可能用到「用户还没答的题的作答值」。

     数学（为什么不用重算一遍 conf_r）—— 本项目「推荐 = 信息增益」的数学依据
     ----------------------------------------------------------------------
     记 A = 已答集合、Σ_AA 为其相关矩阵、某维度 f 的系数向量 b（与 `dimConfidence`
     同一套：本维度题用 corrected 去本题、别的维度用与完整分的普通相关），则
         conf_r_f(A)² = bᵀ Σ_AA⁻¹ b
     对候选题 x 做秩一增广 Σ' = [[Σ_AA, s], [sᵀ, 1]]、c' = [b; c_x]，其中 s = Σ_{A,x}。
     用分块求逆（Schur 补）：
         Σ'⁻¹ = [[Σ_AA⁻¹ + (Σ_AA⁻¹s)(Σ_AA⁻¹s)ᵀ/d,  −Σ_AA⁻¹s/d],
                 [−(Σ_AA⁻¹s)ᵀ/d,                    1/d]]
         d = 1 − sᵀΣ_AA⁻¹s
     于是（记 t = Σ_AA⁻¹s）：
         c'ᵀΣ'⁻¹c' = bᵀΣ_AA⁻¹b + (bᵀt − c_x)² / d
     即 **精确**增量
         Δ_f(x) = (bᵀt − c_x)² / (1 − sᵀt)
         conf_r_f(A∪{x})² = conf_r_f(A)² + Δ_f(x)
     Δ 就是「答这题对该维度最终分的**解释方差增量**」——标准意义的信息增益（R² 增量），
     而不是另立一套评分。注意它**不需要**知道 x 的作答值：相关系数结构已经决定了
     「答 x 能把该维度估多准」。

     ⚠️ 退化分支的含义：`1 − sᵀt → 0` 表示候选题与已答集合**近乎共线**（x 的作答几乎
     能被已答题线性预测出来）⇒ 该题的信息增量在数学上趋近 0（分子分母同时趋零，
     极限存在但浮点上 0/0 不可靠）。此时**不猜**：把该候选题退回 `dimConfidence`
     精确重算，并计入 `n_slow`；调用方据此判断长表会不会退化成慢路径。

     复杂度（task-18 提的慢点：m 大时逐题重算会明显卡）
     ------------------------------------------------
     朴素做法：每个 (候选题 × 维度) 重解一次 m×m 方程 → O(C·D·m³)。
     本实现：**Σ_AA 只 Cholesky 分解一次** O(m³/3)，之后
       · 每个维度一次回代 → v_d = Σ_AA⁻¹b_d（O(D·m²)）
       · 每个候选题一次回代 → t = Σ_AA⁻¹s（O(C·m²)）
       · 每个 (候选题, 维度) 只剩 O(m) 点积 + O(1) 公式
     总计 O(m³/3 + (D+C)·m² + C·D·m)。

     效用（**产品层**，与上面的统计量分开写清）
     ----------------------------------------
       utility = Σ_d w_d·Δ_d(x) + UNLOCK_BONUS · #{d : d 因答 x 首次跨过门槛}
       · w_d = 1（还不可给结果）/ W_AVAILABLE（已可给结果，只留一点精度价值）
       · 只统计「答完 x 后**本维度 ≥1 题已答**（即可展示）」的维度；本维度一题未答的
         维度**不参与打分、也不进推荐文案**（硬规则，与面板一致）。
     纯信息增益排序另存 `order_by_gain`，便于对照/审计。
  */
  var W_AVAILABLE = 0.25;
  var UNLOCK_BONUS = 1.0;

  function _now() {
    return (typeof performance !== 'undefined' && performance.now)
      ? performance.now() : Date.now();
  }

  function forecastGains(payload, answeredIds, candidateIds, opts) {
    opts = opts || {};
    var t0 = _now();
    var params = payload.params || payload;
    var n = params.ids.length;
    var pairs = payload.__pairs || (payload.__pairs = buildPairMatrix(payload.pairB64, n));
    var idxOf = {};
    for (var i = 0; i < n; i++) idxOf[params.ids[i]] = i;
    var ratio = (params.ceiling_ratio != null) ? params.ceiling_ratio : 0.9;
    var dims = Object.keys(params.dim_items);

    var seen = {};
    var answeredIn = answeredIds || [];
    for (i = 0; i < answeredIn.length; i++) seen[answeredIn[i]] = true;
    var A = params.ids.filter(function (id) { return seen[id] === true; });   // 量表内、去重、保序
    var m = A.length;
    var ia = A.map(function (id) { return idxOf[id]; });

    var M = new Float64Array(m * m);
    for (i = 0; i < m; i++) {
      for (var j = 0; j < m; j++) M[i * m + j] = pairs[ia[i] * n + ia[j]];
    }
    var L = m ? (cholesky(M, m, 0) || cholesky(M, m, RIDGE)) : null;

    /* 基线：逐维度的 b 向量（前瞻要复用它算 bᵀt）、r²、门槛、可展示性 */
    var base = {}, Bd = {};
    for (var di = 0; di < dims.length; di++) {
      var dim = dims[di];
      var st = ownStats(params, pairs, n, dim, idxOf);
      var b = new Float64Array(m);
      for (var k = 0; k < m; k++) b[k] = cForItem(params, pairs, n, dim, A[k], idxOf, st);
      var v = m ? (L ? cholSolveVec(L, m, b) : solveWithRidge(M, b, m)) : new Float64Array(0);
      var r2 = (m && v) ? Math.max(0, Math.min(1, dot(b, v))) : 0;
      var ownBefore = 0;
      for (k = 0; k < m; k++) if (st.ownSet[A[k]]) ownBefore++;
      var need = (params.dim_ceiling && params.dim_ceiling[dim] != null)
        ? params.dim_ceiling[dim] * ratio : null;
      var confBefore = Math.sqrt(r2);
      Bd[dim] = b;
      base[dim] = {
        dim: dim, st: st, r2: r2, need: need, n_own: ownBefore,
        n_items: st.own.length, display_allowed: ownBefore > 0,
        available: ownBefore > 0 && need != null && confBefore >= need,
      };
    }

    var out = {
      n_answered: m, n_candidates: 0, n_slow: 0, ratio: ratio,
      weights: { available: W_AVAILABLE, unlock_bonus: UNLOCK_BONUS },
      dims: {}, candidates: [], order_by_gain: [], order: [], best: null, ms: 0,
    };
    for (di = 0; di < dims.length; di++) {
      var b0 = base[dims[di]];
      out.dims[b0.dim] = {
        n_own_answered: b0.n_own, n_items: b0.n_items, need: b0.need,
        display_allowed: b0.display_allowed, available: b0.available,
        // 硬规则：本维度一题未答 → 连基线置信度都不给（面板同样是 blocked）
        confidence_r: b0.display_allowed ? Math.sqrt(b0.r2) : null,
      };
    }

    var tCache = {};
    for (var ci = 0; ci < candidateIds.length; ci++) {
      var cand = candidateIds[ci];
      var x = idxOf[cand];
      if (x === undefined || seen[cand]) continue;         // 量表外 / 已答 → 不进候选池
      var gains = {}, after = {}, unlocks = [], total = 0, ok = true;
      var dispTotal = 0, nDispTargets = 0, blockedBenef = 0;
      var tv = null, denom = 1;
      if (m) {
        if (!(cand in tCache)) {
          var s = new Float64Array(m);
          for (k = 0; k < m; k++) s[k] = pairs[x * n + ia[k]];
          var tt = L ? cholSolveVec(L, m, s) : solveWithRidge(M, s, m);
          tCache[cand] = (!tt || !(1 - dot(s, tt) > 1e-9)) ? null
            : { t: tt, stt: dot(s, tt) };
        }
        tv = tCache[cand];
        if (!tv) ok = false;                      // 退化：1−sᵀt 不可靠 → 走精确重算
        else denom = 1 - tv.stt;
      }
      for (di = 0; ok && di < dims.length; di++) {
        var dm = dims[di], bd = base[dm];
        if (!(bd.n_own + (bd.st.ownSet[cand] ? 1 : 0))) continue;  // 答完仍不可展示 → 不进推荐
        var cx = cForItem(params, pairs, n, dm, cand, idxOf, bd.st);
        var delta = m ? (Math.pow(dot(Bd[dm], tv.t) - cx, 2) / denom) : (cx * cx);
        if (!isFinite(delta) || delta < 0) { ok = false; break; }
        var r2A = Math.max(0, Math.min(1, bd.r2 + delta));
        var cA = Math.sqrt(r2A);
        var av = bd.need != null && cA >= bd.need;
        gains[dm] = delta;
        after[dm] = { confidence_r: cA, available: av,
                      newly_available: av && !bd.available };
        total += (bd.available ? W_AVAILABLE : 1) * delta;
        /* 展示口径：本维度**当前**可展示时，它的增益才可以显示成数字；
           未答维度（硬规则）只记「会受益」，不进可显示的合计 */
        if (bd.display_allowed) {
          dispTotal += delta;
          if (delta > 0) nDispTargets++;
        } else if (bd.st.ownSet[cand]) blockedBenef++;
        if (av && !bd.available) unlocks.push(dm);
      }
      var slow = !ok;
      if (slow) {
        /* 退化/分解失败：用**同一个** dimConfidence 精确重算（慢路径，逐题计数） */
        gains = {}; after = {}; unlocks = []; total = 0;
        dispTotal = 0; nDispTargets = 0; blockedBenef = 0;
        var coded = {};
        for (k = 0; k < m; k++) coded[A[k]] = 0;
        coded[cand] = 0;
        for (di = 0; di < dims.length; di++) {
          var dm2 = dims[di], bd2 = base[dm2];
          if (!(bd2.n_own + (bd2.st.ownSet[cand] ? 1 : 0))) continue;
          var cf = dimConfidence(params, pairs, n, dm2, coded, idxOf);
          var d2 = Math.max(0, cf.confidence_r * cf.confidence_r
                             - (bd2.display_allowed ? bd2.r2 : 0));
          var av2 = bd2.need != null && cf.confidence_r >= bd2.need;
          gains[dm2] = d2;
          after[dm2] = { confidence_r: cf.confidence_r, available: av2,
                         newly_available: av2 && !bd2.available };
          total += (bd2.available ? W_AVAILABLE : 1) * d2;
          if (bd2.display_allowed) {
            dispTotal += d2;
            if (d2 > 0) nDispTargets++;
          } else if (bd2.st.ownSet[cand]) blockedBenef++;
          if (av2 && !bd2.available) unlocks.push(dm2);
        }
      }
      var candDim = null, target = null, bestGain = -1;
      for (di = 0; di < dims.length; di++) {
        var dd = dims[di];
        if (base[dd].st.ownSet[cand]) candDim = candDim == null ? dd : candDim;
        if (gains[dd] != null && gains[dd] > bestGain) { bestGain = gains[dd]; target = dd; }
      }
      out.candidates.push({
        item: cand, own_dim: candDim, target: target,
        utility: total + UNLOCK_BONUS * unlocks.length,
        gain_total: total, unlocks: unlocks, gains: gains, after: after, slow: slow,
        // 展示口径（页面据此决定能不能把数字写出来，见 recoWhy）
        gain_total_displayable: dispTotal, n_displayable_targets: nDispTargets,
        blocked_beneficiaries: blockedBenef,
      });
    }
    out.n_candidates = out.candidates.length;
    var posOf = {};
    for (i = 0; i < params.ids.length; i++) posOf[params.ids[i]] = i;
    var byUtil = out.candidates.slice().sort(function (p, q) {
      return (q.utility - p.utility) || (q.gain_total - p.gain_total)
        || (posOf[p.item] - posOf[q.item]);
    });
    var byGain = out.candidates.slice().sort(function (p, q) {
      return (q.gain_total - p.gain_total) || (posOf[p.item] - posOf[q.item]);
    });
    out.order = byUtil.map(function (c) { return c.item; });
    out.order_by_gain = byGain.map(function (c) { return c.item; });
    out.best = byUtil.length ? byUtil[0] : null;

    /* 冠军的**显示数值**用同一个 dimConfidence 复核（排名用闭式，显示同口径） */
    if (out.best) {
      var codedW = {};
      for (k = 0; k < m; k++) codedW[A[k]] = 0;
      codedW[out.best.item] = 0;
      var confirmed = {};
      for (di = 0; di < dims.length; di++) {
        var dmw = dims[di], bdw = base[dmw];
        if (!(bdw.n_own + (bdw.st.ownSet[out.best.item] ? 1 : 0))) continue;
        var cw = dimConfidence(params, pairs, n, dmw, codedW, idxOf);
        var avw = bdw.need != null && cw.confidence_r >= bdw.need;
        confirmed[dmw] = {
          before: bdw.display_allowed ? Math.sqrt(bdw.r2) : null,
          after: cw.confidence_r, available: avw,
          gain: Math.max(0, cw.confidence_r * cw.confidence_r
                           - (bdw.display_allowed ? bdw.r2 : 0)),
          newly_available: avw && !bdw.available,
        };
      }
      out.best.confirmed = confirmed;
    }
    out.n_slow = out.candidates.filter(function (c) { return c.slow; }).length;
    out.ms = _now() - t0;
    return out;
  }

  /* 预测入口。payload = Python `--dump-params` 的输出：
   *   { params: <params_for_js()>, pairB64: "..." }
   * responses = { item_id: 原始作答 }
   * 返回结构与 engine/inference.py 的 predict() 同形（official / inferred 分层）。 */
  function predict(payload, responses) {
    var params = payload.params, n = params.ids.length;
    var pairs = payload.__pairs || buildPairMatrix(payload.pairB64, n);
    var idxOf = {};
    for (var i = 0; i < n; i++) idxOf[params.ids[i]] = i;
    var coded = codeResponses(params, responses);
    var answered = params.ids.filter(function (id) { return coded[id] !== undefined; });

    var out = {
      instrument: params.instrument,
      model: {
        name: 'conditional_expectation_linear',
        formula: 'pred = T_known + Σ_{j∈U} μⱼ + σ_TU·cᵀΣ_AA⁻¹z_A；' +
                 'confidence_r = √(cᵀΣ_AA⁻¹c)（本维度已答用 corrected，' +
                 '别的维度用与完整分的普通相关——可借跨维度）',
        pair_source: 'pairB64（与 engine/inference.py 同一份）',
        item_stats_source: params.stats_source,
        ceiling_ratio: params.ceiling_ratio,
        ceiling_ratio_cross: params.ceiling_ratio_cross,
        ceiling_ratio_source: params.ceiling_ratio_source,
        ceiling_display_denominator: 'ceiling_cross（跨维度上限＝整卷答满时的 conf_r）；' +
          '判据仍用 ceiling × ceiling_ratio，两者不要混',
        hard_rule: '本维度一道题都没答 → 不展示置信度与预测（NO_OWN_ITEM_ANSWERED）'
      },
      n_answered_total: answered.length,
      official: {},
      inferred: {},
      layer_note: 'official 为量表声明口径（缺失不猜）；inferred 为预测，带 inferred=true，不得混算。'
    };

    var dims = Object.keys(params.dim_items);
    for (var d = 0; d < dims.length; d++) {
      var dim = dims[d];
      var own = params.dim_items[dim];
      var scale = params.dim_agg[dim] === 'mean' ? own.length : 1;
      var aOwn = own.filter(function (id) { return coded[id] !== undefined; });
      var u = own.filter(function (id) { return coded[id] === undefined; });
      var ceil = params.dim_ceiling[dim];
      var ceilCross = params.dim_ceiling_cross ? params.dim_ceiling_cross[dim] : null;
      var need = (typeof ceil === 'number') ? ceil * params.ceiling_ratio : null;
      var q;

      out.official[dim] = officialScore(params, dim, coded);

      var block = {
        dimension: dim, name: params.dim_name[dim], inferred: true,
        n_answered: answered.length, n_answered_own: aOwn.length,
        n_items: own.length, ceiling: ceil, ceiling_cross: ceilCross,
        ceiling_ratio: params.ceiling_ratio,
        ceiling_need: need, aggregation: params.dim_agg[dim]
      };
      // ---- 硬规则：本维度一道题都没答 → 禁止展示置信度与预测 ----
      if (!aOwn.length) {
        block.status = answered.length ? 'NO_OWN_ITEM_ANSWERED' : 'NO_DATA';
        block.display_allowed = false;
        block.available = false;
        block.predicted_score = null;
        block.confidence_r = null; block.r2 = null;
        block.r2_total = null; block.confidence_r_total = null;
        block.confidence_r_own_only = null; block.cross_dim_borrow_gain = null;
        block.ceiling_pct = null; block.ceiling_pct_raw = null;
        block.prediction_se = null;
        block.known_own_sum = null; block.predicted_unknown_sum = null;
        block.predicted_unknown_shift = null; block.sigma_unknown = null;
        block.contributions = []; block.top_contributions = [];
        block.cross_dim_share_of_prediction = null;
        block.reason = answered.length
          ? '本维度一道题都没答：按用户硬规则禁止展示该维度的置信度与预测（与其它维度答了多少、r_total 多高无关）'
          : '整卷尚未作答：无任何输入';
        out.inferred[dim] = block;
        continue;
      }

      // ---- 判据口径 confidence_r：可借跨维度信息（本维度已答部分用 corrected） ----
      var conf = dimConfidence(params, pairs, n, dim, coded, idxOf);
      var confR = conf.confidence_r;

      // ---- 预测（已知的本维度题 + BLUP 未答题） ----
      var known = aOwn.reduce(function (s, id) { return s + coded[id]; }, 0);
      var contributions = [], predUSum = 0, predUShift = 0;
      var sigmaU = 0, sigmaSum = 0, r2Total = 1;
      if (!u.length) {
        // 答满：预测 = 官方分（精确）
      } else {
        var sSd = new Float64Array(own.length), sMu = new Float64Array(own.length);
        for (q = 0; q < own.length; q++) {
          sSd[q] = params.item_sd[own[q]];
          sMu[q] = params.item_mean[own[q]];
        }
        var pos = {};
        for (q = 0; q < own.length; q++) pos[own[q]] = q;
        var ku = u.map(function (id) { return pos[id]; });
        // cov = diag(s) · pairs[own,own] · diag(s)
        var covSum = 0, covUUSum = 0;
        for (var x1 = 0; x1 < own.length; x1++) {
          for (var x2 = 0; x2 < own.length; x2++) {
            var cv = sSd[x1] * sSd[x2] * pairs[idxOf[own[x1]] * n + idxOf[own[x2]]];
            covSum += cv;
            if (u.indexOf(own[x1]) >= 0 && u.indexOf(own[x2]) >= 0) covUUSum += cv;
          }
        }
        sigmaU = Math.sqrt(Math.max(covUUSum, TOL));
        var c = new Float64Array(answered.length);
        for (q = 0; q < answered.length; q++) {
          var ii = idxOf[answered[q]];
          var num = 0;
          for (var kk = 0; kk < ku.length; kk++) {
            num += sSd[ku[kk]] * pairs[ii * n + idxOf[own[ku[kk]]]];
          }
          c[q] = num / sigmaU;
        }
        var iaAll = answered.map(function (id) { return idxOf[id]; });
        var mAll = subMatrix(pairs, n, iaAll);
        var beta = solveWithRidge(mAll, c, answered.length);
        var aList = answered, cList = c;
        if (!beta) {                       // 保底：退回只用本维度已答题
          aList = aOwn;
          cList = new Float64Array(aOwn.length);
          for (q = 0; q < aOwn.length; q++) cList[q] = c[answered.indexOf(aOwn[q])];
          var iaOwn2 = aOwn.map(function (id) { return idxOf[id]; });
          beta = solveWithRidge(subMatrix(pairs, n, iaOwn2), cList, aOwn.length);
          if (!beta) beta = new Float64Array(aOwn.length);
        }
        var z = new Float64Array(aList.length);
        for (q = 0; q < aList.length; q++) {
          z[q] = (coded[aList[q]] - params.item_mean[aList[q]]) / (params.item_sd[aList[q]] || 1);
        }
        var muUSum = 0;
        for (q = 0; q < ku.length; q++) muUSum += sMu[ku[q]];
        var bz = 0;
        for (q = 0; q < aList.length; q++) bz += beta[q] * z[q];
        predUShift = sigmaU * bz;
        predUSum = muUSum + predUShift;
        var r2tu = 0;
        for (q = 0; q < aList.length; q++) r2tu += cList[q] * beta[q];
        r2tu = clamp01(r2tu);
        sigmaSum = Math.sqrt(Math.max(covSum, TOL));
        r2Total = clamp01(1 - (sigmaU * sigmaU) * (1 - r2tu) / (sigmaSum * sigmaSum));
        var ownSet = {};
        for (q = 0; q < own.length; q++) ownSet[own[q]] = true;
        for (q = 0; q < aList.length; q++) {
          var contrib = beta[q] * z[q] * sigmaU;
          contributions.push({
            item: aList[q], kind: ownSet[aList[q]] ? 'own' : 'cross',
            coded_value: coded[aList[q]], weight: beta[q], contribution: contrib
          });
        }
        contributions.sort(function (p1, p2) {
          return Math.abs(p2.contribution) - Math.abs(p1.contribution);
        });
      }
      var totalAbs = contributions.reduce(function (s, x) { return s + Math.abs(x.contribution); }, 0);
      var crossAbs = 0;
      for (q = 0; q < contributions.length; q++) {
        contributions[q].share = totalAbs ? Math.abs(contributions[q].contribution) / totalAbs : 0;
        if (contributions[q].kind === 'cross') crossAbs += Math.abs(contributions[q].contribution);
      }
      var predicted = (known + predUSum) / scale;
      var seSum = Math.sqrt(Math.max(0, sigmaSum * sigmaSum * (1 - r2Total)));
      block.status = 'OK';
      block.display_allowed = true;
      block.predicted_score = predicted;
      block.confidence_r = confR;
      block.confidence_r_own_only = conf.confidence_r_own_only;
      block.cross_dim_borrow_gain = conf.cross_dim_borrow_gain;
      // 「已达上限的 X%」：分母用**跨维度上限**（否则跨维度借用会 >100%）
      block.ceiling_pct_raw = ceilCross ? confR / ceilCross : null;
      block.ceiling_pct = ceilCross ? clamp01(confR / ceilCross) : null;
      block.r2 = confR * confR;
      block.r2_total = r2Total;
      block.confidence_r_total = Math.sqrt(r2Total);
      block.prediction_se = seSum / scale;
      block.available = !!(need !== null && confR >= need);
      block.known_own_sum = known;
      block.predicted_unknown_sum = u.length ? predUSum : 0;
      block.predicted_unknown_shift = u.length ? predUShift : 0;
      block.sigma_unknown = u.length ? sigmaU : 0;
      block.contributions = contributions;
      block.top_contributions = contributions.slice(0, 5);
      block.cross_dim_share_of_prediction = totalAbs ? crossAbs / totalAbs : 0;
      block.note = u.length ? '已知部分已计入；未答部分由已答题的条件期望给出'
                            : '答满该维度全部题：预测 = 官方分（精确）';
      out.inferred[dim] = block;
    }
    return out;
  }

  function predictFrom(paramsPayload, responses) {
    // 便于测试：允许传入未解码的参数包，内部缓存矩阵
    if (!paramsPayload.__pairs) {
      paramsPayload.__pairs = buildPairMatrix(paramsPayload.pairB64, paramsPayload.params.ids.length);
    }
    return predict(paramsPayload, responses);
  }

  return {
    decodePairB64: decodePairB64,
    buildPairMatrix: buildPairMatrix,
    triIndex: triIndex,
    codeResponses: codeResponses,
    dimConfidence: dimConfidence,
    forecastGains: forecastGains,
    predict: predict,
    predictFrom: predictFrom
  };
}));
