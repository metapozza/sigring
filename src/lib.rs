//! # sigring
//!
//! Exact, streaming, sliding-window **truncated path signatures** in a wrapping integer ring (`i64` = Z/2^64, `i128` = Z/2^128).
//!
//! * State is the factorial-cleared tensor `t_n = n! * S_n * K^n` of integer increments (scale `K` turns real data into integers).
//! * A window slide is one fused push + pop (`S' = E(-old) (x) S (x) E(new)`), division-free (Horner form), allocation-free.
//! * Wrapping arithmetic is a ring homomorphism, so the result is **exact whenever the true value fits**, even if intermediates wrapped.
//! * [`SigRing::guard_ok`] is a rigorous sufficient test: every entry satisfies `|t_m| <= L1^m` (L1 = window path length), so the state is
//!   valid if `L1^N < 2^(BITS-1)`. The bound is attained, i.e. tight.
//! * Results are bit-identical across platforms/compilers (integer arithmetic only).
//! * [`PrefixStore`] answers arbitrary interval queries exactly via `antipode(prefix_i) (x) prefix_j` (no float inverse).
//! * Log-signatures: [`LogSigPrep`] + [`SigRing::log_signature_expanded`] (log taken in `f64` at readout only, from the exact signature).
//!
//! ```
//! use sigring::SigRing;
//! let (d, n, w) = (2, 3, 4);
//! let mut eng: SigRing<i64> = SigRing::new(d, n);
//! let incs: Vec<[i64; 2]> = vec![[3, -1], [2, 2], [-4, 1], [5, 0], [1, 1], [0, -2]];
//! for i in 0..incs.len() {
//!     let out = if i >= w { incs[i - w] } else { [0, 0] };
//!     eng.ingest(&incs[i], &out);
//! }
//! assert!(eng.guard_ok());
//! let sig = eng.signature(1.0); // levels 1..=3 of the last 4 increments
//! assert_eq!(sig.len(), 2 + 4 + 8);
//! ```
pub trait Ring: Copy + Default + PartialEq {
    fn add(self, o: Self) -> Self; fn mul(self, o: Self) -> Self; fn neg(self) -> Self;
    fn from_i64(v: i64) -> Self; fn to_f64(self) -> f64; fn is_zero(self) -> bool; const BITS: u32;
}
macro_rules! impl_ring { ($t:ty, $u:ty) => { impl Ring for $t {
    #[inline(always)] fn add(self, o: Self) -> Self { self.wrapping_add(o) }
    #[inline(always)] fn mul(self, o: Self) -> Self { self.wrapping_mul(o) }
    #[inline(always)] fn neg(self) -> Self { self.wrapping_neg() }
    #[inline(always)] fn from_i64(v: i64) -> Self { v as $t }
    fn to_f64(self) -> f64 { let neg = self < 0; let m = self.unsigned_abs() as u128; let r = (m as u64) as f64 + ((m >> 64) as u64 as f64) * 18446744073709551616.0; if neg { -r } else { r } }
    #[inline(always)] fn is_zero(self) -> bool { self == 0 }
    const BITS: u32 = <$t>::BITS;
}}}
impl_ring!(i64, u64);
impl_ring!(i128, u128);
pub struct SigRing<T: Ring = i128> {
    pub d: usize,
    pub n: usize,
    off: Vec<usize>,      // offset of level k (1..=n) in flat state; off[k-1]
    pub t: Vec<T>,     // levels 1..=n concatenated
    pw: Vec<T>,        // scratch: powers of vector, same layout
    pl: Vec<T>,
    binom: Vec<Vec<T>>,
    /// sliding sum over the window of the L1 norm of increments (exact, small); used for the rigorous overflow guard
    pub l1: i128,
    hs: HornerScratch<T>,
    lbuf: Vec<T>,
    rbuf: Vec<T>,
}
#[inline(always)]
fn outer_acc<T: Ring>(coef: T, a: &[T], b: &[T], out: &mut [T]) {
    let lb = b.len();
    for (ia, &x) in a.iter().enumerate() {
        let cx = coef.mul(x);
        let row = &mut out[ia * lb..(ia + 1) * lb];
        for (o, &y) in row.iter_mut().zip(b.iter()) { *o = o.add(cx.mul(y)); }
    }
}
impl<T: Ring> SigRing<T> {
    pub fn new(d: usize, n: usize) -> Self {
        let mut off = vec![0usize; n + 1];
        let mut sz = 1usize;
        let mut acc = 0usize;
        for k in 1..=n { sz *= d; off[k - 1] = acc; acc += sz; off[k] = acc; }
        let total = acc;
        let mut binom = vec![vec![T::default(); n + 1]; n + 1];
        for i in 0..=n { binom[i][0] = T::from_i64(1); for j in 1..=i { binom[i][j] = binom[i - 1][j - 1].add(if j <= i - 1 { binom[i - 1][j] } else { T::default() }); } }
        Self { d, n, off, t: vec![T::default(); total], pw: vec![T::default(); total], pl: vec![T::default(); total], binom, l1: 0, hs: HornerScratch::new(d, n), lbuf: vec![T::default(); d], rbuf: vec![T::default(); d] }
    }
    #[inline(always)] fn lvl(&self, k: usize) -> std::ops::Range<usize> { self.off[k - 1]..self.off[k] }
    fn powers(off: &Vec<usize>, n: usize, v: &[T], buf: &mut [T]) {
        let d = v.len();
        buf[off[0]..off[1]].copy_from_slice(v);
        for k in 2..=n {
            let (lo, hi) = buf.split_at_mut(off[k - 1]);
            let prev = &lo[off[k - 2]..off[k - 1]];
            let cur = &mut hi[..off[k] - off[k - 1]];
            for (ia, &x) in prev.iter().enumerate() { for ib in 0..d { cur[ia * d + ib] = x.mul(v[ib]); } }
        }
    }
    /// Push increment `inc` (right-multiply) and pop `out` (left-multiply by inverse) in one fused update (division-free Horner kernel).
    pub fn ingest(&mut self, inc: &[i64], out: &[i64]) {
        let mut l1 = 0i128;
        for c in 0..self.d { self.rbuf[c] = T::from_i64(inc[c]); self.lbuf[c] = T::from_i64(out[c]).neg(); l1 += (inc[c] as i128).abs() - (out[c] as i128).abs(); }
        self.l1 += l1;
        let l = std::mem::take(&mut self.lbuf); let r = std::mem::take(&mut self.rbuf);
        self.ingest_i128(&l, &r);
        self.lbuf = l; self.rbuf = r;
    }
    /// Same update with the lower-level (ring-element) interface: Horner kernel.
    pub fn ingest_i128(&mut self, l: &[T], r: &[T]) {
        let mut sc = std::mem::replace(&mut self.hs, HornerScratch { a: Vec::new(), b: Vec::new() });
        self.ingest_horner(l, r, &mut sc);
        self.hs = sc;
    }
    /// Reference kernel (power tensors + outer products). Kept for cross-checking; bit-identical to the Horner kernel.
    pub fn ingest_naive(&mut self, inc: &[i64], out: &[i64]) {
        let r: Vec<T> = inc.iter().map(|&x| T::from_i64(x)).collect();
        let l: Vec<T> = out.iter().map(|&x| T::from_i64(x).neg()).collect();
        self.l1 += inc.iter().map(|&x| (x as i128).abs()).sum::<i128>() - out.iter().map(|&x| (x as i128).abs()).sum::<i128>();
        self.ingest_naive_ring(&l, &r);
    }
    pub fn ingest_naive_ring(&mut self, l: &[T], r: &[T]) {
        let n = self.n;
        let rz = r.iter().all(|&x| x.is_zero());
        let lz = l.iter().all(|&x| x.is_zero());
        // step 1: T = S (x) E(r)
        if !rz {
            Self::powers(&self.off, n, r, &mut self.pw);
            for m in (1..=n).rev() {
                let rng = self.off[m - 1]..self.off[m];
                // k = 0 term: r^m
                for i in rng.clone() { self.t[i] = self.t[i].add(self.pw[i]); }
                for k in 1..m {
                    let (a, b) = (self.off[k - 1]..self.off[k], self.off[m - k - 1]..self.off[m - k]);
                    let coef = self.binom[m][k];
                    let (s_lo, s_hi) = self.t.split_at_mut(self.off[m - 1]);
                    outer_acc(coef, &s_lo[a], &self.pw[b], &mut s_hi[..rng.len()]);
                }
            }
        }
        // step 2: S' = E(l) (x) T
        if !lz {
            Self::powers(&self.off, n, l, &mut self.pl);
            for m in (1..=n).rev() {
                let rng = self.off[m - 1]..self.off[m];
                for i in rng.clone() { self.t[i] = self.t[i].add(self.pl[i]); }
                for j in 1..m {
                    let (a, b) = (self.off[j - 1]..self.off[j], self.off[m - j - 1]..self.off[m - j]);
                    let coef = self.binom[m][j];
                    let (s_lo, s_hi) = self.t.split_at_mut(self.off[m - 1]);
                    outer_acc(coef, &self.pl[a], &s_lo[b], &mut s_hi[..rng.len()]);
                }
            }
        }
    }
    /// Bits needed for level m given |x|<=xmax, window w:  log2(m!) + m*(log2(w)+log2(xmax)).  Must be < 127.
    pub fn bits_needed(n: usize, w: usize, xmax: f64) -> f64 {
        let mut worst = 0.0f64;
        for m in 1..=n { let f: f64 = (1..=m).map(|i| (i as f64).log2()).sum(); worst = worst.max(f + m as f64 * ((w as f64).log2() + xmax.log2())); }
        worst
    }
    
    /// Signature levels 1..=n in real units (increments / k_scale), flat.
    pub fn signature(&self, k_scale: f64) -> Vec<f64> {
        let mut out = vec![0.0; self.t.len()];
        for m in 1..=self.n {
            let f: f64 = (1..=m).map(|i| i as f64).product();
            let sc = 1.0 / (f * k_scale.powi(m as i32));
            for i in self.lvl(m) { out[i] = self.t[i].to_f64() * sc; }
        }
        out
    }
}

impl<T: Ring> SigRing<T> {
    /// Rigorous sufficient condition: |t_m| <= L1^m for every level m (L1 = path length of the window), so if L1^N < 2^(BITS-1)
    /// the ring state equals the true integer state and `signature()` is valid.  Conservative (never a false 'safe').
    pub fn guard_ok(&self) -> bool {
        let l = self.l1 as f64; let lim = 2f64.powi(T::BITS as i32 - 1);
        let mut p = 1.0f64; for _ in 0..self.n { p *= l; if p >= lim { return false; } } true
    }
}

// ---------------------------------------------------------------------------------------------
// Log-signature readout.  The ring state is exact; the log needs 1/k so it is taken in f64 *at readout only*,
// starting from the correctly-rounded exact signature (no accumulated drift).
// ---------------------------------------------------------------------------------------------
/// Truncated tensor product of two level-1..n tensors (no scalar part): out_m = sum_{i+j=m} a_i (x) b_j.
fn tmul_f64(off: &[usize], n: usize, a: &[f64], b: &[f64], out: &mut [f64]) {
    for o in out.iter_mut() { *o = 0.0; }
    for m in 2..=n {
        for i in 1..m {
            let j = m - i;
            let (ra, rb) = (off[i - 1]..off[i], off[j - 1]..off[j]);
            let lb = rb.len();
            let base = off[m - 1];
            for (ia, &x) in a[ra].iter().enumerate() {
                let row = &mut out[base + ia * lb..base + (ia + 1) * lb];
                for (o, &y) in row.iter_mut().zip(b[rb.clone()].iter()) { *o += x * y; }
            }
        }
    }
}
/// Expanded tensor logarithm of a signature given as levels 1..n (scalar 1 implied). Returns levels 1..n.
pub fn tensor_log(d: usize, n: usize, sig: &[f64]) -> Vec<f64> {
    let mut off = vec![0usize; n + 1]; let mut sz = 1usize; let mut acc = 0usize;
    for k in 1..=n { sz *= d; acc += sz; off[k] = acc; }
    let mut res = sig.to_vec();
    let mut p = sig.to_vec(); let mut tmp = vec![0.0; sig.len()];
    for k in 2..=n {
        tmul_f64(&off, n, &p, sig, &mut tmp);
        std::mem::swap(&mut p, &mut tmp);
        let c = if k % 2 == 0 { -1.0 } else { 1.0 } / k as f64;
        for (r, &v) in res.iter_mut().zip(p.iter()) { *r += c * v; }
    }
    res
}
/// Lyndon words of length 1..=n over alphabet 0..d, ordered by (length, lexicographic).
pub fn lyndon_words(d: usize, n: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut w: Vec<usize> = vec![0];
    loop {
        if w.len() <= n { out.push(w.clone()); }
        // Duval next
        let m = w.len();
        while w.len() < n { let i = w.len() - m; let x = w[i]; w.push(x); }
        while let Some(&last) = w.last() { if last == d - 1 { w.pop(); } else { break; } }
        match w.last_mut() { Some(l) => *l += 1, None => break }
    }
    out.sort_by(|a, b| a.len().cmp(&b.len()).then(a.cmp(b)));
    out.dedup();
    out
}
/// Precomputed projection of the expanded log onto the Lyndon bracket basis (coordinates) and onto Lyndon-word coefficients.
pub struct LogSigPrep {
    pub d: usize, pub n: usize,
    pub words: Vec<Vec<usize>>,
    off: Vec<usize>,
    idx: Vec<usize>,            // flat index (in levels 1..n layout) of each Lyndon word
    brackets: Vec<Vec<f64>>,    // P_w expanded over words of its own length (dense, d^len)
}
fn word_index(d: usize, w: &[usize]) -> usize { w.iter().fold(0usize, |a, &x| a * d + x) }
impl LogSigPrep {
    pub fn new(d: usize, n: usize) -> Self {
        let words = lyndon_words(d, n);
        let mut off = vec![0usize; n + 1]; let mut sz = 1usize; let mut acc = 0usize;
        for k in 1..=n { sz *= d; acc += sz; off[k] = acc; }
        let idx: Vec<usize> = words.iter().map(|w| off[w.len() - 1] + word_index(d, w)).collect();
        // standard factorisation: w = u v with v the longest proper Lyndon suffix
        use std::collections::HashMap;
        let mut pos: HashMap<Vec<usize>, usize> = HashMap::new();
        for (i, w) in words.iter().enumerate() { pos.insert(w.clone(), i); }
        let mut brackets: Vec<Vec<f64>> = vec![Vec::new(); words.len()];
        for (i, w) in words.iter().enumerate() {
            let len = w.len();
            let mut dense = vec![0.0; d.pow(len as u32)];
            if len == 1 { dense[w[0]] = 1.0; } else {
                let mut split = 0;
                for s in 1..len { if pos.contains_key(&w[s..].to_vec()) { split = s; break; } }
                let (u, v) = (&words[pos[&w[..split].to_vec()]], &words[pos[&w[split..].to_vec()]]);
                let (pu, pv) = (&brackets[pos[&u.clone()]], &brackets[pos[&v.clone()]]);
                let (lu, lv) = (u.len(), v.len());
                let (nu, nv) = (d.pow(lu as u32), d.pow(lv as u32));
                for a in 0..nu { if pu[a] == 0.0 { continue; } for b in 0..nv {
                    dense[a * nv + b] += pu[a] * pv[b]; dense[b * nu + a] -= pv[b] * pu[a]; } }
            }
            brackets[i] = dense;
        }
        Self { d, n, words, off, idx, brackets }
    }
    pub fn len(&self) -> usize { self.words.len() }
    /// Lyndon-word coefficients of an expanded log (method "words").
    pub fn lyndon_word_coefs(&self, expanded: &[f64]) -> Vec<f64> { self.idx.iter().map(|&i| expanded[i]).collect() }
    /// Coordinates of the expanded log (a Lie element) in the Lyndon bracket basis.
    pub fn bracket_coords(&self, expanded: &[f64]) -> Vec<f64> {
        let mut resid = expanded.to_vec();
        let mut c = vec![0.0; self.words.len()];
        for (i, w) in self.words.iter().enumerate() {
            let k = resid[self.idx[i]]; c[i] = k;
            if k != 0.0 { let base = self.off[w.len() - 1]; for (j, &p) in self.brackets[i].iter().enumerate() { resid[base + j] -= k * p; } }
        }
        c
    }
}
impl<T: Ring> SigRing<T> {
    /// Expanded log-signature (levels 1..n, tensor basis) in real units.
    pub fn log_signature_expanded(&self, k_scale: f64) -> Vec<f64> { tensor_log(self.d, self.n, &self.signature(k_scale)) }
}

// ---------------------------------------------------------------------------------------------
// Division-free Horner update.  For target level m:  t'_m = M_m,  M_0 = 1,  M_j = M_{j-1} (x) r + C(m,j) t_j   (right push)
//                                                    M_j = l (x) M_{j-1} + C(m,j) t_j                              (left pop)
// Same exact result in the wrapping ring (ring arithmetic is associative/commutative) with ~1.6x fewer multiplications at d=2,N=4
// and no power tensors.  Process m descending so lower levels are still the old ones.
// ---------------------------------------------------------------------------------------------
pub struct HornerScratch<T: Ring> { a: Vec<T>, b: Vec<T> }
impl<T: Ring> HornerScratch<T> { pub fn new(d: usize, n: usize) -> Self { let sz = d.pow(n as u32).max(1); Self { a: vec![T::default(); sz], b: vec![T::default(); sz] } } }
impl<T: Ring> SigRing<T> {
    pub fn ingest_horner(&mut self, l: &[T], r: &[T], sc: &mut HornerScratch<T>) {
        let (d, n) = (self.d, self.n);
        if !r.iter().all(|&x| x.is_zero()) {
            for m in (1..=n).rev() {
                // M_1
                let t1 = self.off[0]..self.off[1];
                let c1 = self.binom[m][1];
                for b in 0..d { sc.a[b] = r[b].add(c1.mul(self.t[t1.start + b])); }
                let mut len = d;
                for j in 2..=m {
                    let tj = self.off[j - 1];
                    let cj = if j == m { T::from_i64(1) } else { self.binom[m][j] };
                    let nl = len * d;
                    for a in 0..len { let x = sc.a[a]; let row = &mut sc.b[a * d..(a + 1) * d]; for b in 0..d { row[b] = x.mul(r[b]).add(cj.mul(self.t[tj + a * d + b])); } }
                    std::mem::swap(&mut sc.a, &mut sc.b); len = nl;
                }
                let lo = self.off[m - 1];
                self.t[lo..lo + len].copy_from_slice(&sc.a[..len]);
            }
        }
        if !l.iter().all(|&x| x.is_zero()) {
            for m in (1..=n).rev() {
                let c1 = self.binom[m][1];
                for b in 0..d { sc.a[b] = l[b].add(c1.mul(self.t[self.off[0] + b])); }
                let mut len = d;
                for j in 2..=m {
                    let tj = self.off[j - 1];
                    let cj = if j == m { T::from_i64(1) } else { self.binom[m][j] };
                    let nl = len * d;
                    for b in 0..d { let lb = l[b]; for a in 0..len { sc.b[b * len + a] = lb.mul(sc.a[a]).add(cj.mul(self.t[tj + b * len + a])); } }
                    std::mem::swap(&mut sc.a, &mut sc.b); len = nl;
                }
                let lo = self.off[m - 1];
                self.t[lo..lo + len].copy_from_slice(&sc.a[..len]);
            }
        }
    }
}

/// Benchmark/baseline only: lets the same Horner kernel run in f64 so float baselines are optimised identically.
#[doc(hidden)]
impl Ring for f64 {
    #[inline(always)] fn add(self, o: Self) -> Self { self + o }
    #[inline(always)] fn mul(self, o: Self) -> Self { self * o }
    #[inline(always)] fn neg(self) -> Self { -self }
    #[inline(always)] fn from_i64(v: i64) -> Self { v as f64 }
    fn to_f64(self) -> f64 { self }
    #[inline(always)] fn is_zero(self) -> bool { self == 0.0 }
    const BITS: u32 = 53;
}

// ---------------------------------------------------------------------------------------------
// Chen product and antipode in the scaled domain, and exact arbitrary-interval queries.
// Product (levels 1..n, scalar 1 implied):  c_m = a_m + b_m + sum_{k=1}^{m-1} C(m,k) a_k (x) b_{m-k}.
// Antipode (= group inverse of a signature): level m, word w -> (-1)^m * word reversed.
// Because the ring is a ring, antipode(prefix_i) (x) prefix_j equals the true state of the increments i..j (mod 2^BITS)
// even when the prefixes themselves overflowed; it is exact iff the interval's own state fits (see `interval_guard_ok`).
// ---------------------------------------------------------------------------------------------
fn level_offsets(d: usize, n: usize) -> Vec<usize> {
    let mut off = vec![0usize; n + 1]; let mut sz = 1usize;
    for k in 1..=n { sz *= d; off[k] = off[k - 1] + sz; }
    off
}
fn binom_u(n: usize, k: usize) -> i64 { let mut r: i128 = 1; for i in 0..k { r = r * (n - i) as i128 / (i + 1) as i128; } r as i64 }

/// Concatenation (Chen product) of two scaled states.
pub fn state_concat<T: Ring>(d: usize, n: usize, a: &[T], b: &[T]) -> Vec<T> {
    let off = level_offsets(d, n);
    let mut c: Vec<T> = (0..a.len()).map(|i| a[i].add(b[i])).collect();
    for m in 2..=n {
        for k in 1..m {
            let coef = T::from_i64(binom_u(m, k));
            let (ra, rb) = (off[k - 1]..off[k], off[m - k - 1]..off[m - k]);
            let lb = rb.len();
            let base = off[m - 1];
            for (ia, &x) in a[ra].iter().enumerate() {
                let cx = coef.mul(x);
                for (ib, &y) in b[rb.clone()].iter().enumerate() { let o = base + ia * lb + ib; c[o] = c[o].add(cx.mul(y)); }
            }
        }
    }
    c
}

/// Antipode (group inverse) of a scaled state.
pub fn state_antipode<T: Ring>(d: usize, n: usize, a: &[T]) -> Vec<T> {
    let off = level_offsets(d, n);
    let mut out = vec![T::default(); a.len()];
    for m in 1..=n {
        let base = off[m - 1];
        for idx in 0..d.pow(m as u32) {
            let (mut rem, mut rev) = (idx, 0usize);
            for _ in 0..m { rev = rev * d + rem % d; rem /= d; }
            let v = a[base + idx];
            out[base + rev] = if m % 2 == 1 { v.neg() } else { v };
        }
    }
    out
}

/// Stores the exact prefix states of a stream so that the signature of any increment interval `[i, j)` can be queried in
/// O(N d^N-ish) with no recomputation and no float inverse.
pub struct PrefixStore<T: Ring = i128> {
    pub d: usize, pub n: usize,
    eng: SigRing<T>,
    prefixes: Vec<Vec<T>>,   // prefixes[i] = state of increments 0..i
    l1cum: Vec<i128>,        // cumulative L1 of increments
}
impl<T: Ring> PrefixStore<T> {
    pub fn new(d: usize, n: usize) -> Self {
        let eng = SigRing::new(d, n);
        Self { d, n, prefixes: vec![eng.t.clone()], l1cum: vec![0], eng }
    }
    pub fn push(&mut self, inc: &[i64]) {
        let zero = vec![0i64; self.d];
        self.eng.ingest(inc, &zero);
        self.prefixes.push(self.eng.t.clone());
        let l: i128 = inc.iter().map(|&x| (x as i128).abs()).sum();
        self.l1cum.push(self.l1cum.last().unwrap() + l);
    }
    pub fn len(&self) -> usize { self.prefixes.len() - 1 }
    /// Scaled state of increments `i..j`.
    pub fn interval_state(&self, i: usize, j: usize) -> Vec<T> {
        state_concat(self.d, self.n, &state_antipode(self.d, self.n, &self.prefixes[i]), &self.prefixes[j])
    }
    /// Sufficient condition for the interval state to be exact (L1 of the interval to the N-th power < 2^(BITS-1)).
    pub fn interval_guard_ok(&self, i: usize, j: usize) -> bool {
        let l = (self.l1cum[j] - self.l1cum[i]) as f64; let lim = 2f64.powi(T::BITS as i32 - 1);
        let mut p = 1.0f64; for _ in 0..self.n { p *= l; if p >= lim { return false; } } true
    }
    /// Signature levels 1..=n in real units for increments `i..j`.
    pub fn interval_signature(&self, i: usize, j: usize, k_scale: f64) -> Vec<f64> {
        let mut e: SigRing<T> = SigRing::new(self.d, self.n);
        e.t = self.interval_state(i, j);
        e.signature(k_scale)
    }
}
