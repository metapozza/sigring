use sigring::*;

struct Rng(u64);
impl Rng { fn next(&mut self, amp: i64) -> i64 { self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((self.0 >> 33) % (2 * amp as u64 + 1)) as i64 - amp } }

fn recompute<T: Ring>(d: usize, n: usize, incs: &[Vec<i64>]) -> Vec<T> {
    let mut e: SigRing<T> = SigRing::new(d, n); let z = vec![0i64; d];
    for x in incs { e.ingest(x, &z); } e.t
}
fn stream(d: usize, len: usize, amp: i64, off: i64, seed: u64) -> Vec<Vec<i64>> {
    let mut r = Rng(seed); (0..len).map(|_| (0..d).map(|_| r.next(amp) + off).collect()).collect()
}

fn window_matches<T: Ring + std::fmt::Debug>(d: usize, n: usize, w: usize) {
    let incs = stream(d, 500, 40, 0, 7);
    let mut e: SigRing<T> = SigRing::new(d, n); let z = vec![0i64; d];
    for i in 0..incs.len() {
        let out = if i >= w { incs[i - w].clone() } else { z.clone() };
        e.ingest(&incs[i], &out);
        if i % 37 == 0 { let lo = (i + 1).saturating_sub(w); assert_eq!(e.t, recompute::<T>(d, n, &incs[lo..=i]), "d={d} n={n} i={i}"); }
    }
}
#[test] fn sliding_equals_recompute_i64() { window_matches::<i64>(2, 4, 32); window_matches::<i64>(3, 3, 16); }
#[test] fn sliding_equals_recompute_i128() { window_matches::<i128>(2, 4, 32); window_matches::<i128>(3, 4, 16); window_matches::<i128>(5, 3, 20); }

#[test] fn horner_equals_naive() {
    let (d, n, w) = (3, 4, 24); let incs = stream(d, 400, 500, 0, 3);
    let (mut a, mut b): (SigRing<i128>, SigRing<i128>) = (SigRing::new(d, n), SigRing::new(d, n)); let z = vec![0i64; d];
    for i in 0..incs.len() { let out = if i >= w { &incs[i - w] } else { &z }; a.ingest(&incs[i], out); b.ingest_naive(&incs[i], out); }
    assert_eq!(a.t, b.t);
}

/// Prefixes overflow i64 (offset 1e3 per step, 3000 steps, depth 4) but a short interval fits: the query must still be exact.
#[test] fn interval_exact_even_when_prefixes_overflow() {
    let (d, n) = (2, 4); let incs = stream(d, 3000, 50, 1_000, 11);
    let mut p: PrefixStore<i64> = PrefixStore::new(d, n);
    for x in &incs { p.push(x); }
    let mut checked = 0;
    for &(i, j) in &[(0usize, 8usize), (500, 509), (1000, 1007), (2990, 3000), (3, 3)] {
        if !p.interval_guard_ok(i, j) { continue; }
        assert_eq!(p.interval_state(i, j), recompute::<i64>(d, n, &incs[i..j]), "i={i} j={j}"); checked += 1;
    }
    assert!(checked >= 3, "guard rejected too many intervals ({checked})");
    // and the whole prefix state really has wrapped (guard says unsafe for the full range)
    assert!(!p.interval_guard_ok(0, 3000));
}
#[test] fn interval_random_i128() {
    let (d, n) = (3, 4); let incs = stream(d, 120, 30, 5, 5);
    let mut p: PrefixStore<i128> = PrefixStore::new(d, n); for x in &incs { p.push(x); }
    let mut r = Rng(99);
    for _ in 0..200 { let i = r.next(59).abs() as usize; let j = i + r.next(30).unsigned_abs() as usize % (120 - i + 1);
        assert_eq!(p.interval_state(i, j), recompute::<i128>(d, n, &incs[i..j])); }
}
#[test] fn antipode_is_inverse() {
    let (d, n) = (2, 4); let incs = stream(d, 50, 20, 0, 1);
    let mut e: SigRing<i128> = SigRing::new(d, n); let z = vec![0i64; d]; for x in &incs { e.ingest(x, &z); }
    let id = state_concat(d, n, &e.t, &state_antipode(d, n, &e.t));
    assert!(id.iter().all(|&v| v == 0));
}
/// |t_m| <= L1^m for every level (exact-integer check of the guard bound; i128, no overflow in this regime).
#[test] fn guard_bound_holds() {
    let (d, n, w) = (3, 4, 12); let incs = stream(d, 300, 25, 3, 21);
    let mut e: SigRing<i128> = SigRing::new(d, n); let z = vec![0i64; d];
    let off: Vec<usize> = { let mut o = vec![0]; let mut s = 1; for _ in 1..=n { s *= d; o.push(o.last().unwrap() + s); } o };
    for i in 0..incs.len() {
        let out = if i >= w { &incs[i - w] } else { &z }; e.ingest(&incs[i], out);
        let l1 = e.l1; assert!(e.guard_ok());
        for m in 1..=n { for k in off[m - 1]..off[m] { assert!(e.t[k].unsigned_abs() <= (l1 as u128).pow(m as u32)); } }
    }
}
/// A single segment has log-signature equal to the segment.
#[test] fn logsig_single_segment() {
    let (d, n) = (3, 4); let inc = [7i64, -3, 5];
    let mut e: SigRing<i128> = SigRing::new(d, n); e.ingest(&inc, &[0, 0, 0]);
    let prep = LogSigPrep::new(d, n); let b = prep.bracket_coords(&e.log_signature_expanded(1.0));
    assert_eq!(b.len(), prep.len());
    for (i, w) in prep.words.iter().enumerate() {
        let want = if w.len() == 1 { inc[w[0]] as f64 } else { 0.0 };
        assert!((b[i] - want).abs() < 1e-9, "word {w:?}: {} vs {}", b[i], want);
    }
}
#[test] fn lyndon_counts_match_witt() {
    // number of Lyndon words of length k over d letters = (1/k) sum_{e|k} mu(k/e) d^e ; totals for d=2: 2,1,2,3,6
    let l = lyndon_words(2, 5); let mut c = [0usize; 6]; for w in &l { c[w.len()] += 1; }
    assert_eq!(&c[1..], &[2, 1, 2, 3, 6]);
}
