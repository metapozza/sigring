use sigring::SigRing;
use std::time::Instant;
fn go<T: sigring::Ring + std::fmt::Debug>(name: &str, d: usize, n: usize, w: usize, ticks: usize) {
    let mut incs = vec![0i64; ticks * d]; let mut s: u64 = 77;
    for x in incs.iter_mut() { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); *x = ((s >> 33) % 2001) as i64 - 1000; }
    // equality
    let mut g: SigRing<T> = SigRing::new(d, n); let mut h: SigRing<T> = SigRing::new(d, n);
    let (mut l, mut r) = (vec![T::default(); d], vec![T::default(); d]);
    for i in 0..30000 { let inc = &incs[i * d..(i + 1) * d]; let old: Vec<i64> = if i >= w { incs[(i - w) * d..(i - w + 1) * d].to_vec() } else { vec![0; d] };
        g.ingest_naive(inc, &old);
        for c in 0..d { r[c] = T::from_i64(inc[c]); l[c] = T::from_i64(old[c]).neg(); } h.ingest(inc, &old); }
    assert!(g.t == h.t, "MISMATCH {} d={} n={}", name, d, n);
    let mut best = [f64::MAX; 2];
    for _ in 0..5 {
        let mut g: SigRing<T> = SigRing::new(d, n); let t0 = Instant::now();
        let mut old = vec![0i64; d];
        for i in 0..ticks { let inc = &incs[i * d..(i + 1) * d]; if i >= w { old.copy_from_slice(&incs[(i - w) * d..(i - w + 1) * d]); } g.ingest_naive(inc, &old); }
        best[0] = best[0].min(t0.elapsed().as_secs_f64() * 1e9 / ticks as f64);
        let mut h: SigRing<T> = SigRing::new(d, n); let t0 = Instant::now(); let mut old = vec![0i64; d];
        for i in 0..ticks { let inc = &incs[i * d..(i + 1) * d]; if i >= w { old.copy_from_slice(&incs[(i - w) * d..(i - w + 1) * d]); } h.ingest(inc, &old); }
        best[1] = best[1].min(t0.elapsed().as_secs_f64() * 1e9 / ticks as f64);
        if g.t[0] == T::from_i64(1) && h.t[0] == T::from_i64(1) { println!("!"); }
    }
    println!("{:>5} d={} N={} W={}: naive {:6.0} ns | horner {:6.0} ns | speedup {:.2}x  (bit-identical state)", name, d, n, w, best[0], best[1], best[0] / best[1]);
}
fn main() {
    for &(d, n) in &[(2usize, 3usize), (2, 4), (3, 3), (3, 4), (5, 3), (2, 5), (8, 2), (4, 4)] { go::<i64>("i64", d, n, 64, 400_000); }
    for &(d, n) in &[(2usize, 4usize), (3, 3), (5, 3), (2, 5), (3, 4)] { go::<i128>("i128", d, n, 64, 300_000); }
}
