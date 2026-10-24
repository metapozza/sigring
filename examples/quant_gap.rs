//! Synthetic tick stream with an overnight gap: exact integer ring vs. the float fused update (which uses an inverse).
//!
//! Path = (price in ticks, elapsed milliseconds). Every GAP_EVERY ticks the clock jumps by 6.5 hours (23,400,000 ms).
//! The float state absorbs the huge increment and, after it leaves the window, keeps a rounding residue that is large
//! relative to the true (small) values. The exact ring does not.
//!
//! Run:  cargo run --release --example quant_gap
//! Synthetic data only: it illustrates the mechanism, it is not a measurement on market data.
use sigring::SigRing;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); self.0 >> 33 }
}

fn main() {
    const D: usize = 2; const N: usize = 4; const W: usize = 64;
    const TICKS: usize = 200_000; const GAP_EVERY: usize = 20_000; const GAP_MS: i64 = 23_400_000;
    let mut rng = Lcg(42);
    let mut incs: Vec<[i64; D]> = Vec::with_capacity(TICKS);
    for t in 0..TICKS {
        let dp = (rng.next() % 7) as i64 - 3;                 // price change, in ticks
        let dt = if t > 0 && t % GAP_EVERY == 0 { GAP_MS } else { 50 + (rng.next() % 100) as i64 };
        incs.push([dp, dt]);
    }
    let mut exact: SigRing<i128> = SigRing::new(D, N);
    let mut exact64: SigRing<i64> = SigRing::new(D, N);
    let mut fl: SigRing<f64> = SigRing::new(D, N);
    let zero = [0i64; D];
    let (mut worst_after_gap, mut worst_all) = (0.0f64, 0.0f64);
    let (mut windows, mut guard128_ok, mut guard64_ok) = (0usize, 0usize, 0usize);
    for t in 0..TICKS {
        let out = if t >= W { incs[t - W] } else { zero };
        exact.ingest(&incs[t], &out); exact64.ingest(&incs[t], &out); fl.ingest(&incs[t], &out);
        if t + 1 < W { continue; }
        windows += 1;
        if exact.guard_ok() { guard128_ok += 1; }
        if exact64.guard_ok() { guard64_ok += 1; }
        // compare only windows that contain no gap increment but follow one (the float residue case)
        let has_gap = incs[t + 1 - W..=t].iter().any(|x| x[1] == GAP_MS);
        let (se, sf) = (exact.signature(1.0), fl.signature(1.0));
        let scale = se.iter().fold(0.0f64, |a, &v| a.max(v.abs())).max(1e-300);
        let err = se.iter().zip(&sf).fold(0.0f64, |a, (&x, &y)| a.max((x - y).abs())) / scale;
        worst_all = worst_all.max(err);
        if !has_gap && t > GAP_EVERY { worst_after_gap = worst_after_gap.max(err); }
    }
    println!("d={D} N={N} window={W} ticks={TICKS} gap={GAP_MS} ms every {GAP_EVERY} ticks ({windows} full windows)");
    println!("exact i128 ring: guard certified {guard128_ok}/{windows} windows");
    println!("exact i64  ring: guard certified {guard64_ok}/{windows} windows (the guard refuses windows that contain a gap)");
    println!("float fused update, worst relative error vs exact: {worst_all:.3e} overall, {worst_after_gap:.3e} in gap-free windows after a gap");
}
