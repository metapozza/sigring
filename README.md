# sigring

Exact, streaming, **sliding-window truncated path signatures** in a wrapping integer ring — Rust crate (plus a PyO3 wheel, see `python/`).

Push one increment and drop the oldest one in a single fused update whose cost does **not** depend on the window length. The state is integer, so the result is exact (when it fits), bit-reproducible on every platform, and you can *prove* at runtime that it fits.

```rust
use sigring::SigRing;
let (d, n, window) = (2, 4, 64);
let mut eng: SigRing<i64> = SigRing::new(d, n);      // Z/2^64; use i128 for headroom
// increments are integers: quantise real data with a scale K (e.g. price * 1e4)
eng.ingest(&new_inc, &oldest_inc);                    // push + pop, no allocation
assert!(eng.guard_ok());                              // rigorous overflow certificate
let sig = eng.signature(1e4);                         // levels 1..=n in real units
```

## Install

```bash
cargo add sigring            # Rust (crates.io)
pip install sigring numpy    # Python wheel (CPython >= 3.9, Linux x86_64 wheel; other platforms build from source and need a Rust toolchain)
```

## What it is (and is not)

* **Is:** an exact integer formulation of the sliding signature update `S' = E(-old) ⊗ S ⊗ E(new)` on the factorial-cleared state `t_n = n!·S_n·K^n`, using a division-free Horner kernel; a rigorous guard; exact arbitrary-interval queries (`PrefixStore`, antipode instead of a float inverse); log-signatures (Lyndon words / Lyndon-bracket coordinates, log taken in f64 at readout).
* **Is not:** a new signature algorithm. The Horner "fused multiply-exponentiate" form is Signatory's (Kidger & Lyons 2020); sliding-window aggregation in general is a known problem (Tangwongsan et al., SWAG). What is specific here is the exact integer ring formulation, the guard theorem, determinism, and no per-window recomputation. See `docs/LITERATURE.md`.
* **Is not** a differentiable / GPU / batch-training library. Use `pySigLib`, `signax`, `iisignature` or `pathsig` for that.

## Guard theorem

Every entry of the scaled state satisfies `|t_m| ≤ L1^m`, where `L1` is the sliding sum of the L1 norms of the integer increments (ℓ1 tensor norm is multiplicative and `‖S_m‖₁ ≤ L1^m / m!`). So the state is guaranteed exact when `L1^N < 2^(BITS-1)`. The bound is tight (attained by 1-D monotone paths). Verified exactly on 32,400 (window, level) pairs with 0 violations; `guard_ok()` is conservative (never a false "safe"), but it can say "unsafe" for states that happen to fit. The guard depends on your **units**: quantise coarsely (or use `i128`) if depth × resolution is large. Real data: the `i128` guard certified 100% of the windows on every real dataset we tried; the `i64` guard passes at coarse units and fails at depth 4 with cent / 1e-5 units.

## Measured (read the caveats)

All numbers from one 2-core 2.1 GHz Xeon sandbox, ±25% run-to-run noise. Reproduce with `cd bench/sigbench && cargo build --release && ./target/release/sigbench bench <int64|int|fused_float_h|twostack_h|recompute> <d> <N> <W> <ticks> <kind>` (kind 0 = random walk), `cargo run --release --example bench_kernels` (naive vs Horner kernel), and the Python scripts in `python/bench/` (need `pysiglib`/`iisignature`). Real-market-data inputs and the pipeline-shape script are not included in this repo; those numbers are reported as measured, not as reproducible from here.

| (d, N, W) | i64 ring | i128 ring | f64 fused (inverse, Horner) | f64 two-stack (Horner) |
|---|---|---|---|---|
| 2, 4, 64 | 329 ns | 442 ns | 233 ns | 427 ns |
| 3, 4, 64 | 511 ns | 1151 ns | 436 ns | 743 ns |
| 5, 3, 64 | 382 ns | 958 ns | 386 ns | 651 ns |
| 2, 4, 256 | 220 ns | 386 ns | 175 ns | 394 ns |

* **Be honest about speed:** the i64 ring is *roughly tied* with the best float fused update (which uses an inverse and can be badly wrong, below) and ~1.2–2× faster than a float two-stack at small shapes (2.9× at d=15). It is **not** faster than float in general.
* Versus recomputing every window with a fast library (d=2, N=4): ≈4–5× over pySigLib at W=64, ≈12× at W=256, ≈47× at W=1024, parity at W=16.
* Why bother with exactness: the float fused update with an inverse had max relative error 1.7e-11 on random walks, 8.1e-6 on heavy tails and 0.56 on spike streams (d=2, N=4, W=64); on real ES (price, time) ticks across a 6.5-hour gap it was off by 4.3e9. Float two-stack/recompute stay ≈1e-15 but cost more memory (≈9.7 MB vs ≈69 KB at d=15, N=3, W=336) and have O(W) worst-case flips (per-push p99.9 4.6 ms vs 55 µs there).
* Log-signature (ingest + Lyndon-bracket readout): ≈0.7–0.9 µs/tick vs pySigLib recompute 1.1 / 3.3 / 15.6 µs and iisignature 7 / 30 / 139 µs at W = 64 / 256 / 1024. Matches pySigLib (Lyndon words and brackets) to 6e-16 and `signature-rs` to 1.6e-15.
* Pipeline shape (real EURUSD hourly, 4,832 windows, d=15, N=3, 336 increments): i64 ring 39 µs/window, i128 71 µs, pySigLib 243 µs, iisignature 3,326 µs. Agreement with iisignature: 1.5e-16 at N=3.
* Python wheel: `push()` ≈250–320 ns, `ingest_batch` ≈180–225 ns/tick; pySigLib's `SigWindowStream` costs 18–32 µs/tick from Python.

## Determinism

State hashes are identical across 8 Rust build variants and an independent C port on x86 (gcc/clang), aarch64, riscv64, ppc64 (BE), s390x (BE), arm32 and mips32 (BE, i64 only), and wasm32 (node). This holds by construction for integer arithmetic. Note a plain-f64 *Rust* build was also bit-identical across our 6 x86 variants; the practical advantage is for C/C++ libraries (FMA contraction split a C float control into two bit patterns), threaded/GPU reductions and cross-architecture reproducibility.

Reproducing the hashes: `det/ref.txt` holds the reference hashes. `cargo run --release --example det_hash -- 2 4 64 200000 50000 64` and the independent C port (`gcc -O2 -DWORD=64 det/det.c -o detc && ./detc 2 4 64 200000 50000`) both print the `2 4 64 64 ...` lines of `det/ref.txt`. The cross-architecture runs used qemu-user and cross GCC toolchains (not scripted here); CI re-runs the hash on Linux, macOS (arm64) and Windows runners.

## Limitations

* Increments must be integers after scaling; `K` and units are the user's responsibility. Depth × bit-width bounds the usable resolution (the guard tells you).
* Cost per tick is O(d^N · N) like any dense truncated signature; memory is the signature size. No sparse/Lyndon-basis internal state yet.
* Single-threaded, CPU only. No autodiff.
* Not tested: real embedded targets, a real zkVM, DABA-style baselines, live market feeds. Benchmarks are one noisy machine.
* Competitors that could not be run: `signatory` (does not build against current torch), `pathsig` (needs CUDA).

## Tests

`cargo test --release`: sliding == recompute (i64, i128), Horner == naive (bit-identical), exact interval queries (including when the prefixes themselves overflowed `i64`), antipode is the group inverse, `|t_m| ≤ L1^m`, log-signature of a single segment, Witt-formula Lyndon counts.

## License

Source-available under the [PolyForm Noncommercial License 1.0.0](LICENSE) — free for personal, research, educational and non-profit use; commercial use needs a licence. See [COMMERCIAL.md](COMMERCIAL.md) (includes a free 30-day evaluation permission). This is not an OSI-approved open-source licence.
