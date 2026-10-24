First public release.

**What it is:** exact, streaming, sliding-window truncated path signatures in a wrapping integer ring (i64 / i128). A window slide is one fused push+pop with cost independent of the window length; results are bit-identical across platforms; a rigorous guard (`|t_m| <= L1^m`) certifies when the state is exact. Also: exact arbitrary-interval queries (`PrefixStore`, antipode instead of a float inverse) and log-signatures (Lyndon words / bracket coordinates).

**Verified:** sliding == recompute (i64, i128), Horner == naive (bit-identical), interval queries exact even when prefixes overflowed, guard bound holds, state hashes identical across x86, aarch64, riscv64, ppc64 BE, s390x BE, arm32, mips32 BE, wasm32.

**Honest limits:** speed is roughly equal to the best float fused update (not faster); it is faster than recomputing each window (about 4-5x at W=64, 12x at W=256, 47x at W=1024 vs pySigLib at d=2, N=4). The Horner kernel is not novel (Signatory). Integer increments only; units and depth bound usable resolution (the guard tells you). Single-threaded CPU. Benchmarks are from one noisy 2-core machine.

**Licence:** source-available, PolyForm Noncommercial 1.0.0 (free for personal, research, educational and non-profit use; 30-day evaluation permission and commercial licences: see COMMERCIAL.md). Not OSI open source.

Crate: https://crates.io/crates/sigring
