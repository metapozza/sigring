# Prior art and positioning

Scope: what was checked, what was read in full, and what was only searched. Claims are limited to what was verified in this project; "not verified" is stated explicitly.

## Verified by reading the source / paper

| Work | Relevant fact | Consequence for sigring |
|---|---|---|
| Kidger & Lyons, *Signatory* (2020), §4.1 / App. A.1 | Fused multiply-exponentiate computed Horner-style, with divisions. | The Horner kernel is **not novel**. sigring's variant is division-free on the factorial-scaled integer state. |
| Signatory §4.2 / §5.5 | Arbitrary-interval queries via stored prefix signatures and a (float) inverse; the paper warns of instability for large intervals. | sigring's `PrefixStore` does the same in the exact ring using the antipode (index reversal with sign (−1)^m). We did **not** reproduce a large float instability in our small test (rel. err 1.3e-14 at offset 1e7, level 4), so the exactness claim is about guarantees, not demonstrated float failure there. |
| pySigLib (Apache-2.0) paper and code | Horner-based batch signatures; `SigWindowStream` recomputes per window; no streaming sliding algorithm described. | Main practical comparator for per-window cost (see README table). |
| Sliding-window aggregation (SWAG, arXiv 2009.13768) | Two-stacks: amortised O(1), worst-case O(n) flip, 2n aggregates; DABA: worst-case O(1), 2n aggregates; DABA Lite: n+2 aggregates. Need only associativity. | The signature is a monoid, so the two-stack baseline is the fair float comparator and is what the benchmarks use. DABA itself was **not** benchmarked. |
| pathsig (arXiv 2602.24066) | CUDA-only; states that the Chen inverse for windows is numerically unstable. | Not runnable here (no GPU). |
| iisignature, esig, RoughPy, signax | Run and compared (see README). | |
| signature-rs (AGPL-3.0) | A Rust crate exists: batch BCH log-signatures. Matches within 1.6e-15; 100–260 µs/window. | Earlier claim "no Rust signature crate exists" was wrong and is corrected. License is AGPL, so sigring does not derive from it. |
| sktime 1.2.0 / aeon 1.6.0 | Windowed signature features recompute every window via esig. | Concrete consumers that would benefit from an incremental engine. |

## Not verified

* Whether an exact-integer sliding signature has been published elsewhere. Searches of crates.io, PyPI, GitHub and arXiv found none, which is absence of evidence, not proof. Claim in the README is therefore "specific to this crate" only in the sense of "not found".
* Rolling-hash precedent (Rabin–Karp mod 2^64) and cumulative on-chain accumulators (Uniswap v2 price accumulators) use the same modular-subtraction principle for scalars; only search results were seen, not read in full.
* DABA / DABA Lite performance on signatures; embedded or zkVM targets; live market feeds.

## Academic value (assessment)

Not a new theorem-level result. Defensible contributions: (1) the guard bound `|t_m| ≤ L1^m` with tightness and its use as a runtime certificate in a wrapping ring, (2) the exact integer sliding/interval formulation, (3) an empirical determinism and real-data study. Suitable for a short note or workshop paper at best; the primary value is engineering.
