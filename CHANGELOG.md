# Changelog

## 0.1.2
- Docs: new "For quant and systematic-trading use" README section (units, guard sizing table, checklist), contact block, absolute links so the crates.io page renders correctly.
- Examples: `examples/quant_gap.rs` (float vs exact on a stream with clock gaps), `python/examples/quant_ticks.py` (per-window guard, i64 vs i128).
- Project hygiene: `SECURITY.md`, `CONTRIBUTING.md`, issue templates, Python and examples jobs in CI.
- No change to the library code or results.

## 0.1.1
- README: install section, reproducible benchmark and determinism instructions (`bench/sigbench`, `det/`).
- Python: `WindowSig` reports module `sigring`; removed an unused import; wheel depends on `sigring = "0.1"`.
- Package no longer ships `det/`, `bench/` or CI files.
- No functional change to the Rust API or results.

## 0.1.0
First public release.
