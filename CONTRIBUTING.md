# Contributing

Bug reports, questions, benchmarks that disagree with the README, and counter-examples to the guard are the most valuable contributions. Open an issue first for anything non-trivial.

- Run `cargo test --release` (and `cargo run --release --example det_hash -- 2 4 64 200000 50000 64`, which must print lines present in `det/ref.txt`) before sending changes.
- Keep exactness guarantees intact: integer arithmetic only in the state, no float in the update path.
- Licensing of contributions: the project is source-available under the PolyForm Noncommercial License 1.0.0 and the maintainer also grants commercial licences (see `COMMERCIAL.md`). By submitting a contribution you confirm you wrote it (or have the right to submit it) and you grant the maintainer a perpetual, worldwide, irrevocable licence to use, modify, sublicense and relicense it, including under commercial terms. If you cannot agree to that, open an issue describing the change instead of a pull request.

## Releasing (maintainer)
1. Bump the version in `Cargo.toml`, `python/Cargo.toml`, `python/pyproject.toml` and `CITATION.cff`; add a `## X.Y.Z` section to `CHANGELOG.md`; merge to `main` with CI green.
2. `git tag -a vX.Y.Z -m "sigring X.Y.Z" && git push origin vX.Y.Z`.
3. The `release` workflow builds and tests wheels, then waits for approval of the `release` environment (crates.io) and the `pypi` environment (PyPI). No API tokens are stored; both registries trust this repository's workflow via OIDC.
4. A manual run of the workflow (Actions -> release -> Run workflow) builds everything and publishes nothing.
