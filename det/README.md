# Determinism check

- `ref.txt`: reference state hashes (FNV-1a over the ring state every 50,000 ticks) for four (d, N, W) configs at 64 and 128 bits. Line format: `d N W bits tick hash`.
- Rust: `cargo run --release --example det_hash -- <d> <N> <W> <ticks> <every> <bits>`.
- `det.c`: independent C port of the same update (`gcc -O2 -DWORD=64|128 det.c`; args `d N W ticks every`).
- `det_wasm.c`, `run_wasm.js`: wasm32 build/run harness used with node.
- Cross-architecture runs (aarch64, riscv64, ppc64 BE, s390x BE, arm32, mips32 BE (64-bit only), wasm32) were done with qemu-user and cross GCC toolchains and compared against `ref.txt`.
