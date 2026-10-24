# sigring (Python)

Exact, streaming, sliding-window truncated path signatures. Python bindings for the Rust crate [sigring](https://github.com/metapozza/sigring).

Source-available under the PolyForm Noncommercial License 1.0.0 (free for personal, research, educational and non-profit use; commercial use needs a licence, see the repository's `COMMERCIAL.md`).

```python
import numpy as np
from sigring import WindowSig

K = 1024                                   # scale: real increments * K -> integers
inc = np.random.default_rng(0).integers(-100, 101, (5000, 2)).astype(np.int64)
ws = WindowSig(d=2, n=4, window=64, scale=float(K), bits=64)   # bits=128 for headroom
sigs = ws.push_batch(inc, "sig", 1)        # one signature row per tick
print(sigs.shape, ws.guard_ok())           # guard_ok(): rigorous overflow certificate
```

Readouts: `"sig"`, `"logsig_expanded"`, `"logsig_words"`, `"logsig_brackets"`. See the repository README for guarantees, benchmarks and limitations.

## Quant use, support, licensing

- Example with a per-window exactness check: `python/examples/quant_ticks.py` in the repository.
- Check `guard_ok()` for every window you use; `bits=128` gives headroom for clock gaps and large increments.
- Questions: https://github.com/metapozza/sigring/discussions . Evaluation help, commercial licences, roles: metapozza@proton.me
