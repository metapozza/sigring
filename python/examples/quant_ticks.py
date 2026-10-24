"""Windowed path-signature features from a tick stream, with a per-window exactness check.
Synthetic data, so it runs without any market data.

Path = (price in ticks, elapsed milliseconds), both integers, so scale=1.0.
For real data, quantise first: price -> round(price / tick_size), time -> integer milliseconds.

    pip install sigring numpy
    python quant_ticks.py
"""
import numpy as np
from sigring import WindowSig

rng = np.random.default_rng(7)
T, D, N, W, STRIDE = 20_000, 2, 3, 128, 16

price_ticks = rng.integers(-3, 4, T)           # price change per tick, in ticks
dt_ms = rng.integers(50, 150, T)               # time between ticks, ms
dt_ms[10_000] = 23_400_000                     # one overnight gap (6.5 hours)
inc = np.stack([price_ticks, dt_ms], axis=1).astype(np.int64)


def features(bits):
    """Push tick by tick; keep a row every STRIDE ticks and record whether the guard certifies that window."""
    ws = WindowSig(d=D, n=N, window=W, scale=1.0, bits=bits)
    rows, certified = [], []
    for t in range(T):
        ws.push(inc[t])
        if ws.full() and (t - W + 1) % STRIDE == 0:
            rows.append(ws.read("logsig_brackets"))
            certified.append(ws.guard_ok())
    return np.array(rows), np.array(certified), ws


# 64-bit state is fastest but cannot hold windows that contain the 6.5-hour gap; the guard says so.
f64, ok64, _ = features(64)
f128, ok128, ws = features(128)
print(f"feature matrix {f128.shape}: {f128.shape[1]} Lyndon-basis coordinates per window")
print(f"i64  guard certifies {ok64.sum()}/{len(ok64)} windows -> drop or recompute the rest")
print(f"i128 guard certifies {ok128.sum()}/{len(ok128)} windows")
assert ok128.all()
# Where both are certified they are bit-for-bit the same features.
assert np.array_equal(f64[ok64], f128[ok64])

# Built-in sanity check: level-1 of the signature is the plain sum of increments in the window.
sig = ws.read("sig")
assert np.allclose(sig[:D], inc[-W:].sum(axis=0))
print("level-1 check against window sums: ok")
# Columns of f128: price displacement, time displacement, then higher-order (area-type) terms. Feed them to your model.
