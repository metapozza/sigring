import numpy as np, time, pysiglib, iisignature, sys
from numpy.lib.stride_tricks import sliding_window_view
def best(f,rep=5):
    ts=[]
    for _ in range(rep):
        t=time.perf_counter(); f(); ts.append(time.perf_counter()-t)
    return min(ts)
rows=[]
for (d,n,W) in [(2,4,64),(2,3,64),(3,3,64),(3,4,64),(5,3,64),(2,5,64),(2,4,256),(2,4,1024)]:
    B=max(64,200000//W//4)
    rng=np.random.default_rng(0); X=np.cumsum(rng.integers(-1000,1001,(B+W,d))/1024.0,axis=0)
    wins=np.ascontiguousarray(sliding_window_view(X,(W+1,d))[:,0][:B])   # (B, W+1, d)
    out=[f"{d},{n},{W}"]
    for m in (1,2):
        pysiglib.prepare_log_sig(d,n,method=m)
        t=best(lambda: pysiglib.log_sig(wins,n,method=m,n_jobs=1)); out.append(f"{t/B*1e9:.0f}")
    s=iisignature.prepare(d,n,"S")
    t=best(lambda: iisignature.logsig(wins,s,"S")); out.append(f"{t/B*1e9:.0f}")
    print("d,N,W | pySigLib m1 | pySigLib m2 | iisignature S  (ns/window):"," | ".join(out),flush=True)
