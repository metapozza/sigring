import numpy as np, time, sigring, pysiglib
from sigring import WindowSig
def tm(f,rep=5):
    b=1e9
    for _ in range(rep):
        t=time.perf_counter(); f(); b=min(b,time.perf_counter()-t)
    return b
rng=np.random.default_rng(1)
# correctness: batch sig == per-window pySigLib sig
d,n,W=2,4,64; K=1024; T=5000
inc=rng.integers(-1000,1001,(T,d)).astype(np.int64)
ws=WindowSig(d,n,W,float(K),64)
out=ws.push_batch(inc,"sig",1)
path=np.vstack([np.zeros(d),np.cumsum(inc/K,axis=0)])
ref=np.array([pysiglib.sig(path[t:t+W+1]-path[t],n) for t in (0,100,2000,T-W)])
idx=[0,100,2000,T-W]
print("rows",out.shape,"max abs err vs pySigLib sig:",max(np.max(np.abs(out[i]-ref[j])/np.max(np.abs(ref[j]))) for j,i in enumerate(idx)))
# overhead
for (d,n,W) in [(2,4,64),(3,3,64),(2,4,1024)]:
    T=200000; inc=rng.integers(-1000,1001,(T,d)).astype(np.int64)
    ws=WindowSig(d,n,W,1024.0,64)
    t_ing=tm(lambda: ws.ingest_batch(inc))/T*1e9
    ws=WindowSig(d,n,W,1024.0,64)
    t_sig=tm(lambda: ws.push_batch(inc,"sig",1))/T*1e9
    ws=WindowSig(d,n,W,1024.0,64)
    t_lsb=tm(lambda: ws.push_batch(inc,"logsig_brackets",1))/T*1e9
    ws=WindowSig(d,n,W,1024.0,64); rows=[inc[i] for i in range(20000)]
    def loop():
        for r in rows: ws.push(r)
    t_push=tm(loop,3)/20000*1e9
    ws=WindowSig(d,n,W,1024.0,64)
    def loop2():
        for r in rows: ws.push(r); ws.read("sig")
    t_pr=tm(loop2,3)/20000*1e9
    print(f"d={d} N={n} W={W}: ingest_batch {t_ing:.0f} ns/tick | push_batch sig {t_sig:.0f} | push_batch logsig_brackets {t_lsb:.0f} | python loop push() {t_push:.0f} | push()+read() {t_pr:.0f} ns/tick")
# pySigLib streaming per tick from Python for reference
d,n,W=2,4,64
pysiglib.prepare_log_sig(d,n,method=2)
ss=pysiglib.SigWindowStream(d,n,W)
pts=np.cumsum(rng.standard_normal((3000,d)),axis=0)
t=time.perf_counter()
for p in pts[:2000]: ss.push(p)
print("pySigLib SigWindowStream python loop ns/tick:",(time.perf_counter()-t)/2000*1e9)
t=time.perf_counter(); ss2=pysiglib.SigWindowStream(d,n,W); ss2.push_batch(pts); print("pySigLib SigWindowStream push_batch ns/tick:",(time.perf_counter()-t)/3000*1e9)
