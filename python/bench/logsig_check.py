import random, subprocess, numpy as np, pysiglib, iisignature
from fractions import Fraction
def run(d,n,W,K,stride,incs):
    inp=f"{d} {n} {W} {K} {stride}\n"+"\n".join(" ".join(map(str,r)) for r in incs)+"\n"
    out=subprocess.run(["sigring/target/release/logsig"],input=inp,capture_output=True,text=True).stdout.strip().split("\n")
    res=[]
    for i in range(0,len(out),3):
        res.append(tuple(np.array(list(map(float,l.split()))) for l in out[i:i+3]))
    return res
# exact oracle: signature of piecewise-linear path via exp(inc) concatenation in Fractions, then tensor log in Fractions
def tmul(a,b,d,n):
    out=[[Fraction(0)]*(d**m) for m in range(n+1)]
    for i in range(n+1):
        for j in range(n+1-i):
            for x,u in enumerate(a[i]):
                if u==0: continue
                for y,v in enumerate(b[j]): out[i+j][x*len(b[j])+y]+=u*v
    return out
def expo(v,d,n):
    out=[[Fraction(1)]]
    for m in range(1,n+1):
        prev=out[-1]; out.append([p*Fraction(v[c],1)/m for p in prev for c in range(d)])
    return out
def exact_log(win,d,n):
    S=[[Fraction(1)]]+[[Fraction(0)]*(d**m) for m in range(1,n+1)]
    for v in win: S=tmul(S,expo(v,d,n),d,n)
    X=[[Fraction(0)]]+S[1:]; P=[[Fraction(1)]]+[[Fraction(0)]*(d**m) for m in range(1,n+1)]
    L=[[Fraction(0)]*(d**m) for m in range(n+1)]
    for k in range(1,n+1):
        P=tmul(P,X,d,n)
        for m in range(n+1):
            for i,v in enumerate(P[m]): L[m][i]+=Fraction((-1)**(k+1),k)*v
    return L
random.seed(3)
for (d,n,W,K) in [(2,4,16,4096),(3,3,12,4096),(2,5,10,2048),(4,3,8,1024)]:
    T=W+9; incs=[[random.randint(-3*K,3*K) for _ in range(d)] for _ in range(T)]
    res=run(d,n,W,K,1,incs)
    # pySigLib/iisignature on the same real-valued path
    worst_ex=0; worst_pl=[0,0,0]; worst_ii=0
    for t,(ex,wc,bc) in enumerate(res):
        win=incs[t:t+W]
        L=exact_log(win,d,n)
        exp=np.array([float(c)/K**m for m in range(1,n+1) for c in L[m]])
        sc=np.array([max(abs(float(c))/K**m for c in L[m]) or 1 for m in range(1,n+1) for _ in range(d**m)])
        worst_ex=max(worst_ex,np.max(np.abs(ex-exp)/sc))
        path=np.vstack([np.zeros(d),np.cumsum(np.array(win,float)/K,axis=0)])
        e0=pysiglib.log_sig(path,n,method=0)  # expanded, no scalar
        worst_pl[0]=max(worst_pl[0],np.max(np.abs(e0-exp)/sc))
        # iisignature expanded
        s=iisignature.prepare(d,n,"x"); ii=iisignature.logsig(path,s,"x")
        worst_ii=max(worst_ii,np.max(np.abs(ii-exp)/sc))
    print(f"d={d} N={n} W={W}: ring-expanded-log vs exact {worst_ex:.2e} | pySigLib m0 vs exact {worst_pl[0]:.2e} | iisignature x vs exact {worst_ii:.2e}")
# basis conventions: compare Lyndon projections to pySigLib methods 1 and 2
d,n,W,K=2,4,16,4096
incs=[[random.randint(-3*K,3*K) for _ in range(d)] for _ in range(W)]
res=run(d,n,W,K,1,incs); ex,wc,bc=res[0]
path=np.vstack([np.zeros(d),np.cumsum(np.array(incs,float)/K,axis=0)])
for m in (1,2):
    pysiglib.prepare_log_sig(d,n,method=m)
    o=pysiglib.log_sig(path,n,method=m)
    print("method",m,"len",len(o),"vs words",np.max(np.abs(o-wc)/(abs(wc).max())),"vs brackets",np.max(np.abs(o-bc)/(abs(bc).max())))
print("--- d=3,N=4 basis discrimination")
d,n,W,K=3,4,10,4096
incs=[[random.randint(-3*K,3*K) for _ in range(d)] for _ in range(W)]
res=run(d,n,W,K,1,incs); ex,wc,bc=res[0]
path=np.vstack([np.zeros(d),np.cumsum(np.array(incs,float)/K,axis=0)])
print("mine lens",len(wc),len(bc),"pysig len",pysiglib.log_sig_length(d,n))
for m in (1,2):
    pysiglib.prepare_log_sig(d,n,method=m)
    o=pysiglib.log_sig(path,n,method=m)
    print("method",m,"vs words",np.max(np.abs(o-wc))/abs(wc).max(),"vs brackets",np.max(np.abs(o-bc))/abs(bc).max())
