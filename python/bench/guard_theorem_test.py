# Test of the guard theorem:  for any window of integer increments x_1..x_W in Z^d,  every entry of the factorial-cleared state t_m = m! * S_m
# satisfies |t_m[i]| <= (sum_i ||x_i||_1)^m.  (Exact Python integers, no rings.)  Also reports tightness and the rate of equality.
import random, itertools, math, sys
from math import comb
def tmul_exp_right(t, r, d, N):
    # t <- t (x) E(r) in the scaled domain: t'_m = sum_k C(m,k) t_k (x) r^{(m-k)}  (t_0 = 1)
    pw=[None,[v for v in r]]
    for k in range(2,N+1): pw.append([a*b for a in pw[k-1] for b in r])
    new=[None]+[[0]*(d**m) for m in range(1,N+1)]
    for m in range(1,N+1):
        for i,v in enumerate(pw[m]): new[m][i]+=v                 # k=0
        for k in range(1,m+1):
            lb=d**(m-k)
            c=comb(m,k)
            if k==m:
                for i,v in enumerate(t[m]): new[m][i]+=v
            else:
                for a,x in enumerate(t[k]):
                    if x==0: continue
                    for b,y in enumerate(pw[m-k][1:] if False else pw[m-k]): new[m][a*lb+b]+=c*x*y
    return new
def window_state(win,d,N):
    t=[None]+[[0]*(d**m) for m in range(1,N+1)]
    for r in win: t=tmul_exp_right(t,r,d,N)
    return t
random.seed(11)
trials=0; viol=0; eq=0; worst_ratio=0.0
gens=[lambda: random.randint(-5,5), lambda: random.choice([-1,0,1]), lambda: int(random.paretovariate(1.2))*random.choice([-1,1]), lambda: random.randint(-1000,1000), lambda: 3]
for d in (1,2,3):
    for N in (2,3,4):
        for W in (1,2,5,9):
            for g in range(len(gens)):
                for _ in range(60):
                    win=[[gens[g]() for _ in range(d)] for _ in range(W)]
                    L1=sum(abs(x) for r in win for x in r)
                    t=window_state(win,d,N)
                    for m in range(1,N+1):
                        mx=max(abs(v) for v in t[m]); bound=L1**m; trials+=1
                        if mx>bound: viol+=1; print("VIOLATION",d,N,W,m,mx,bound)
                        if mx==bound: eq+=1
                        if bound>0: worst_ratio=max(worst_ratio,mx/bound)
print(f"checked {trials} (window,level) pairs over d in 1..3, N in 2..4, W in 1..9, 5 increment distributions")
print("violations of |t_m| <= L1^m :",viol,"| cases attaining equality:",eq,"| max observed ratio |t_m|/L1^m:",worst_ratio)
