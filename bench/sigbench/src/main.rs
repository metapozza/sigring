use sigring::SigRing;
use std::time::Instant;
// ---------- float implementation of the same scaled-domain algebra ----------
struct Lay { d:usize, n:usize, off:Vec<usize>, binom:Vec<Vec<f64>> }
impl Lay {
    fn new(d:usize,n:usize)->Self{ let mut off=vec![0usize;n+1]; let mut sz=1; let mut acc=0; for k in 1..=n{sz*=d; off[k-1]=acc; acc+=sz; off[k]=acc;}
        let mut b=vec![vec![0f64;n+1];n+1]; for i in 0..=n{b[i][0]=1.0; for j in 1..=i{b[i][j]=b[i-1][j-1]+if j<=i-1{b[i-1][j]}else{0.0};}} Self{d,n,off,binom:b} }
    fn total(&self)->usize{ self.off[self.n] }
}
fn powers(l:&Lay,v:&[f64],buf:&mut [f64]){ let d=l.d; buf[l.off[0]..l.off[1]].copy_from_slice(v);
    for k in 2..=l.n{ let (lo,hi)=buf.split_at_mut(l.off[k-1]); let prev=&lo[l.off[k-2]..l.off[k-1]]; let cur=&mut hi[..l.off[k]-l.off[k-1]];
        for (ia,&x) in prev.iter().enumerate(){ for ib in 0..d{ cur[ia*d+ib]=x*v[ib]; } } } }
fn outer_acc(coef:f64,a:&[f64],b:&[f64],out:&mut [f64]){ let lb=b.len(); for (ia,&x) in a.iter().enumerate(){ let cx=coef*x; let row=&mut out[ia*lb..(ia+1)*lb]; for (o,&y) in row.iter_mut().zip(b){ *o+=cx*y; } } }
// right-multiply S by E(r):  S <- S (x) E(r)
fn rmul_exp(l:&Lay,s:&mut [f64],pw:&mut [f64],r:&[f64]){ powers(l,r,pw);
    for m in (1..=l.n).rev(){ let len=l.off[m]-l.off[m-1]; for i in l.off[m-1]..l.off[m]{ s[i]+=pw[i]; }
        for k in 1..m{ let (a,b)=(l.off[k-1]..l.off[k],l.off[m-k-1]..l.off[m-k]); let (lo,hi)=s.split_at_mut(l.off[m-1]); outer_acc(l.binom[m][k],&lo[a],&pw[b],&mut hi[..len]); } } }
// left-multiply S by E(v):  S <- E(v) (x) S
fn lmul_exp(l:&Lay,s:&mut [f64],pw:&mut [f64],v:&[f64]){ powers(l,v,pw);
    for m in (1..=l.n).rev(){ let len=l.off[m]-l.off[m-1]; for i in l.off[m-1]..l.off[m]{ s[i]+=pw[i]; }
        for j in 1..m{ let (a,b)=(l.off[j-1]..l.off[j],l.off[m-j-1]..l.off[m-j]); let (lo,hi)=s.split_at_mut(l.off[m-1]); outer_acc(l.binom[m][j],&pw[a],&lo[b],&mut hi[..len]); } } }
// general product C = A (x) B  (scaled domain)
fn gmul(l:&Lay,a:&[f64],b:&[f64],c:&mut [f64]){ for m in 1..=l.n{ let len=l.off[m]-l.off[m-1]; let cm=&mut c[l.off[m-1]..l.off[m]]; cm.iter_mut().for_each(|x|*x=0.0);
        for i in 0..len{ let _=i; }
        for k in 0..=m{ let coef=l.binom[m][k];
            if k==0{ for (o,&y) in cm.iter_mut().zip(&b[l.off[m-1]..l.off[m]]){*o+=y;} }
            else if k==m{ for (o,&x) in cm.iter_mut().zip(&a[l.off[m-1]..l.off[m]]){*o+=x;} }
            else{ outer_acc(coef,&a[l.off[k-1]..l.off[k]],&b[l.off[m-k-1]..l.off[m-k]],cm); } } } }
// ---------- stream generator ----------
struct Gen{s:u64,d:usize,amp:i64,drift:i64,kind:u8,cnt:u64,file:Vec<i64>}
impl Gen{ fn u(&mut self)->u64{ self.s=self.s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); self.s>>11 }
  fn next(&mut self,out:&mut [i64]){ self.cnt+=1; if self.kind==9 { let p=((self.cnt-1) as usize)*self.d; for c in 0..self.d{ out[c]=self.file[p+c]; } return; } for c in 0..self.d{ let r=self.u(); let mut v=((r%(2*self.amp as u64+1)) as i64)-self.amp;
        if self.kind==2{ // heavy tail: pareto-ish
            let uu=((self.u()%1_000_000)+1) as f64/1e6; let m=(1.0/uu.powf(0.7)).min(2000.0); v=(v as f64*m/ 3.0).round() as i64; }
        if c==0{v+=self.drift;} out[c]=v; }
    if self.kind==3 && self.cnt%1_000_003==0 { out[0]=self.amp*20000; out[1]=-self.amp*15000; } } }
fn main(){
    let a:Vec<String>=std::env::args().collect();
    let mode=a[1].as_str(); let method=a[2].as_str();
    let d:usize=a[3].parse().unwrap(); let n:usize=a[4].parse().unwrap(); let w:usize=a[5].parse().unwrap(); let ticks:u64=a[6].parse().unwrap();
    let kind:u8=a[7].parse().unwrap(); // 0 rw,1 trend,2 heavy,3 spikes
    let amp:i64=a.get(8).map(|x|x.parse().unwrap()).unwrap_or(5000); let drift:i64=if kind==1{2000}else{0};
    let lay=Lay::new(d,n); let tot=lay.total();
    let file:Vec<i64>= if kind==9 { std::fs::read_to_string(std::env::var("SIG_FILE").unwrap()).unwrap().split_ascii_whitespace().map(|x|x.parse().unwrap()).collect() } else {vec![]};
    let ticks= if kind==9 { ticks.min((file.len()/d) as u64) } else { ticks };
    let mut g=Gen{s:99,d,amp,drift,kind,cnt:0,file};
    let mut ring=vec![0i64;w*d]; let mut h=0usize; let mut inc=vec![0i64;d];
    let mut ig:SigRing<i128>=SigRing::new(d,n); let mut ig64:SigRing<i64>=SigRing::new(d,n); let mut l64=vec![0i64;d]; let mut r64=vec![0i64;d];
    let mut l128=vec![0i128;d]; let mut r128=vec![0i128;d];
    let mut fh:SigRing<f64>=SigRing::new(d,n); let mut bh:SigRing<f64>=SigRing::new(d,n); let mut th:SigRing<f64>=SigRing::new(d,n); let zero=vec![0f64;d]; let mut front_h:Vec<Vec<f64>>=vec![]; let mut back_h:Vec<Vec<f64>>=vec![]; let mut qh=0usize; let mut outh=vec![0f64;tot];
    let mut sf=vec![0f64;tot]; let mut pw=vec![0f64;tot]; let mut lf=vec![0f64;d]; let mut rf=vec![0f64;d];
    // two-stack state
    let mut front:Vec<Vec<f64>>=vec![]; let mut back:Vec<Vec<f64>>=vec![]; let mut bagg=vec![0f64;tot]; let mut qlen=0usize; let mut outv=vec![0f64;tot];
    let mut incs_hist:Vec<f64>=vec![]; // for twostack flip: store increments
    let mut back_incs:Vec<Vec<f64>>=vec![];
    let checkpoints:Vec<u64>= if mode=="acc"{ (1..=8).map(|i| ticks*i/8).collect() } else {vec![]};
    let mut ci=0usize;
    let t0=Instant::now();
    let mut lat:Vec<u32>=Vec::with_capacity(ticks as usize);
    for t in 0..ticks{
        let tt=Instant::now();
        g.next(&mut inc);
        let oldo=h*d; let has_old=(t as usize)>=w;
        match method{
            "int64"=>{ for c in 0..d{ r64[c]=inc[c]; l64[c]= if has_old {-ring[oldo+c]} else {0}; } ig64.ingest_i128(&l64,&r64); }
            "int"=>{ for c in 0..d{ r128[c]=inc[c] as i128; l128[c]= if has_old {-(ring[oldo+c] as i128)} else {0}; } ig.ingest_i128(&l128,&r128); }
            "fused_float"=>{ for c in 0..d{ rf[c]=inc[c] as f64; lf[c]= if has_old {-(ring[oldo+c] as f64)} else {0.0}; } rmul_exp(&lay,&mut sf,&mut pw,&rf); if has_old{ lmul_exp(&lay,&mut sf,&mut pw,&lf);} }
            "twostack"=>{ for c in 0..d{ rf[c]=inc[c] as f64; }
                // push: bagg = bagg (x) E(x)
                rmul_exp(&lay,&mut bagg,&mut pw,&rf); back_incs.push(rf.clone()); qlen+=1;
                if has_old{ qlen-=1;
                    if front.is_empty(){ // flip: suffix aggregates
                        let mut agg=vec![0f64;tot]; let mut stack:Vec<Vec<f64>>=Vec::with_capacity(back_incs.len());
                        for x in back_incs.iter().rev(){ lmul_exp(&lay,&mut agg,&mut pw,x); stack.push(agg.clone()); }
                        front=stack; back_incs.clear(); bagg.iter_mut().for_each(|x|*x=0.0); }
                    front.pop(); }
                if mode!="acc" || true { let ft=front.last(); match ft{ Some(f)=>gmul(&lay,f,&bagg,&mut outv), None=>outv.copy_from_slice(&bagg)} } }
            "fused_float_h"=>{ for c in 0..d{ rf[c]=inc[c] as f64; lf[c]= if has_old {-(ring[oldo+c] as f64)} else {0.0}; } fh.ingest_i128(&lf,&rf); }
            "twostack_h"=>{ for c in 0..d{ rf[c]=inc[c] as f64; }
                bh.ingest_i128(&zero,&rf); back_h.push(rf.clone()); qh+=1;
                if has_old{ qh-=1;
                    if front_h.is_empty(){
                        for x in th.t.iter_mut(){*x=0.0;} let mut stack:Vec<Vec<f64>>=Vec::with_capacity(back_h.len());
                        for x in back_h.iter().rev(){ th.ingest_i128(x,&zero); stack.push(th.t.clone()); }
                        front_h=stack; back_h.clear(); for x in bh.t.iter_mut(){*x=0.0;} }
                    front_h.pop(); }
                match front_h.last(){ Some(f)=>gmul(&lay,f,&bh.t,&mut outh), None=>outh.copy_from_slice(&bh.t)} }
            "recompute"=>{ if t as usize+1>=w { sf.iter_mut().for_each(|x|*x=0.0); // window = last w increments incl current
                    for k in 0..w{ let idx=(h+1+k)%w; for c in 0..d{ rf[c]= if k==w-1 {inc[c] as f64} else {ring[idx*d+c] as f64}; } rmul_exp(&lay,&mut sf,&mut pw,&rf); } } }
            _=>panic!("method"),
        }
        ring[oldo..oldo+d].copy_from_slice(&inc); h=(h+1)%w;
        lat.push(tt.elapsed().as_nanos() as u32);
        if ci<checkpoints.len() && t+1==checkpoints[ci]{ ci+=1;
            // print window oldest first then state (float or int)
            let mut line=String::new(); for i in 0..w{ let idx=(h+i)%w; for c in 0..d{ line.push_str(&format!("{} ",ring[idx*d+c])); } }
            println!("{}",line.trim());
            let st:Vec<String>=match method{ "int"=>ig.t.iter().map(|v|v.to_string()).collect(), "int64"=>ig64.t.iter().map(|v|v.to_string()).collect(), "fused_float"=>sf.iter().map(|v|format!("{:e}",v)).collect(), "twostack"=>outv.iter().map(|v|format!("{:e}",v)).collect(), "fused_float_h"=>fh.t.iter().map(|v|format!("{:e}",v)).collect(), "twostack_h"=>outh.iter().map(|v|format!("{:e}",v)).collect(), _=>sf.iter().map(|v|format!("{:e}",v)).collect() };
            println!("{}",st.join(" ")); }
    }
    let el=t0.elapsed().as_secs_f64();
    if std::env::var("SIG_LAT").is_ok() { let mut v=lat.clone(); v.sort_unstable(); let q=|p:f64| v[((v.len()-1) as f64*p) as usize]; println!("LAT ns: p50 {} p99 {} p99.9 {} p99.99 {} max {}",q(0.5),q(0.99),q(0.999),q(0.9999),v[v.len()-1]); }
    if mode=="bench"{ // consume
        let s=match method{"int"=>ig.t[0] as f64,"int64"=>ig64.t[0] as f64,"fused_float"=>sf[0],"twostack"=>outv[0],"fused_float_h"=>fh.t[0],"twostack_h"=>outh[0],_=>sf[0]}; println!("{:.1} ns/tick (chk {:e})",el*1e9/ticks as f64,s); }
}
