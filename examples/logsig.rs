use sigring::{SigRing, LogSigPrep};
use std::io::Read;
// stdin: d n w K stride, then T lines of d integer increments (already scaled by K). Prints log-sigs of each full window every `stride` ticks.
fn main(){
    let mut s=String::new(); std::io::stdin().read_to_string(&mut s).unwrap();
    let mut it=s.split_ascii_whitespace();
    let d:usize=it.next().unwrap().parse().unwrap(); let n:usize=it.next().unwrap().parse().unwrap(); let w:usize=it.next().unwrap().parse().unwrap();
    let k:f64=it.next().unwrap().parse().unwrap(); let stride:usize=it.next().unwrap().parse().unwrap();
    let vals:Vec<i64>=it.map(|x|x.parse().unwrap()).collect();
    let t=vals.len()/d;
    let mut e:SigRing<i128>=SigRing::new(d,n); let prep=LogSigPrep::new(d,n);
    let mut ring=vec![vec![0i64;d];w]; let mut h=0usize;
    for i in 0..t{
        let inc=&vals[i*d..(i+1)*d]; let old=ring[h].clone(); ring[h].copy_from_slice(inc); h=(h+1)%w; e.ingest(inc,&old);
        if i+1>=w && (i+1-w)%stride==0 {
            let ex=e.log_signature_expanded(k);
            let wc=prep.lyndon_word_coefs(&ex); let bc=prep.bracket_coords(&ex);
            let f=|v:&Vec<f64>|v.iter().map(|x|format!("{:.17e}",x)).collect::<Vec<_>>().join(" ");
            println!("{}\n{}\n{}",f(&ex),f(&wc),f(&bc));
        }
    }
}
