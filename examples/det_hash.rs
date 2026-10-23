use sigring::SigRing;
fn fnv<I: Iterator<Item = u128>>(it: I, nb: usize) -> u64 { let mut h: u64 = 1469598103934665603; for mut v in it { for _ in 0..nb { h ^= (v & 0xff) as u64; h = h.wrapping_mul(1099511628211); v >>= 8; } } h }
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let d: usize = a[1].parse().unwrap(); let n: usize = a[2].parse().unwrap(); let w: usize = a[3].parse().unwrap();
    let t: u64 = a[4].parse().unwrap(); let every: u64 = a[5].parse().unwrap(); let bits: u32 = a[6].parse().unwrap();
    let mut e64: SigRing<i64> = SigRing::new(d, n); let mut e128: SigRing<i128> = SigRing::new(d, n);
    let mut ring = vec![0i64; w * d]; let mut h = 0usize; let mut s: u64 = 12345; let amp: i64 = 1000;
    let mut inc = vec![0i64; d];
    for tk in 1..=t {
        for c in 0..d { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); inc[c] = ((s >> 33) % (2 * amp as u64 + 1)) as i64 - amp; }
        let old: Vec<i64> = ring[h * d..(h + 1) * d].to_vec(); ring[h * d..(h + 1) * d].copy_from_slice(&inc); h = (h + 1) % w;
        if bits == 64 { e64.ingest(&inc, &old); } else { e128.ingest(&inc, &old); }
        if tk % every == 0 {
            let hh = if bits == 64 { fnv(e64.t.iter().map(|&v| v as u64 as u128), 8) } else { fnv(e128.t.iter().map(|&v| v as u128), 16) };
            println!("{} {} {} {} {} {:016x}", d, n, w, bits, tk, hh);
        }
    }
}
