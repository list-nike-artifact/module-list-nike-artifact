//! Test vectors for checking this port against the original x86-64 build.
//!
//! Everything printed here is a deterministic function of fixed inputs, so the
//! same binary built from upstream `ref0` (with the asm field arithmetic and the
//! AES-NI keystream) must print exactly the same bytes.  Run
//!
//!     cargo run --release --bin kat
//!
//! here and in the x86 tree, and diff the two outputs.

use pswoosh::*;
use pswoosh::arithmetic::{fq::*, poly::*};

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{:02x}", x)).collect()
}

fn elem(e: Elem) -> String {
    format!("{:016x}{:016x}{:016x}{:016x}", e[3], e[2], e[1], e[0])
}

/* A few elements that exercise carries and reductions rather than small values */
fn probes() -> Vec<Elem> {
    let mut v = vec![
        [0, 0, 0, 0],
        [1, 0, 0, 0],
        QQ, HQ, TQQ,
        [Q[0]-1, Q[1], Q[2], Q[3]],
        [0xffffffffffffffff, 0, 0, 0],
        [0, 0xffffffffffffffff, 0xffffffffffffffff, 0x1234],
        [0x0123456789abcdef, 0xfedcba9876543210, 0x5555555555555555, 0x2aaaaa],
    ];
    /* a deterministic spread over the whole range */
    let mut x: u64 = 0x9e3779b97f4a7c15;
    for _ in 0..8 {
        let mut e: Elem = [0; NLIMBS];
        for i in 0..NLIMBS {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            e[i] = x;
        }
        e[NLIMBS-1] &= 0x3fffff;
        if cmp(e, Q) != 0x80 {
            sub(&mut e.clone(), e, Q);
        }
        v.push(e);
    }
    v
}

fn main() {
    pswoosh::with_stack(run);
}

fn run() {
    println!("# field arithmetic");
    let p = probes();
    for a in &p {
        for b in &p {
            let mut c = fp_init();
            add(&mut c, *a, *b);
            print!("{} ", elem(c));
            sub(&mut c, *a, *b);
            print!("{} ", elem(c));
            mul(&mut c, *a, *b);
            println!("{}", elem(c));
        }
        let mut m = fp_init();
        let mut r = fp_init();
        toM(&mut m, *a);
        fromM(&mut r, m);
        println!("toM {} fromM {}", elem(m), elem(r));
    }

    println!("# matrix");
    let a: Matrix = genmatrix(&[0x42; SYMBYTES], true);
    let at: Matrix = genmatrix(&[0x42; SYMBYTES], false);
    println!("A[0][0]   {}", elem(a[0][0][0]));
    println!("A[31][31] {}", elem(a[31][31][255]));
    println!("A row0    {}", hex(&poly_tobytes(a[0][0])[..32]));

    println!("# keystream");
    let mut buf = [0u8; NOISE_BYTES];
    expand_seed_aes(&[0x5a; SYMBYTES], 7, &mut buf);
    println!("aes {} .. {}", hex(&buf[..32]), hex(&buf[NOISE_BYTES-16..]));

    println!("# noise");
    let s = getnoise(&[0x5a; SYMBYTES], 0);
    println!("s[0][0..4] {} {} {} {}", elem(s[0][0]), elem(s[0][1]), elem(s[0][2]), elem(s[0][3]));

    println!("# scheme");
    let (sk1, pk1) = pswoosh_keygen_seeded(&a,  &[0x01; SYMBYTES], true);
    let (sk2, pk2) = pswoosh_keygen_seeded(&at, &[0x02; SYMBYTES], false);
    println!("sk1 {}", hex(&sk1[..32]));
    println!("pk1 {}", hex(&pk1[..32]));
    println!("pk1 tail {}", hex(&pk1[PUBLICKEY_BYTES-32..]));
    println!("pk2 {}", hex(&pk2[..32]));
    let ss1 = pswoosh_skey_deriv(&pk1, &pk2, &sk1, true);
    let ss2 = pswoosh_skey_deriv(&pk2, &pk1, &sk2, false);
    println!("ss1 {}", hex(&ss1));
    println!("ss2 {}", hex(&ss2));
    println!("agree {}", ss1 == ss2);
    println!("|pk| {} |sk| {} |ss| {}", PUBLICKEY_BYTES, SECRETKEY_BYTES, SYMBYTES);
}
