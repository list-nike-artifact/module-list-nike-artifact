//! Prints, in one page, everything about this implementation that can be
//! checked against the paper by eye.

use mlnike::listrec;
use mlnike::nike;
use mlnike::params::*;
use mlnike::{fq, poly};

fn main() {
    println!("parameters (tab:mod-params, recommended row)");
    println!("  d      = {D}");
    println!("  rho    = {RHO}");
    println!("  c      = {C}");
    println!("  beta   = {BETA}");
    println!("  n      = {N}");
    println!("  m      = {M}");
    println!("  q      = 2^{LOG_Q} - {} = {Q}", (1u32 << LOG_Q) - Q);
    println!("  eta    = {ETA}");
    println!("  log l  = {LOG_ELL}   (list size l = {ELL})");

    println!("\nhypothesis eq:mod-H:  2 m d eta beta <= q / 2^(rho+1)");
    let lhs = 2 * M as u64 * D as u64 * ETA as u64 * BETA as u64;
    let rhs = Q as u64 / (1 << (RHO + 1));
    println!("  {lhs} <= {rhs}   {}", if lhs <= rhs { "ok" } else { "FAILS" });

    let pp = nike::stp(&[0x5au8; 16]);
    let (ska, pka) = nike::gen(&pp, &[1u8; 16]);
    let (skb, pkb) = nike::gen(&pp, &[2u8; 16]);

    println!("\nsizes: serialise the keys just generated and measure");
    let pk_len = nike::pk_bytes(&pka).len();
    let sk_len = nike::sk_bytes(&ska).len();
    let shk_bits = nike::sdk_first(&ska, &pkb).len() * 8;
    println!("  |pk|  = {pk_len} bytes   table says {}   {}", table::PK_BYTES, ok(pk_len == table::PK_BYTES));
    println!("  |sk|  = {sk_len} bytes   table says {}   {}", table::SK_BYTES, ok(sk_len == table::SK_BYTES));
    println!("  |shk| = {shk_bits} bits    table says {}     {}", table::SHK_BITS, ok(shk_bits == table::SHK_BITS));

    println!("\nraw keys (Sdk steps 3 and 4)");
    let ka = nike::raw_key(&ska, &pkb);
    let kb = nike::raw_key(&skb, &pka);
    let gap = (0..D)
        .map(|i| fq::to_centred(fq::sub(ka[i], kb[i])).unsigned_abs())
        .max()
        .unwrap();
    println!("  max |k_1 - k_2| over the {D} coefficients = {gap}");
    println!("  worst case 2 m d eta beta = {lhs}, radius q/2^(rho+1) = {rhs}  {}", ok(gap as u64 <= rhs));

    println!("\nlist reconciliation (fig:adv-list-rec)");
    let pa = listrec::plan(&ka);
    let pb = listrec::plan(&kb);
    println!("  hedged coefficients: {} (must equal log l = {LOG_ELL})", pa.hedged_indices().len());
    println!("  A's first entry     : {}", hex(&pa.first()));
    println!("  B's first entry     : {}", hex(&pb.first()));
    match pb.scan(&pa.first()) {
        Some(t) => println!("  A's key is entry {t} of B's list of {ELL}   ok"),
        None => println!("  A's key is NOT in B's list                 FAILS"),
    }
    match pa.scan(&pb.first()) {
        Some(t) => println!("  B's key is entry {t} of A's list of {ELL}   ok"),
        None => println!("  B's key is NOT in A's list                 FAILS"),
    }

    println!("\nrepeat over 32 independent key pairs");
    let mut bad = 0;
    let mut worst = 0u32;
    for t in 0..32u8 {
        let (s1, p1) = nike::gen(&pp, &[0x10 + t; 16]);
        let (s2, p2) = nike::gen(&pp, &[0x80 + t; 16]);
        let k1 = nike::raw_key(&s1, &p2);
        let k2 = nike::raw_key(&s2, &p1);
        let g = (0..D)
            .map(|i| fq::to_centred(fq::sub(k1[i], k2[i])).unsigned_abs())
            .max()
            .unwrap();
        worst = worst.max(g);
        let l1 = listrec::plan(&k1);
        let l2 = listrec::plan(&k2);
        if !l2.contains(&l1.first()) || !l1.contains(&l2.first()) {
            bad += 1;
        }
    }
    println!("  failures: {bad}/32  {}", ok(bad == 0));
    println!("  largest gap seen: {worst}");

    println!("\narithmetic spot checks");
    let r = poly::mul_ref(&[1u32; D], &{
        let mut p = [0u32; D];
        p[D - 1] = 1;
        p
    });
    println!("  (sum X^i) * X^(d-1) has constant coefficient {} = -1 mod q  {}",
        r[0], ok(r[0] == Q - 1));
    println!("  q * q mod q = {}  {}", fq::mul(Q - 1, Q - 1), ok(fq::mul(Q - 1, Q - 1) == 1));
}

fn ok(b: bool) -> &'static str {
    if b { "ok" } else { "FAILS" }
}

fn hex(k: &listrec::ShKey) -> String {
    k.iter().map(|b| format!("{b:02x}")).collect()
}
