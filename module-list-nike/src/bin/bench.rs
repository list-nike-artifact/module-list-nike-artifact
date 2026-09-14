//! Cycle counts for the three algorithms of `fig:mod-list-nike`, with `ListRec`
//! broken out inside `Sdk`, where it is step 5.
//!
//! Every number is a median over the stated number of independent calls, timed
//! with `Instant` on a pinned performance core and multiplied by the core clock.

use mlnike::cycles::{self, measure, Stat};
use mlnike::listrec;
use mlnike::nike;
use mlnike::params::*;

fn main() {
    let full = std::env::args().any(|a| a == "--full-list");
    // Cycles are times multiplied by the core clock.  We use the published
    // P-core frequency, and print what `calibrate_ghz` independently measures
    // on the same core as a check on it.
    const GHZ: f64 = 4.05;
    // The same number of independent calls behind every row.
    const RUNS: usize = 1000;
    let ghz = GHZ;
    println!("clock: {GHZ} GHz (Apple M3 P-core), calibration on this core: {:.3} GHz\n",
        cycles::calibrate_ghz());

    let pp = nike::stp(&[0x5au8; 16]);
    let (ska, _) = nike::gen(&pp, &[1u8; 16]);
    let (_, pkb) = nike::gen(&pp, &[2u8; 16]);
    let ka = nike::raw_key(&ska, &pkb);
    let plan = listrec::plan(&ka);

    let stp = measure("Stp", ghz, RUNS, || nike::stp(&[0x5au8; 16]));
    let gen = measure("Gen", ghz, RUNS, || nike::gen(&pp, &[3u8; 16]));
    let sdk = measure("Sdk", ghz, RUNS, || nike::sdk_plan(&ska, &pkb));
    let lr = measure("of which ListRec", ghz, RUNS, || listrec::plan(&ka));
    let walk = measure("of which ListRec", ghz, RUNS, || plan.scan_all(&[0xffu8; SHK_BYTES]));

    println!("{:<40} {:>7} {:>11} {:>11} {:>10}", "", "runs", "median cyc", "mean cyc", "us");
    line(&stp, ghz);
    line(&gen, ghz);

    println!("\nSdk for a party that needs only its own key.  ListRec is step 5 of Sdk;");
    println!("its output is the list, a compact description of all {ELL} entries");
    println!("that already holds entry 1, so a party that wants that entry stops here:");
    line(&sdk, ghz);
    line2(&lr, ghz);

    // The honest responder, over RUNS independent pairs rather than RUNS
    // repetitions of one pair, since what it pays depends on the pair.
    let pairs: Vec<(listrec::Plan, listrec::ShKey)> = (0..RUNS as u32)
        .map(|t| {
            let (mut sa, mut sb) = ([0u8; 16], [0u8; 16]);
            sa[0..4].copy_from_slice(&t.to_le_bytes());
            sb[0..4].copy_from_slice(&t.to_le_bytes());
            sa[8] = 1;
            sb[8] = 2;
            let (ska, pka) = nike::gen(&pp, &sa);
            let (skb, pkb) = nike::gen(&pp, &sb);
            (nike::sdk_plan(&skb, &pka), nike::sdk_plan(&ska, &pkb).first())
        })
        .collect();
    let mut i = 0usize;
    let honest = measure("Sdk", ghz, RUNS, || {
        let (p, target) = &pairs[i % RUNS];
        i += 1;
        p.scan(target).expect("peer entry not in the list")
    });
    let examined: Vec<u64> = pairs.iter().map(|(p, t)| p.scan(t).unwrap()).collect();
    let mut sorted = examined.clone();
    sorted.sort_unstable();

    println!("\nSdk for the same party on an honest pair, which does not have to walk");
    println!("the whole list: the entry it wants is among the first few, so it looks at");
    println!("{} of the {ELL} entries at the median and {} at the worst of {RUNS} pairs:",
        sorted[RUNS / 2], sorted[RUNS - 1]);
    line(&Stat { label: "Sdk".into(), runs: honest.runs,
        median: sdk.median + honest.median, mean: 0.0 }, ghz);
    line2(&Stat { label: "of which ListRec".into(), runs: honest.runs,
        median: lr.median + honest.median, mean: lr.mean + honest.mean }, ghz);

    println!("\nSdk for a party that must find the peer's entry in the list.  The");
    println!("remaining {ELL} - 1 entries are the rest of step 5, hence of Sdk:");
    let total = Stat { label: "Sdk".into(), runs: walk.runs, median: sdk.median + walk.median, mean: 0.0 };
    line(&total, ghz);
    line2(&Stat { label: "of which ListRec".into(), runs: walk.runs,
        median: lr.median + walk.median, mean: lr.mean + walk.mean }, ghz);
    println!("  {:<62} {:>11.1}", "cycles per list entry", walk.median / ELL as f64);
    if full {
        let mut out = vec![0u8; ELL as usize * SHK_BYTES];
        let w = measure("of which ListRec, list written out", ghz, 3, || plan.enumerate_into(&mut out));
        line2(&w, ghz);
        println!("  {:<62} {:>11.1}", "cycles per list entry", w.median / ELL as f64);
    } else {
        println!("  (pass --full-list to also time writing all {} MiB out)",
            ELL as usize * SHK_BYTES / (1 << 20));
    }

    println!("\nWhere Gen and Sdk spend their time (schoolbook products in R_q, d = {D}):");
    println!("  Gen: {} products, {:.0} cycles each, {} coefficient products at {:.2}/cycle",
        2 * N * M, gen.per(2 * N * M), 2 * N * M * D * D,
        (2 * N * M * D * D) as f64 / gen.median);
    let k_stat = measure("k_stat", ghz, RUNS, || nike::k_stat(&ska, &pkb));
    let k_comp = measure("k_comp", ghz, RUNS, || nike::k_comp(&ska, &pkb));
    println!("  Sdk: k_stat = s^T x, {N} products, {:.0} cycles each, {} total",
        k_stat.per(N), fmt(k_stat.median));
    println!("       k_comp = b^T u, {M} products, {:.0} cycles each, {} total",
        k_comp.per(M), fmt(k_comp.median));
}

/// One measured row.
fn line(s: &Stat, ghz: f64) {
    let mean = if s.mean > 0.0 { fmt(s.mean) } else { String::new() };
    println!("  {:<38} {:>7} {:>11} {:>11} {:>10.1}", s.label, s.runs, fmt(s.median), mean, s.micros(ghz));
}

/// A row for a subroutine of the row above it.
fn line2(s: &Stat, ghz: f64) {
    println!("    {:<36} {:>7} {:>11} {:>11} {:>10.1}", s.label, s.runs, fmt(s.median), fmt(s.mean), s.micros(ghz));
}

fn fmt(c: f64) -> String {
    if c >= 1e6 {
        format!("{:.2}M", c / 1e6)
    } else if c >= 1e4 {
        format!("{:.1}k", c / 1e3)
    } else {
        format!("{c:.0}")
    }
}
