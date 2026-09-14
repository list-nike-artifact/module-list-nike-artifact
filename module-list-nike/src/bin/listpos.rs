//! How much of the responder's list does it actually have to look at?
//!
//! Correctness only asks that the initiator's entry occur *somewhere* in the
//! responder's list, so `ELL` is a bound on what the responder pays and not its
//! expected cost.  What decides the cost is how many of the hedged coefficients
//! have to take their alternative value, and that is a property of the pair, not
//! of `ELL`.  This measures it over honest pairs; the paper quotes the figures
//! it prints.
use mlnike::listrec::Plan;
use mlnike::nike::{gen, sdk_plan, stp};
use mlnike::params::{ELL, LOG_ELL};

const TRIALS: u32 = 1000;

fn main() {
    let pp = stp(&[0x5au8; 16]); // the same parameters as bin/bench
    let mut weight = [0u32; LOG_ELL + 1];
    let mut looks: Vec<u64> = Vec::with_capacity(TRIALS as usize);
    let mut index: Vec<u64> = Vec::with_capacity(TRIALS as usize);

    for t in 0..TRIALS {
        let (mut sa, mut sb) = ([0u8; 16], [0u8; 16]);
        sa[0..4].copy_from_slice(&t.to_le_bytes());
        sb[0..4].copy_from_slice(&t.to_le_bytes());
        sa[8] = 1;
        sb[8] = 2;
        let (ska, pka) = gen(&pp, &sa);
        let (skb, pkb) = gen(&pp, &sb);
        let pa: Plan = sdk_plan(&ska, &pkb); // initiator: names the key
        let pb: Plan = sdk_plan(&skb, &pka); // responder: searches for it
        let target = pa.first();
        let Some(i) = pb.scan_all(&target) else {
            println!("trial {t}: the initiator's key is not in the responder's list");
            std::process::exit(1);
        };
        // Entry `i` deviates at the hedged ranks given by the set bits of `i`.
        weight[i.count_ones() as usize] += 1;
        index.push(i);
        looks.push(pb.scan(&target).expect("found by index but not by weight"));
    }

    let n = TRIALS as usize;
    let stats = |v: &mut Vec<u64>| {
        v.sort_unstable();
        let mean = v.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
        (v[n / 2], v[n * 99 / 100], v[n - 1], mean)
    };

    println!("honest pairs                          {n}");
    println!("list length l                          {ELL}  (2^{LOG_ELL})\n");

    println!("hedged coefficients that had to deviate");
    for (w, &c) in weight.iter().enumerate() {
        if c > 0 {
            println!("  {w:2} of {LOG_ELL}                             {c:5} pairs");
        }
    }
    let heavy = weight.iter().enumerate().filter(|&(w, _)| w > 3).map(|(_, &c)| c).sum::<u32>();
    println!("  more than 3, so beyond the prefix  {heavy:5} pairs\n");

    println!("{:<38} {:>10} {:>10} {:>10} {:>12}", "entries examined", "median", "p99", "max", "mean");
    let (m, p, x, a) = stats(&mut looks);
    println!("{:<38} {m:>10} {p:>10} {x:>10} {a:>12.1}", "  lightest first, as implemented");
    let (m, p, x, a) = stats(&mut index);
    println!("{:<38} {m:>10} {p:>10} {x:>10} {a:>12.1}", "  in index order, for comparison");
}
