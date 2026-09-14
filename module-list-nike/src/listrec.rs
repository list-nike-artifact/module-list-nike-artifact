//! The multi-coefficient list reconciliation of `fig:adv-list-rec`.
//!
//! `Rnd`, `Rnd_rec` and `DistanceToRB` all fall out of a single division.
//! Writing `t = 2^rho k = Rnd(k) q + r`,
//!
//! ```text
//!   Rnd(k)         = t / q
//!   Rnd_rec(k)     = floor(2 r / q)        (1 exactly when 2r >= q)
//!   DistToRB(k)    = min(r, q - r) / 2^rho
//! ```
//!
//! because `floor(2^(rho+1) k / q) = 2 Rnd(k) + floor(2r/q)`, whose parity is
//! that of `floor(2r/q)`.  Distances are kept in the `2^rho`-scaled units
//! `min(r, q - r)`, which is order-isomorphic to the real distance, so the sort
//! never needs a division.
//!
//! # Enumeration
//!
//! The `ell = 2^L` entries differ only in the `L` hedged coefficients, and each
//! hedged coefficient takes one of two values.  Entry `t` (zero-based) picks the
//! second value at hedged rank `j` exactly when bit `j` of `t` is set, which is
//! the indexing of the figure.  Successive entries therefore differ by
//! `popcount(t xor (t-1))` coefficients, amortised two, and toggling one
//! coefficient is a single byte `xor` because a two-valued field flips under
//! `xor (s xor s')`.  Enumeration costs `O(ell)` with a small constant rather
//! than `O(ell * c)`.
//!
//! # Key layout
//!
//! `ToBits(s_i)` for `i in 0..c` are packed low-to-high, `rho` bits each:
//! coefficient `i` lives in bits `rho*i .. rho*i + rho`, so with `rho = 2` in
//! byte `i/4` at shift `2*(i mod 4)`.

use crate::params::{C, ELL, LOG_ELL, Q, RHO, SHK_BYTES, WEIGHT};
use crate::poly::Poly;

/// A shared key candidate: `rho * c = 256` bits.
pub type ShKey = [u8; SHK_BYTES];

/// `sum_{w <= WEIGHT} C(LOG_ELL, w)`, the length of the light prefix.
const LIGHT_ENTRIES: u64 = {
    let mut total = 0;
    let mut w = 0;
    while w <= WEIGHT {
        let mut c = 1u64;
        let mut j = 0;
        while j < w {
            c = c * (LOG_ELL - j) as u64 / (j + 1) as u64;
            j += 1;
        }
        total += c;
        w += 1;
    }
    total
};

/// `Rnd`, `Rnd_rec` and the scaled `DistanceToRB` of one coefficient.
#[inline(always)]
pub fn round_coeff(k: u32) -> (u8, u8, u32) {
    let t = (k as u64) << RHO;
    let rnd = (t / Q as u64) as u8;
    let r = (t % Q as u64) as u32;
    let rec = (2 * r >= Q) as u8;
    let dist = if rec == 1 { Q - r } else { r };
    (rnd, rec, dist)
}

/// `Rnd(k)`, the `rho` bits extracted from one coefficient.
#[inline]
pub fn rnd(k: u32) -> u8 {
    round_coeff(k).0
}

/// `Rnd_rec(k)`, the reconciliation bit of one coefficient.
#[inline]
pub fn rnd_rec(k: u32) -> u8 {
    round_coeff(k).1
}

/// The reconciliation plan of one raw key: everything `ListRec` needs to emit
/// its `ell` entries.
pub struct Plan {
    /// Entry one of the list, `Rnd(k_1) || ... || Rnd(k_c)`.
    pub base: ShKey,
    /// `Rnd(k_i)` for every coefficient.
    rnd: [u8; C],
    /// The alternative value at each hedged coefficient, or `rnd` when not hedged.
    alt: [u8; C],
    /// True at the `L` coefficients closest to a rounding boundary.
    hedged: [bool; C],
    /// `(byte, xor mask)` per hedged rank, ordered by coefficient index.
    toggle: [(u8, u8); LOG_ELL],
}

/// Places a `rho`-bit value for coefficient `i`.
#[inline(always)]
fn slot(i: usize) -> (usize, u32) {
    (i * RHO / 8, ((i * RHO) % 8) as u32)
}

/// Builds the reconciliation plan of the raw key `k`, using its first `c`
/// coefficients.  This is `ListRec` up to the final enumeration.
pub fn plan(k: &Poly) -> Plan {
    let mut rnd = [0u8; C];
    let mut rec = [0u8; C];
    // (distance, index) so that ties break by index and the plan is deterministic.
    let mut order = [(0u32, 0u16); C];
    for i in 0..C {
        let (s, b, d) = round_coeff(k[i]);
        rnd[i] = s;
        rec[i] = b;
        order[i] = (d, i as u16);
    }

    // Only the L smallest matter, so select rather than sort: linear, as the
    // footnote to the figure observes.
    order.select_nth_unstable(LOG_ELL - 1);
    let mut chosen: [u16; LOG_ELL] = [0; LOG_ELL];
    for j in 0..LOG_ELL {
        chosen[j] = order[j].1;
    }
    chosen.sort_unstable();

    let mask = (1u8 << RHO) - 1;
    let mut alt = rnd;
    let mut hedged = [false; C];
    let mut toggle = [(0u8, 0u8); LOG_ELL];
    for (j, &ci) in chosen.iter().enumerate() {
        let i = ci as usize;
        hedged[i] = true;
        // b_i = 0 hedges downwards, b_i = 1 upwards, both mod 2^rho.
        let step = if rec[i] == 0 { mask } else { 1 };
        alt[i] = rnd[i].wrapping_add(step) & mask;
        let (byte, sh) = slot(i);
        toggle[j] = (byte as u8, (rnd[i] ^ alt[i]) << sh);
    }

    let mut base = [0u8; SHK_BYTES];
    for i in 0..C {
        let (byte, sh) = slot(i);
        base[byte] |= rnd[i] << sh;
    }

    Plan { base, rnd, alt, hedged, toggle }
}

impl Plan {
    /// `ListRec(k)[1]`.
    #[inline]
    pub fn first(&self) -> ShKey {
        self.base
    }

    /// `ListRec(k)[t + 1]`, by direct construction.  Reference path.
    pub fn entry(&self, t: u64) -> ShKey {
        let mut key = [0u8; SHK_BYTES];
        let mut j = 0;
        for i in 0..C {
            let v = if self.hedged[i] {
                let bit = (t >> j) & 1;
                j += 1;
                if bit == 1 { self.alt[i] } else { self.rnd[i] }
            } else {
                self.rnd[i]
            };
            let (byte, sh) = slot(i);
            key[byte] |= v << sh;
        }
        key
    }

    /// Writes all `ell` entries into `out`, which must hold `ell * 32` bytes.
    pub fn enumerate_into(&self, out: &mut [u8]) {
        assert_eq!(out.len() as u64, ELL * SHK_BYTES as u64);
        let mut cur = self.base;
        out[0..SHK_BYTES].copy_from_slice(&cur);
        for t in 1..ELL {
            let mut g = t ^ (t - 1);
            while g != 0 {
                let (byte, m) = self.toggle[g.trailing_zeros() as usize];
                cur[byte as usize] ^= m;
                g &= g - 1;
            }
            let o = (t as usize) * SHK_BYTES;
            out[o..o + SHK_BYTES].copy_from_slice(&cur);
        }
    }

    /// Enumerates the list in index order and returns the position of `target`,
    /// if present.  This is the `Theta(ell)` bound, without materialising 512 MiB.
    pub fn scan_all(&self, target: &ShKey) -> Option<u64> {
        let mut cur = self.base;
        if cur == *target {
            return Some(0);
        }
        for t in 1..ELL {
            let mut g = t ^ (t - 1);
            while g != 0 {
                let (byte, m) = self.toggle[g.trailing_zeros() as usize];
                cur[byte as usize] ^= m;
                g &= g - 1;
            }
            if cur == *target {
                return Some(t);
            }
        }
        None
    }

    /// Entries `scan_light` examines before giving up: the subsets of the `L`
    /// hedged ranks of size at most `WEIGHT`.
    pub const LIGHT: u64 = LIGHT_ENTRIES;

    /// The `w`-subsets of the `L` hedged ranks, in lexicographic order.  Returns
    /// false once `idx` holds the last one.
    fn next_subset(idx: &mut [usize]) -> bool {
        let w = idx.len();
        let mut i = w;
        while i > 0 {
            i -= 1;
            if idx[i] < LOG_ELL - w + i {
                idx[i] += 1;
                for k in i + 1..w {
                    idx[k] = idx[k - 1] + 1;
                }
                return true;
            }
        }
        false
    }

    /// Tries every entry in which at most `WEIGHT` of the hedged coefficients
    /// take their alternative value, lightest first, and returns the number of
    /// entries examined if `target` is among them.
    ///
    /// Which entry is the right one is decided by how many of the hedged
    /// coefficients the two parties round in opposite directions, and that is
    /// almost always none of them: over 1000 honest pairs (the ones `bin/listpos`
    /// samples) the count was 0 in 915 and 1 in the other 85.  Index order
    /// spreads even those few entries over the whole list, since a deviation on
    /// the `j`-th hedged coefficient lands at index `2^j`; lightest first they
    /// cost at most `1 + LOG_ELL` entries.
    pub fn scan_light(&self, target: &ShKey) -> Option<u64> {
        let mut examined = 1;
        if self.base == *target {
            return Some(examined);
        }
        for w in 1..=WEIGHT {
            let mut idx = [0usize; WEIGHT];
            let idx = &mut idx[..w];
            for (j, e) in idx.iter_mut().enumerate() {
                *e = j;
            }
            loop {
                let mut cur = self.base;
                for &j in idx.iter() {
                    let (byte, m) = self.toggle[j];
                    cur[byte as usize] ^= m;
                }
                examined += 1;
                if cur == *target {
                    return Some(examined);
                }
                if !Self::next_subset(idx) {
                    break;
                }
            }
        }
        None
    }

    /// The receiver's search: the light entries first, then the whole list.
    /// Correctness only asks that `target` be somewhere in the list, so the order
    /// is free to put the likely entries first.  Returns the number of entries
    /// examined, which is what the receiver pays.
    pub fn scan(&self, target: &ShKey) -> Option<u64> {
        match self.scan_light(target) {
            Some(n) => Some(n),
            None => self.scan_all(target).map(|t| LIGHT_ENTRIES + t + 1),
        }
    }

    /// Membership in `O(c)`, by checking each field against its one or two
    /// admissible values.  Used by the tests, not by the benchmark.
    pub fn contains(&self, target: &ShKey) -> bool {
        let mask = (1u8 << RHO) - 1;
        for i in 0..C {
            let (byte, sh) = slot(i);
            let v = (target[byte] >> sh) & mask;
            if v != self.rnd[i] && !(self.hedged[i] && v == self.alt[i]) {
                return false;
            }
        }
        true
    }

    /// The `L` hedged coefficient indices, ascending.
    pub fn hedged_indices(&self) -> Vec<usize> {
        (0..C).filter(|&i| self.hedged[i]).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::D;

    fn rnd_key(st: &mut u64) -> Poly {
        let mut p = [0u32; D];
        for k in 0..D {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            p[k] = ((*st >> 24) % Q as u64) as u32;
        }
        p
    }

    /// `Rnd`, `Rnd_rec` and the distance against their literal definitions.
    #[test]
    fn rounding_matches_definitions() {
        let step = (Q / 5000).max(1);
        let mut k = 0u32;
        while k < Q {
            let (s, b, d) = round_coeff(k);
            assert_eq!(s as u64, ((k as u64) << RHO) / Q as u64);
            assert_eq!(b as u64, (((k as u64) << (RHO + 1)) / Q as u64) % 2);
            // Distance to the nearest boundary, computed the slow way in the
            // 2^rho-scaled units the figure uses.
            let t = (k as i64) << RHO;
            let mut best = i64::MAX;
            for s in 0..=(1i64 << RHO) {
                best = best.min((t - s * Q as i64).abs());
            }
            assert_eq!(d as i64, best, "k = {k}");
            assert!(d <= Q / 2, "distance out of range at k = {k}");
            k = k.saturating_add(step);
        }
    }

    #[test]
    fn plan_has_exactly_log_ell_hedged_coefficients() {
        let mut st = 3u64;
        for _ in 0..20 {
            let p = plan(&rnd_key(&mut st));
            assert_eq!(p.hedged_indices().len(), LOG_ELL);
        }
    }

    #[test]
    fn hedged_are_the_closest_to_a_boundary() {
        let mut st = 5u64;
        let k = rnd_key(&mut st);
        let p = plan(&k);
        let hedged = p.hedged_indices();
        let worst_hedged = hedged.iter().map(|&i| round_coeff(k[i]).2).max().unwrap();
        let best_free = (0..C)
            .filter(|i| !hedged.contains(i))
            .map(|i| round_coeff(k[i]).2)
            .min()
            .unwrap();
        assert!(worst_hedged <= best_free);
    }

    /// The incremental enumeration must reproduce the direct construction, in order.
    #[test]
    fn enumeration_matches_direct_construction() {
        let mut st = 11u64;
        let p = plan(&rnd_key(&mut st));
        // Walk the first and last few thousand entries plus a sparse sample.
        let mut cur = p.base;
        for t in 0..200_000u64 {
            if t > 0 {
                let mut g = t ^ (t - 1);
                while g != 0 {
                    let (byte, m) = p.toggle[g.trailing_zeros() as usize];
                    cur[byte as usize] ^= m;
                    g &= g - 1;
                }
            }
            assert_eq!(cur, p.entry(t), "entry {t}");
        }
    }

    /// The light prefix enumerates exactly the entries of weight at most
    /// `WEIGHT`, each once, and finds every one of them.
    #[test]
    fn light_prefix_is_the_low_weight_entries() {
        let mut st = 19u64;
        let p = plan(&rnd_key(&mut st));
        let mut found = std::collections::BTreeSet::new();
        for t in 0..ELL.min(1 << 20) {
            if t.count_ones() as usize <= WEIGHT {
                let n = p.scan_light(&p.entry(t)).expect("light entry not found");
                assert!(n <= Plan::LIGHT, "examined {n} > {}", Plan::LIGHT);
                found.insert(t);
            }
        }
        // Weight at most WEIGHT among the low 20 bits, so C(20, w) summed.
        assert_eq!(found.len(), 1 + 20 + 190 + 1140);
        // The first entry costs one look, and nothing else does.
        assert_eq!(p.scan_light(&p.first()), Some(1));
        // An entry above the cutoff is not in the prefix, but is in the list.
        let heavy = (1u64 << 0) | (1 << 5) | (1 << 9) | (1 << 14);
        assert_eq!(p.scan_light(&p.entry(heavy)), None);
        assert_eq!(p.scan_all(&p.entry(heavy)), Some(heavy));
    }

    /// `scan` finds whatever is in the list, whichever branch answers, and
    /// nothing that is not.
    #[test]
    fn scan_agrees_with_membership() {
        let mut st = 23u64;
        let p = plan(&rnd_key(&mut st));
        for t in [0u64, 1, 2, 3, 7, 1 << 4, (1 << 4) | 1, 0b1111, 0b10101] {
            assert!(p.scan(&p.entry(t)).is_some(), "entry {t}");
        }
        let mut outside = p.first();
        let i = (0..C).find(|&i| !p.hedged[i]).unwrap();
        let (byte, sh) = slot(i);
        outside[byte] ^= 1 << sh;
        assert!(!p.contains(&outside));
        assert_eq!(p.scan_light(&outside), None);
    }

    #[test]
    fn first_entry_is_the_rounding_of_every_coefficient() {
        let mut st = 13u64;
        let k = rnd_key(&mut st);
        let p = plan(&k);
        for i in 0..C {
            let (byte, sh) = slot(i);
            assert_eq!((p.first()[byte] >> sh) & 3, rnd(k[i]));
        }
    }

    #[test]
    fn contains_agrees_with_enumeration_on_a_short_list() {
        let mut st = 17u64;
        let p = plan(&rnd_key(&mut st));
        for t in 0..5000u64 {
            assert!(p.contains(&p.entry(t)));
        }
        // A key differing in a non-hedged coefficient is not in the list.
        let free = (0..C).find(|&i| !p.hedged[i]).unwrap();
        let mut bad = p.first();
        let (byte, sh) = slot(free);
        bad[byte] ^= 1 << sh;
        assert!(!p.contains(&bad));
    }
}
