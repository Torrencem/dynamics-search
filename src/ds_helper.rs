
use crate::util::*;
use std::collections::{HashSet, BTreeMap};
use std::hash::BuildHasherDefault;
use fnv::FnvHasher;
use arrayvec::ArrayVec;
use crate::math::*;

type FNVHashSet = HashSet<usize, BuildHasherDefault<FnvHasher>>;

// In general: for a polynomial in Q, find the possible periods
// greater than goal
#[allow(unused)]
pub fn possible_periods_search(f: PolynomialInQ, goal: usize) -> Option<FNVHashSet> {
    let mut res = FNVHashSet::default();
    let mut first = true;
    for p in 2..=100 {
        if prime(p) && f.has_good_reduction(p) {
            if first {
                res = fast_possible_periods(f.do_reduction(p));
                first = false;
            } else {
                let pers = fast_possible_periods(f.do_reduction(p));
                res = res.intersection(&pers).map(|&x| x).collect();
            }
            // Check if our set contains anything
            // large enough to be interesting
            let mut found = false;
            for possible in &res {
                if *possible > goal {
                    found = true;
                    break;
                }
            }
            if !found {
                return None;
            }
        }
    }

    // Remove everything not in the goal
    let not_interesting: FNVHashSet = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res)
}

// Find the possible periods of a polynomial
// in Q(w)
pub fn possible_periods_search_qw(f: PolynomialInQw, goal: usize) -> Option<FNVHashSet> {
    let mut res = FNVHashSet::default();
    let mut first = true;
    for p in 2..=300 {
        if prime(p as usize) && has_qw_homomorphism(p) {
            if first {
                let (red1, red2) = f.reductions(p);
                match (red1, red2) {
                    (None, None) => continue,
                    (Some(r1), None) => res = fast_possible_periods(r1),
                    (None, Some(r2)) => res = fast_possible_periods(r2),
                    (Some(r1), Some(r2)) => {
                        res = fast_possible_periods(r1);
                        res = res.intersection(&fast_possible_periods(r2)).map(|&x| x).collect();
                    }
                }
                first = false;
            } else {
                let (red1, red2) = f.reductions(p);
                if let Some(r1) = red1 {
                    res = res.intersection(&fast_possible_periods(r1)).map(|&x| x).collect();
                }
                if let Some(r2) = red2 {
                    res = res.intersection(&fast_possible_periods(r2)).map(|&x| x).collect();
                }
                // Check if our set contains anything
                // large enough to be interesting
                let mut found = false;
                for possible in &res {
                    if *possible > goal {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return None;
                }
            }
        }
    }

    // Remove everything not in the goal
    let not_interesting: FNVHashSet = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res)
}

pub fn fast_possible_periods(f: Polynomial) -> FNVHashSet {
    let p = f.p_mod.unwrap();

    // We don't need the point at infinity, so we
    // can use an array of size p instead of p^2
    //    (also, this way, hash_point = id)
    let mut point_table = vec![(0, 0); p as usize];
    let mut index = 1;
    let mut periods = FNVHashSet::default();

    for p_start in 0..p {
        let mut P = p_start;
        let mut hash_p = P as usize;
        if point_table[hash_p].1 == 0 {
            let startindex = index;
            while point_table[hash_p].1 == 0 {
                point_table[hash_p].1 = index;
                let Q = f.eval(P);
                let hash_q = Q as usize;
                point_table[hash_p].0 = hash_q;
                P = Q;
                hash_p = hash_q;
                index += 1;
            }

            if point_table[hash_p].1 >= startindex {
                let period = index - point_table[hash_p].1;
                periods.insert(period);
                let charpoly_constant = f.multiplier(period, P);
                if charpoly_constant == 0 {
                    continue; // Exclude 0
                }
                // lrorder is both lorder and rorder from sage
                let lrorder = multiplicative_order(charpoly_constant, p);
                
                let r = lrorder as usize;
                periods.insert(period * r);
                if p == 2 || p == 3 { // Over Q or Q(w), we need to consider e=1 for p=2 and p=3
                    periods.insert(period * r * (p as usize));
                }
            }
        }
    }

    periods
}

pub fn z4c_possible_periods_search(c: Rational, goal: usize) -> Option<FNVHashSet> {
    let mut res = FNVHashSet::default();
    let mut first = true;
    for p in 2..=100 {
        if prime(p) && c.denom % p as i64 != 0 {
            if first {
                res = z4_table_possible_periods(p, c.reduce(p)).clone();
                first = false;
            } else {
                let pers = z4_table_possible_periods(p, c.reduce(p));
                res = res.intersection(&pers).map(|&x| x).collect();
            }
            // Check if our set contains anything
            // large enough to be interesting
            let mut found = false;
            for possible in &res {
                if *possible > goal {
                    found = true;
                    break;
                }
            }
            if !found {
                return None;
            }
        }
    }

    // Remove everything not in the goal
    let not_interesting: FNVHashSet = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res.clone())
}

pub fn z3c_possible_periods_search(c: QwElement, goal: usize) -> Option<FNVHashSet> {
    let mut res = FNVHashSet::default();
    let mut first = true;
    for p in 2..=100 {
        if prime(p) && has_qw_homomorphism(p as i64) {
            if first {
                let (red1, red2) = c.reductions(p as i64);
                match (red1, red2) {
                    (None, None) => continue,
                    (Some(r1), None) => res = z3_table_possible_periods(p, r1 as usize).clone(),
                    (None, Some(r2)) => res = z3_table_possible_periods(p, r2 as usize).clone(),
                    (Some(r1), Some(r2)) => {
                        res = z3_table_possible_periods(p, r1 as usize).clone();
                        res = res.intersection(&z3_table_possible_periods(p, r2 as usize)).map(|&x| x).collect();
                    }
                }
                first = false;
            } else {
                let (red1, red2) = c.reductions(p as i64);
                if let Some(r1) = red1 {
                    res = res.intersection(&z3_table_possible_periods(p, r1 as usize)).map(|&x| x).collect();
                }
                if let Some(r2) = red2 {
                    res = res.intersection(&z3_table_possible_periods(p, r2 as usize)).map(|&x| x).collect();
                }
                // Check if our set contains anything
                // large enough to be interesting
                let mut found = false;
                for possible in &res {
                    if *possible > goal {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return None;
                }
            }
        }
    }

    // Remove everything not in the goal
    let not_interesting: FNVHashSet = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res)
}

pub fn z4_table_possible_periods(p: usize, c: usize) -> &'static FNVHashSet {
    &Z4_TABLE[p - 2][c]
}

pub fn z3_table_possible_periods(p: usize, c: usize) -> &'static FNVHashSet {
    &Z3_TABLE[p - 2][c]
}

lazy_static! {
    static ref Z4_TABLE: Vec<Vec<FNVHashSet>> = {
        let mut res = Vec::with_capacity(101);
        for p in 2..=100 {
            let mut interm = Vec::with_capacity(101);
            if !prime(p) {
                res.push(interm);
                continue;
            }
            for c in 0..p {
                let fc = Polynomial::new(
                    vec![1, 0, 0, 0, c as i64], Some(p as i64)
                );
                let res2 = fast_possible_periods(fc);
                interm.push(res2);
            }
            res.push(interm);
        }
        res
    };

    static ref Z3_TABLE: Vec<Vec<FNVHashSet>> = {
        let mut res = Vec::with_capacity(101);
        for p in 2..=100 {
            let mut interm = Vec::with_capacity(101);
            if !prime(p) {
                res.push(interm);
                continue;
            }
            for c in 0..p {
                let fc = Polynomial::new(
                    vec![1, 0, 0, c as i64], Some(p as i64)
                );
                let res2 = fast_possible_periods(fc);
                interm.push(res2);
            }
            res.push(interm);
        }
        res
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduction_example() {
        let num = EisensteinInteger::new(-8, 5);
        let den = EisensteinInteger::new(6, 7);
        let nd = QwElement::new(num, den);

        assert_eq!(nd.reductions(7), (Some(2), Some(5)));

        println!("{:?}", num.gcd(&den));
    }
}