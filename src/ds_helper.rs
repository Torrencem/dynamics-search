
use crate::util::*;
use std::collections::HashSet;
use crate::math::*;

pub fn possible_periods_search(f: PolynomialInQ, goal: usize) -> Option<HashSet<usize>> {
    let mut res = HashSet::new();
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
    let not_interesting: HashSet<_> = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res)
}

#[allow(unused)]
pub fn possible_periods(f: PolynomialInQ, prime_bound: usize) -> HashSet<usize> {
    let mut res = HashSet::new();
    let mut first = true;
    for p in 2..=prime_bound {
        if small_prime(p) && f.has_good_reduction(p) {
            if first {
                res = fast_possible_periods(f.do_reduction(p));
                first = false;
            } else {
                let pers = fast_possible_periods(f.do_reduction(p));
                res = res.intersection(&pers).map(|&x| x).collect();
            }
        }
    }

    res
}

pub fn fast_possible_periods(f: Polynomial) -> HashSet<usize> {
    let p = f.p_mod.unwrap();

    // We don't need the point at infinity, so we
    // can use an array of size p instead of p^2
    //    (also, this way, hash_point = id)
    let mut point_table = vec![(0, 0); p as usize];
    let mut index = 1;
    let mut periods = HashSet::new();

    for p_start in 0..p {
        let mut P = p_start;
        let mut hash_p = P as usize;
        // Unsafe is used for get_unchecked and get_unchecked_mut,
        // which give a small performance boost over array
        // index checking
        unsafe {
            if point_table.get_unchecked(hash_p).1 == 0 {
                let startindex = index;
                while point_table.get_unchecked(hash_p).1 == 0 {
                    let mut pt_hash_p = point_table.get_unchecked_mut(hash_p);
                    pt_hash_p.1 = index;
                    let Q = f.eval(P);
                    let hash_q = Q as usize;
                    pt_hash_p.0 = hash_q;
                    P = Q;
                    hash_p = hash_q;
                    index += 1;
                }

                let pt_hash_p = point_table.get_unchecked(hash_p);
                if pt_hash_p.1 >= startindex {
                    let period = index - pt_hash_p.1;
                    periods.insert(period);
                    let charpoly_constant = f.multiplier(period, P);
                    if charpoly_constant == 0 {
                        continue; // Exclude 0
                    }
                    // lrorder is both lorder and rorder from sage
                    let lrorder = c_multiplicative_order(charpoly_constant, p);
                    
                    let r = lrorder as usize;
                    periods.insert(period * r);
                    if p == 2 || p == 3 {
                        periods.insert(period * r * (p as usize));
                    }
                }
            }
        }
    }

    periods
}