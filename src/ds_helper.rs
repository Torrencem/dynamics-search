
use crate::util::*;
use std::collections::HashSet;
use crate::math::*;

pub fn possible_periods_search(f: PolynomialInQ, goal: usize) -> Option<HashSet<usize>> {
    let mut res = HashSet::new();
    let mut first = true;
    for p in 2..=100 {
        if small_prime(p) && f.has_reduction(p) {
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

pub fn possible_periods(f: PolynomialInQ, goal: Option<usize>) -> (bool, HashSet<usize>) {
    let mut res = HashSet::new();
    let mut first = true;
    for p in 2..=100 {
        if small_prime(p) && f.has_reduction(p) {
            if first {
                res = fast_possible_periods(f.do_reduction(p));
                first = false;
            } else {
                let pers = fast_possible_periods(f.do_reduction(p));
                res = res.intersection(&pers).map(|&x| x).collect();
            }
            if let Some(max_cycle_len) = goal {
                // Check if our set contains anything
                // large enough to be interesting
                let mut found = false;
                for possible in &res {
                    if *possible > max_cycle_len {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return (false, res);
                }
            }
        }
    }

    // Remove everything not in the goal
    if let Some(max_cycle_len) = goal {
        let not_interesting: HashSet<_> = (0..=max_cycle_len).collect();
        res = res.difference(&not_interesting).map(|&x| x).collect();
    }

    (true, res)
}

pub fn fast_possible_periods(f: Polynomial) -> HashSet<usize> {
    let p = f.p_mod.unwrap();

    let mut point_table = vec![(0, 0); (p*p) as usize];
    let mut index = 2; // 2 to account for the point at infinity
    let mut periods = HashSet::new();
    periods.insert(1); // again for the point at infinity

    for p_start in 0..p {
        let mut P = p_start;
        let mut hash_p = (P + p) as usize;
        if point_table[hash_p].1 == 0 {
            let startindex = index;
            while point_table[hash_p].1 == 0 {
                point_table[hash_p].1 = index;
                let Q = f.eval(P).rem_euclid(p);
                let hash_q = Q as usize;
                point_table[hash_p].0 = hash_q;
                P = Q;
                hash_p = hash_q;
                index += 1;
            }

            if point_table[hash_p].1 >= startindex {
                let period = index - point_table[hash_p].1;
                periods.insert(period);
                let l = f.multiplier(period, P);
                let charpoly_constant = l;
                if charpoly_constant == 0 {
                    continue; // Exclude 0
                }
                // It's both
                let lrorder = multiplicative_order(charpoly_constant, p);
                
                let r = lrorder as usize;
                periods.insert(period * r);
                if p == 2 || p == 3 {
                    periods.insert(period * r * (p as usize));
                }
            }
        }
    }

    periods
}