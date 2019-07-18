
use crate::util::*;
use std::collections::{HashSet, HashMap};
use crate::math::*;

// In general: for a polynomial in Q, find the possible periods
// greater than goal
#[allow(unused)]
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

// Find the possible periods of a polynomial
// in Q(w)
pub fn possible_periods_search_qw(f: PolynomialInQw, goal: usize) -> Option<HashSet<usize>> {
    let mut res = HashSet::new();
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
    let not_interesting: HashSet<_> = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res)
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
                if p == 2 || p == 3 {
                    periods.insert(period * r * (p as usize));
                }
            }
        }
    }

    periods
}

pub fn z4c_possible_periods_search(c: Rational, goal: usize) -> Option<HashSet<usize>> {
    let mut res = HashSet::new();
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
    let not_interesting: HashSet<_> = (0..=goal).collect();
    res = res.difference(&not_interesting).map(|&x| x).collect();

    Some(res.clone())
}

pub fn z4_table_possible_periods(p: usize, c: usize) -> &'static HashSet<usize> {
    Z4_TABLE.get(&p).unwrap().get(&c).unwrap()
}

lazy_static! {
    static ref Z4_TABLE: HashMap<usize, HashMap<usize, HashSet<usize>>> = {
        let mut res = HashMap::new();
        for p in 2..=100 {
            if !prime(p) {
                continue;
            }
            let mut at_p = HashMap::new();
            for c in 0..p {
                let fc = Polynomial::new(
                    vec![1, 0, 0, 0, c as i64], Some(p as i64)
                );
                let res = fast_possible_periods(fc);
                at_p.insert(c, res);
            }
            res.insert(p, at_p);
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