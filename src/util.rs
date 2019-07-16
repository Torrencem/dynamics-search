use std::fmt;

use crate::math::*;

#[derive(Debug)]
pub struct Polynomial {
    pub coeffs: Vec<i64>,
    pub p_mod: Option<i64>,
}

impl Polynomial {
    pub fn new(v: Vec<i64>, p_mod: Option<i64>) -> Polynomial {
        Polynomial { coeffs: v, p_mod: p_mod }
    }

    pub fn eval(&self, x: i64) -> i64 {
        let e = self.coeffs.len();
        let mut res = 0;
        for indx in 0..e-1 {
            res += self.coeffs[indx];
            if let Some(p) = self.p_mod {
                let r = res.checked_mul(x);
                match r {
                    None => {
                        res %= p;
                        res *= x;
                    },
                    Some(r) => {
                        res = r;
                    }
                }
            } else {
                res *= x;
            }
        }
        res += self.coeffs[e-1];
        if let Some(p) = self.p_mod {
            res.rem_euclid(p)
        } else {
            res
        }
    }

    pub fn derivative(&self) -> Polynomial {
        let l = self.coeffs.len();
        Polynomial::new((1..l)
            .map(|i| {
                self.coeffs[i-1] * ((l-i) as i64)
            }).collect(), self.p_mod)
    }

    pub fn n_orbit(&self, x: i64, n: usize) -> Vec<i64> {
        let mut res = Vec::with_capacity(n);
        let mut curr = x;
        res.push(curr);
        for _ in 0..n-1 {
            curr = self.eval(curr);
            res.push(curr);
        }
        res
    }

    pub fn multiplier(&self, period: usize, x: i64) -> i64 {
        let orbit = self.n_orbit(x, period);
        let s_der = self.derivative();
        if let Some(p_mod) = self.p_mod {
            orbit.into_iter().map(|a| {
                s_der.eval(a)
            }).fold(1i64, |prod, i| {
                let r = prod.checked_mul(i);
                match r {
                    None => {
                        (prod % p_mod) * i
                    },
                    Some(r) => {
                        r
                    }
                }
            }).rem_euclid(p_mod)
        } else {
            orbit.into_iter().map(|a| {
                s_der.eval(a)
            }).product()
        }
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = self.coeffs.clone().into_iter().map(|x| x.to_string() + " ").collect();
        write!(f, "{}", s.trim())
    }
}

pub struct Rational {
    pub numer: i64,
    pub denom: i64,
}

impl Rational {
    pub fn new(a: i64, b: i64) -> Rational {
        Rational {numer: a, denom: b}
    }

    pub fn zero() -> Rational {
        Rational {numer: 0, denom: 1}
    }

    pub fn one() -> Rational {
        Rational {numer: 1, denom: 1}
    }
}

pub struct PolynomialInQ {
    pub coeffs: Vec<Rational>,
}

impl PolynomialInQ {
    pub fn from(coeffs: Vec<Rational>) -> PolynomialInQ {
        PolynomialInQ {coeffs: coeffs}
    }

    pub fn has_good_reduction(&self, p: usize) -> bool {
        for c in &self.coeffs {
            if c.denom % (p as i64) == 0 {
                return false;
            }
        }
        return true;
    }

    pub fn do_reduction(&self, p: usize) -> Polynomial {
        let coeffs = self.coeffs.iter().map(|c| {
            (c.numer % (p as i64)) * mod_inverse(c.denom % (p as i64), p as i64) % (p as i64)
        }).collect();

        Polynomial::new(coeffs, Some(p as i64))
    }
}