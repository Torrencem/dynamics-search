use num_integer::{Integer};

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

    pub fn from(v: Vec<i64>) -> Polynomial {
        Polynomial { coeffs: v, p_mod: None }
    }

    pub fn slow_eval(&self, x: i64) -> i64 {
        let mut e = self.coeffs.len() as u32;
        let mut res = 0;

        for i in 0usize..self.coeffs.len() {
            let c = self.coeffs[i];
            e -= 1;
            if c != 0 {
                res += x.pow(e) * c;
                if let Some(p_mod) = self.p_mod {
                    res %= p_mod; // TODO: Use something faster
                }
            }
        }

        res
    }

    pub fn eval(&self, x: i64) -> i64 {
        let mut e = self.coeffs.len();
        let mut res = 0;
        for indx in 0..e-1 {
            res += self.coeffs[indx];
            res *= x;
            if let Some(p) = self.p_mod {
                res %= p;
            }
        }
        res += self.coeffs[e-1];
        res
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
            }).fold(1, |sum, i| (sum * i) % p_mod ).rem_euclid(p_mod)
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
    pub numer: i128,
    pub denom: i128,
}

impl Rational {
    pub fn new(a: i128, b: i128) -> Rational {
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

    pub fn has_reduction(&self, p: usize) -> bool {
        for c in &self.coeffs {
            if c.denom % (p as i128) == 0 {
                return false;
            }
        }
        return true;
    }

    pub fn do_reduction(&self, p: usize) -> Polynomial {
        let coeffs = self.coeffs.iter().map(|c| {
            ((c.numer % (p as i128)) as i64) * mod_inverse((c.denom % (p as i128)) as i64, p as i64)
        }).collect();

        Polynomial::new(coeffs, Some(p as i64))
    }
}