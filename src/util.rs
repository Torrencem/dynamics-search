use std::fmt;

use crate::math::*;

use num_rational::Rational64;
use num_complex::Complex32;

use std::f32::consts::PI;

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
        let mut exp = self.coeffs.len() - 1;
        let mut first = false;
        for c in &self.coeffs {
            if first {
                write!(f, "{}x^{}", c, exp)?;
                first = false;
            } else {
                if exp == 0 {
                    write!(f, " + {}", c)?;
                } else if exp == 1 {
                    write!(f, " + {}x", c)?;
                } else {
                    write!(f, " + {}x^{}", c, exp)?;
                }
            }
            exp -= 1;
        }
        Ok(())
    }
}

#[derive(PartialEq, Debug, Copy, Clone)]
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

    pub fn reduce(&self, p: usize) -> usize {
        ((self.numer % (p as i64)) * mod_inverse(self.denom % (p as i64), p as i64)).rem_euclid(p as i64) as usize
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denom == 1 {
            write!(f, "{}", self.numer)
        } else {
            write!(f, "{}/{}", self.numer, self.denom)
        }
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

impl fmt::Display for PolynomialInQ {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut exp = self.coeffs.len() - 1;
        let mut first = false;
        for c in &self.coeffs {
            if first {
                write!(f, "{}x^{}", c, exp)?;
                first = false;
            } else {
                if exp == 0 {
                    write!(f, " + {}", c)?;
                } else if exp == 1 {
                    write!(f, " + {}x", c)?;
                } else {
                    write!(f, " + {}x^{}", c, exp)?;
                }
            }
            exp -= 1;
        }
        Ok(())
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct EisensteinInteger {
    a: i64,
    b: i64
}

impl EisensteinInteger {
    pub fn new(a: i64, b: i64) -> EisensteinInteger {
        EisensteinInteger {a, b}
    }

    // Evaluate the natural homomorphisms from Z[w]
    // to F_p sending 1 to 1 and w to sqrt(-3) in F_p
    pub fn reductions(&self, p: i64) -> (i64, i64) {
        debug_assert!(p > 2);
        let d2 = mod_inverse(2, p);
        let (s1, s2) = cipolla(p - 3, p).unwrap();
        let (w1, w2) = ((-1 + s1)*d2, (-1 + s2)*d2);
        ((self.a + w1 * self.b).rem_euclid(p), (self.a + w2 * self.b).rem_euclid(p))
    }

    pub fn one() -> EisensteinInteger {
        EisensteinInteger {a:1, b:0}
    }

    pub fn zero() -> EisensteinInteger {
        EisensteinInteger {a:0, b:0}
    }

    pub fn is_unit(&self) -> bool {
        if self.a.abs() == 1 && self.b == 0 {
            true
        } else if self.b.abs() == 1 && self.a == 0 {
            true
        } else if self.a == -1 && self.b == -1 {
            true
        } else {
            false
        }
    }

    pub fn is_zero(&self) -> bool {
        self.a == 0 && self.b == 0
    }

    pub fn product(&self, other: EisensteinInteger) -> EisensteinInteger {
        EisensteinInteger { a: self.a * other.a - self.b * other.b, b: self.a * other.b + self.b * other.a - self.b * other.b }
    }

    pub fn conjugate(&self) -> EisensteinInteger {
        EisensteinInteger { a: self.a - self.b, b: -self.b }
    }

    // Returns |x|^2 in C
    pub fn norm_sq(&self) -> i64 {
        self.a * self.a - self.a * self.b + self.b * self.b
    }

    pub fn difference(&self, other: &EisensteinInteger) -> EisensteinInteger {
        EisensteinInteger {a: self.a - other.a, b: self.b - other.b}
    }

    pub fn division(&self, other: &EisensteinInteger) -> EisensteinInteger {
        // alpha / beta = 1/(|beta|^2)alpha beta_conjugate
        // in complex absolute value
        let num = self.product(other.conjugate());
        let normsq = other.norm_sq();
        let a = Rational64::new(num.a, normsq);
        let b = Rational64::new(num.b, normsq);
        // Round a and b to the nearest integer
        let a = a.round().to_integer();
        let b = b.round().to_integer();
        EisensteinInteger { a, b }
    }

    pub fn gcd(&self, other: &EisensteinInteger) -> EisensteinInteger {
        if self.norm_sq() < other.norm_sq() {
            return other.gcd(self);
        }
        assert!(!other.is_zero());
        let q = self.division(other);
        let remainder = self.difference(&other.product(q));
        if remainder.is_zero() {
            *other
        } else {
            other.gcd(&remainder)
        }
    }

    pub fn approx_coords(&self) -> Complex32 {
        let a = self.a as f32;
        let b = self.b as f32;
        Complex32::new(a - (b/2.0), b * (3.0f32).sqrt() / 2.0)
    }

    pub fn phase_angle(&self) -> f32 {
        let raw = self.approx_coords().to_polar().1;
        if raw < 0.0 {
            raw + 2.0*PI
        } else if raw >= 2.0*PI {
            raw - 2.0*PI
        } else {
            raw
        }
    }
}

impl fmt::Display for EisensteinInteger {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.b == 0 {
            write!(f, "{}", self.a)
        } else if self.a == 0 {
            write!(f, "{}*w", self.b)
        } else {
            write!(f, "{} + {}*w", self.a, self.b)
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct QwElement {
    numer: EisensteinInteger,
    denom: EisensteinInteger
}

impl QwElement {
    pub fn new(numer: EisensteinInteger, denom: EisensteinInteger) -> QwElement {
        QwElement { numer, denom }
    }
    // Evaluate the natural homomorphisms from Z[w]
    // to F_p sending 1 to 1 and w to sqrt(-3) in F_p
    // Returns None if either choice has bad reduction
    pub fn reductions(&self, p: i64) -> (Option<i64>, Option<i64>) {
        let (a1, a2) = self.numer.reductions(p);
        let (b1, b2) = self.denom.reductions(p);

        let r1 = {
            if b1 == 0 {
                None
            } else {
                Some((a1 * mod_inverse(b1, p)).rem_euclid(p))
            }
        };

        let r2 = {
            if b2 == 0 {
                None
            } else {
                Some((a2 * mod_inverse(b2, p)).rem_euclid(p))
            }
        };

        (r1, r2)
    }

    pub fn one() -> QwElement {
        QwElement::new(EisensteinInteger::one(), EisensteinInteger::one())
    }

    pub fn zero() -> QwElement {
        QwElement::new(EisensteinInteger::zero(), EisensteinInteger::one())
    }
    
    pub fn approx_coords(&self) -> Complex32 {
        self.numer.approx_coords() / self.denom.approx_coords()
    }

    pub fn phase_angle(&self) -> f32 {
        let raw = self.approx_coords().to_polar().1;
        if raw < 0.0 {
            raw + 2.0*PI
        } else if raw > 2.0*PI {
            raw - 2.0*PI
        } else {
            raw
        }
    }
}

impl fmt::Display for QwElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denom == EisensteinInteger::one() {
            write!(f, "{}", self.numer)
        } else if self.numer == EisensteinInteger::zero() {
            write!(f, "0")
        } else if self.numer.b == 0 {
            write!(f, "{}/({})", self.numer, self.denom)
        } else {
            write!(f, "({})/({})", self.numer, self.denom)
        }
    }
}

pub struct PolynomialInQw {
    pub coeffs: Vec<QwElement>,
}

impl PolynomialInQw {
    pub fn from(coeffs: Vec<QwElement>) -> PolynomialInQw {
        PolynomialInQw {coeffs}
    }

    pub fn reductions(&self, p: i64) -> (Option<Polynomial>, Option<Polynomial>) {
        let mut p1 = Vec::with_capacity(self.coeffs.len());
        let mut p2 = Vec::with_capacity(self.coeffs.len());

        let mut has_p1 = true;
        let mut has_p2 = true;

        for c in &self.coeffs {
            let (red1, red2) = c.reductions(p);
            if has_p1 {
                match red1 {
                    None => {
                        has_p1 = false;
                        continue
                    },
                    Some(red1) => {
                        p1.push(red1);
                    }
                }
            }
            if has_p2 {
                match red2 {
                    None => {
                        has_p2 = false;
                        continue
                    },
                    Some(red2) => {
                        p2.push(red2);
                    }
                }
            }
            if !(has_p1 || has_p2) {
                break;
            }
        }

        let p1 = {
            if !has_p1 {
                None
            } else {
                Some(Polynomial::new(p1, Some(p)))
            }
        };

        let p2 = {
            if !has_p2 {
                None
            } else {
                Some(Polynomial::new(p2, Some(p)))
            }
        };

        (p1, p2)
    }
}

impl fmt::Display for PolynomialInQw {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut exp = self.coeffs.len() - 1;
        let mut first = false;
        for c in &self.coeffs {
            if first {
                write!(f, "{}x^{}", c, exp)?;
                first = false;
            } else {
                if exp == 0 {
                    write!(f, " + {}", c)?;
                } else if exp == 1 {
                    write!(f, " + {}x", c)?;
                } else {
                    write!(f, " + {}x^{}", c, exp)?;
                }
            }
            exp -= 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn another_example() {
        let a = EisensteinInteger::new(2, 3);
        let b = EisensteinInteger::new(0, 1);
        println!("{:?}", a.division(&b));
    }
}