#![feature(euclidean_division, test)]
mod math;
mod util;
mod ds_helper;

use math::*;
use util::*;
use ds_helper::*;
use rayon::prelude::*;

use num_integer::{Integer};

// use flint::fmpz_poly::*;
// use flint::arith::*;
// use flint::fmpz::Fmpz;
// use flint::traits::*;

extern crate num_integer;
extern crate rayon;
extern crate test;
//extern crate numeric_literals;

pub fn run_large_search() {
    (-100_000_000..=100_000_000i64).into_par_iter().for_each(|a| {
        for b in 1..=100i64 {
            let b = b*b*b*b;
            // Analytic bound (should be above or below gcd?)
            let flo:f32 = (a as f32) / (b as f32);
            if flo > -0.913942 {
                continue;
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                continue;
            }
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(a, b)]
            );

            let res = possible_periods_search(fc, 2);

            if let Some(possibilities) = res {
                println!("Manually check {}/{}\n(for: {:?})", a, b, possibilities);
            }
        }
    });

    println!("Completed Search!");
}

pub fn run_smaller_search() {
    (-1_000_000..=1_000_000i64).into_par_iter().for_each(|a| {
        for b in 1..=37i64 {
            let b = b*b*b*b;
            // Analytic bound (should be above or below gcd?)
            let flo:f32 = (a as f32) / (b as f32);
            if flo > -0.913942 {
                continue;
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                continue;
            }
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(a, b)]
            );

            let res = possible_periods_search(fc, 2);

            if let Some(possibilities) = res {
                println!("Manually check {}/{}\n(for: {:?})", a, b, possibilities);
            }
        }
    });

    println!("Completed search!");
}

fn main() {
    let mut output = "".to_string();
            (-3_000..=3_000i64).into_iter().for_each(|a| {
                for b in 1..=26i64 {
                    let b = b*b;
                    // Analytic bound (should be above or below gcd?)
                    let flo:f32 = (a as f32) / (b as f32);
                    if flo > -0.913942 {
                        continue;
                    }
                    let g = a.gcd(&b);
                    if g != 1 && g != -1 {
                        continue;
                    }
                    let fc = PolynomialInQ::from(
                        vec![Rational::one(), Rational::zero(), Rational::new(a, b)]
                    );

                    let res = possible_periods_search(fc, 2);

                    if let Some(possibilities) = res {
                        output += format!("Manually check {}/{}\n(for: {:?})\n", a, b, possibilities).as_ref();
                    }
                }
            });

            let correct_output = 
r"Manually check -2689/576
(for: {3})
Manually check -1849/576
(for: {3})
Manually check -421/144
(for: {3})
Manually check -301/144
(for: {3})
Manually check -29/16
(for: {3})
";
            println!("\n\n{}\n\n{}",output, correct_output);
            assert!(correct_output == output);

    println!("Completed search!");
}



#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn test_functionality() {
        let mut output = "".to_string();
        (-100_000..=100_000i64).into_iter().for_each(|a| {
            for b in 1..=60i64 {
                let b = b*b;
                // Analytic bound (should be above or below gcd?)
                let flo:f32 = (a as f32) / (b as f32);
                if flo > -0.913942 {
                    continue;
                }
                let g = a.gcd(&b);
                if g != 1 && g != -1 {
                    continue;
                }
                let fc = PolynomialInQ::from(
                    vec![Rational::one(), Rational::zero(), Rational::new(a, b)]
                );

                let res = possible_periods_search(fc, 2);

                if let Some(possibilities) = res {
                    output += format!("Manually check {}/{}\n(for: {:?})\n", a, b, possibilities).as_ref();
                }
            }
        });

        let correct_output = 
r"Manually check -34861/3600
(for: {3})
Manually check -25621/3600
(for: {3})
Manually check -11081/1600
(for: {3})
Manually check -8149/3600
(for: {3})
Manually check -7841/1600
(for: {3})
Manually check -6469/3600
(for: {3})
Manually check -2689/576
(for: {3})
Manually check -1849/576
(for: {3})
Manually check -421/144
(for: {3})
Manually check -301/144
(for: {3})
Manually check -29/16
(for: {3})
";
        // println!("\n\n{}\n\n{}",output, correct_output);
        assert!(output == correct_output);
    }

    #[bench]
    fn bench_possible_periods(b: &mut Bencher) {
        b.iter(|| {
            let fc = PolynomialInQ::from(
                    vec![Rational::one(), Rational::zero(), Rational::new(-29, 16)]
                );

            let res = possible_periods_search(fc, 2);

            assert!(res.unwrap().contains(&3));
        });
    }

    #[bench]
    fn bench_search_speed_test(b: &mut Bencher) {
        b.iter(|| {
            let mut output = "".to_string();
            (-3_000..=3_000i64).into_iter().for_each(|a| {
                for b in 1..=26i64 {
                    let b = b*b;
                    // Analytic bound (should be above or below gcd?)
                    let flo:f32 = (a as f32) / (b as f32);
                    if flo > -0.913942 {
                        continue;
                    }
                    let g = a.gcd(&b);
                    if g != 1 && g != -1 {
                        continue;
                    }
                    let fc = PolynomialInQ::from(
                        vec![Rational::one(), Rational::zero(), Rational::new(a, b)]
                    );

                    let res = possible_periods_search(fc, 2);

                    if let Some(possibilities) = res {
                        output += format!("Manually check {}/{}\n(for: {:?})\n", a, b, possibilities).as_ref();
                    }
                }
            });

            let correct_output = 
r"Manually check -2689/576
(for: {3})
Manually check -1849/576
(for: {3})
Manually check -421/144
(for: {3})
Manually check -301/144
(for: {3})
Manually check -29/16
(for: {3})
";
            println!("\n\n{}\n\n{}",output, correct_output);
            assert!(correct_output == output);
        })
    }
}