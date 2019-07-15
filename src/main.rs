#![feature(euclidean_division, test)]
#![allow(non_snake_case)]
mod math;
mod util;
mod ds_helper;

use math::*;
use util::*;
use ds_helper::*;
use rayon::prelude::*;

use num_integer::{Integer};

extern crate num_integer;
extern crate rayon;
extern crate test;

pub fn run_large_search() {
    (-10_000_000_000..=10_000_000_000i128).into_par_iter().for_each(|a| {
        for b in 1..=158i128 {
            let b = 2*b;
            let b = b*b*b*b;
            // Did this come up in our old search?
            if a.abs() < 1_000_000_000 && b < 1_000_000_000 {
                continue;
            }
            let flo:f32 = (a as f32) / (b as f32);
            if flo > -0.913942 {
                continue;
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                continue;
            }
            // Check for 2 mod 3 or 4 mod 5 cases
            // Does a/b have reduction mod 3?
            if b % 3 != 0 {
                // is a/b = 2 mod 3?
                if (a * mod_inverse(b % 3, 3)).rem_euclid(3) == 2 {
                    continue;
                }
            }
            // Does a/b have reduction mod 5?
            if b % 5 != 0 {
                // is a/b = 4 mod 5?
                if (a * mod_inverse(b % 5, 5)).rem_euclid(5) == 4 {
                    continue;
                }
            }
            // Create z^4 + c from it's coefficients (1, 0, 0, 0, a/b)
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(a, b)]
            );

            let res = possible_periods_search(fc, 2);

            if let Some(possibilities) = res {
                println!("Manually check {}/{}\n(for: {:?})", a, b, possibilities);
            }
        }
    });
}

fn main() {
    run_large_search();

    println!("Completed search!");
}



#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn test_functionality() {
        let mut output = "".to_string();
        (-100_000..=100_000i128).into_iter().for_each(|a| {
            for b in 1..=60i128 {
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
    fn bench_large_example(b: &mut Bencher) {
        b.iter(|| {
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(-5649488755,639128961)]
            );

            let res = possible_periods_search(fc, 1);

            assert!(res.unwrap().contains(&2));

            // Wrong example

            let fc = PolynomialInQ::from(                                                                  //       v--- notice the 4
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(-5649488754,639128961)]
            );

            let res = possible_periods_search(fc, 1);

            assert!(res.is_none());
        })
    }

    #[bench]
    fn bench_search_speed_test(b: &mut Bencher) {
        b.iter(|| {
            let mut output = "".to_string();
            (-3_000..=3_000i128).into_iter().for_each(|a| {
                for b in 1..=26i128 {
                    let b = b*b;
                    // Random Analytic bound (doesn't actually apply for z^2 case?)
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