#![feature(euclidean_division, test)]
#![allow(non_snake_case)]
mod math;
mod util;
mod ds_helper;

use math::*;
use util::*;
use ds_helper::*;
use rayon::prelude::*;

use std::fmt;

use num_integer::{Integer};
use std::f32::consts::PI;

extern crate num_complex;
extern crate num_integer;
extern crate num_rational;
extern crate rayon;
extern crate test;
extern crate arrayvec;
#[macro_use] extern crate lazy_static;
extern crate clap;
extern crate fnv;
use clap::{Arg, App, SubCommand};

pub fn format_search_result<T: fmt::Display, V: fmt::Debug>(ch: T, set: V) -> String {
    format!("Check {}, since I can't rule out periods in: {:?}", ch, set)
}

// Search through a given parameter space
// (uses z^4 + c, with appropriate optimizations)
pub fn search_z4_opt(height_max: i64, height_min: i64) {
    let bmin = ((height_max as f32).sqrt().sqrt() / 2.0).floor() as i64;
    (1..=bmin).into_par_iter().for_each(|b| {
        println!("{}", b);
        for a in -height_max..=height_max {
            let b = 2*b;
            let b = b*b*b*b;
            if a <= height_min && b <= height_min {
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
            
            let c = Rational::new(a, b);

            let res = z4c_possible_periods_search(c, 2);

            if let Some(possibilities) = res {
                println!("{}", format_search_result(c, possibilities));
            }
        }
    });
    println!("Completed search!");
}

pub fn search_z3_opt(height_max: i64, _height_min: i64) {
    let bmax = ((height_max as f32).cbrt()).floor() as i64;
    (-height_max..=height_max).into_par_iter().for_each(|num_a| {
            // println!("{}", num_a);
            for num_b in -height_max..=height_max {
                for denom_a in -bmax..=bmax {
                    for denom_b in -bmax..=bmax {
                        let denom = EisensteinInteger::new(denom_a, denom_b);
                        let denom = denom.product(denom.product(denom));
                        let numer = EisensteinInteger::new(num_a, num_b);
                        let c = QwElement::new(numer, denom);
                        if c.approx_coords().to_polar().1 > (PI / 6.0) {
                            continue;
                        }
                        if denom.is_zero() || numer.is_zero() {
                            continue;
                        }
                        if !numer.gcd(&denom).is_unit() {
                            continue;
                        }

                        if let Some(set) = z3c_possible_periods_search(c, 1) {
                            println!("{}", format_search_result(c, set));
                        }
                    }
                }
            }
        });
}

// Setup the command line interface
fn main() {
    let hmax_arg = Arg::with_name("height_max")
                        .short("hmax")
                        .long("height_max")
                        .help("Maximum height of c values to check")
                        .takes_value(true)
                        .required(true);
    let hmin_arg = Arg::with_name("height_min")
                        .short("hmin")
                        .long("height_min")
                        .help("Minimum height of c values to check")
                        .takes_value(true)
                        .default_value("0");
    let matches = App::new("Large Period Searcher")
            .version("0.1")
            .author("Matt Torrence <torrma01@gettysburg.edu>")
            .about("Arithmetic Dynamics tool for finding periodic points")
            .subcommand(SubCommand::with_name("z4c")
                .about("Search z^4 + c with standard optimizations / reductions")
                .arg(hmax_arg.clone())
                .arg(hmin_arg.clone()))
            .subcommand(SubCommand::with_name("z3c")
                .about("Search z^3 + c over Q(w) with standard optimizations / reductions")
                .arg(hmax_arg)
                .arg(hmin_arg))
            .get_matches();
    
    if let Some(matches) = matches.subcommand_matches("z4c") {
        let hmax: i64 = matches
                            .value_of("height_max")
                            .unwrap()
                            .parse()
                            .unwrap_or_else(|a| panic!("Error parsing, expected integer: {}", a));
        let hmin: i64 = matches
                            .value_of("height_min")
                            .unwrap()
                            .parse()
                            .unwrap_or_else(|a| panic!("Error parsing, expected integer: {}", a));
        search_z4_opt(hmax, hmin);
    }
    if let Some(matches) = matches.subcommand_matches("z3c") {
        let hmax: i64 = matches
                            .value_of("height_max")
                            .unwrap()
                            .parse()
                            .unwrap_or_else(|a| panic!("Error parsing, expected integer: {}", a));
        let hmin: i64 = matches
                            .value_of("height_min")
                            .unwrap()
                            .parse()
                            .unwrap_or_else(|a| panic!("Error parsing, expected integer: {}", a));
        search_z3_opt(hmax, hmin);
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;
    use test::black_box;

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
    fn bench_faster_conditions(b: &mut Bencher) {
        b.iter(|| {
            let a = -5649488755i64;
            let b = 639128961i64;
            let flo:f32 = (a as f32) / (b as f32);
            if flo > -0.913942 {
                panic!();
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                panic!();
            }
            black_box(g);

            let a = -5649488753i64;
            if flo > -0.913942 {
                panic!();
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                panic!();
            }
            g
        })
    }

    #[bench]
    fn bench_faster_search(b: &mut Bencher) {
        b.iter(|| {

            let res = z4c_possible_periods_search(Rational::new(-5649488755,639128961), 1);
            
            assert!(!res.is_none());
            assert!(res.unwrap().contains(&2));

            let res = z4c_possible_periods_search(Rational::new(-5649488753,639128961), 1);

            assert!(res.is_none());
        })
    }

    #[bench]
    fn bench_conditions(ben: &mut Bencher) {
        let (a, b) = (-3749999571i128, 3906250000i128);
        ben.iter(|| {
            let flo:f32 = (a as f32) / (b as f32);
            if flo > -0.913942 {
                panic!();
            }
            let g = a.gcd(&b);
            if g != 1 && g != -1 {
                panic!();
            }
            // Check for 2 mod 3 or 4 mod 5 cases
            // Does a/b have reduction mod 3?
            if b % 3 != 0 {
                // is a/b = 2 mod 3?
                if (a * mod_inverse(b % 3, 3)).rem_euclid(3) == 2 {
                    panic!();
                }
            }
            // Does a/b have reduction mod 5?
            if b % 5 != 0 {
                // is a/b = 4 mod 5?
                if (a * mod_inverse(b % 5, 5)).rem_euclid(5) == 4 {
                    panic!();
                }
            }
        });
    }

    #[bench]
    fn bench_after_conditions(ben: &mut Bencher) {
        let (a, b) = (-3749999571i64, 3906250000i64);
        ben.iter(|| {
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(a, b)]
            );

            let res = possible_periods_search(fc, 2);

            assert!(res.is_none());
        });
    }

    #[bench]
    fn bench_reduction(ben: &mut Bencher) {
        let (a, b) = (-3749999571i64, 3906250000i64);
        ben.iter(|| {
            let fc = PolynomialInQ::from(
                vec![Rational::one(), Rational::zero(), Rational::zero(), Rational::zero(), Rational::new(a, b)]
            );

            assert!(fc.has_good_reduction(3));

            let fc_red = fc.do_reduction(3);

            black_box(fc_red);
        })
    }

    #[bench]
    fn bench_search_speed_test(b: &mut Bencher) {
        b.iter(|| {
            let mut output = "".to_string();
            (-3_000..=3_000i64).into_iter().for_each(|a| {
                for b in 1..=26i64 {
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