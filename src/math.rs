use num_integer::{Integer};

pub fn mod_inverse<T: Integer + Copy>(num: T, prime: T) -> T {
    let mut a : T = prime;
    let mut b : T = num;
    let mut x : T = T::one();
    let mut y : T = T::zero();
    while b != T::zero() {
        let t = b;
        let q = a / t;
        b = a - q*t;
        a = t;
        let t = x;
        x = y - q*t;
        y = t;
    }
    
    if y < T::zero() {
        y + prime
    } else {
        y
    }
}

#[allow(unused)]
pub fn multiplicative_order<T: Integer + Copy>(a: T, p: T) -> u32 {
    let mut res = 1;
    let mut curr = a;
    debug_assert!(a > T::zero());
    while curr != T::one() {
        curr = curr * a;
        curr = curr % p;
        res += 1;
    }
    res
}

lazy_static! {
    static ref PRIMES: [bool; 301] = {
        let mut p = [false; 301];
        for i in 0..301 {
            p[i] = small_prime(i);
        }
        p
    };
}

#[inline]
pub fn prime(n: usize) -> bool {
    PRIMES[n]
}

pub fn small_prime(n: usize) -> bool {
    for x in 2..=(n as f32).sqrt().ceil() as usize {
        if n % x == 0 {
            return false;
        }
    }
    return true;
}

const fn num_bits<T>() -> usize { std::mem::size_of::<T>() * 8 }

pub fn log_2(x: i64) -> u32 {
    assert!(x > 0);
    num_bits::<i64>() as u32 - x.leading_zeros() - 1
}

// Return a^b mod p
pub fn mod_power(a: i64, b: i64, p: i64) -> i64 {
    let mut res = 1i64;
    for _ in 0..b {
        let r = res.checked_mul(a);
        match r {
            Some(val) => res = val,
            None => {
                res %= p;
                res *= a;
            }
        }
    }
    res.rem_euclid(p)
}

fn cipolla_mult(ab: (i64, i64), cd: (i64, i64), w: i64, p: i64) -> (i64, i64) {
    let (a, b) = ab;
    let (c, d) = cd;
    ((a*c + b*d*w)%p, (a*d + b*c)%p)
}

// Take in an integer n and odd prime p
// return both square roots of n mod p as (a, b)
// or None if no roots exist
pub fn cipolla(n: i64, p: i64) -> Option<(i64, i64)> {
    debug_assert!(p > 1);
    let n = n.rem_euclid(p);
    if n == 0 || n == 1 {
        return Some((n, (-n).rem_euclid(p)));
    }
    let phi = p - 1;
    if mod_power(n, phi/2, p) != 1 {
        return None;
    }
    if p % 4 == 3 {
        let ans = mod_power(n, (p + 1)/4, p);
        return Some((ans, (-ans).rem_euclid(p)));
    }
    let mut aa = 0;
    for i in 1..p {
        let temp = mod_power((i*i - n)%p, phi/2, p);
        if temp == phi {
            aa = i;
            break;
        }
    }
    let exponent = (p + 1)/2;
    let mut x1 = (aa, 1);
    let mut x2 = cipolla_mult(x1, x1, aa*aa - n, p);
    let l = log_2(exponent);
    for i in (0..l).rev() {
        if (exponent & (1 << i)) == 0 {
            x2 = cipolla_mult(x2, x1, aa*aa-n, p);
            x1 = cipolla_mult(x1, x1, aa*aa-n, p); 
        } else {
            x1 = cipolla_mult(x1, x2, aa*aa-n, p);
            x2 = cipolla_mult(x2, x2, aa*aa-n, p);
        }
    }
    Some((x1.0, (-x1.0).rem_euclid(p)))
}

// Check if -3 is a quadratic residue mod p
pub fn has_qw_homomorphism(p: i64) -> bool {
    p == 3 || mod_power(p-3, (p - 1)/2, p) == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cipolla() {
        println!("(2, 7)");
        assert_eq!(cipolla(2, 7), Some((4, 3)));
        println!("(8218, 10007)");
        assert_eq!(cipolla(8218, 10007), Some((9872, 135)));
        println!("(56, 101)");
        assert_eq!(cipolla(56, 101), Some((37, 64)));
        println!("(1, 11)");
        assert_eq!(cipolla(1, 11), Some((1, 10)));
        println!("(8219, 10007)");
        assert_eq!(cipolla(8219, 10007), None);

        println!("{:?}", cipolla(0, 3));
    }
}