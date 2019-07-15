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

// TODO: Maybe use algorithm 4.79 from http://cacr.uwaterloo.ca/hac/about/chap4.pdf ?
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

pub fn small_prime(n: usize) -> bool {
    for x in 2..=(n as f32).sqrt().ceil() as usize {
        if n % x == 0 {
            return false;
        }
    }
    return true;
}