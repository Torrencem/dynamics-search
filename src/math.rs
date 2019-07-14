use num_integer::{Integer};
use std;

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

pub fn point_from_hash<T: Integer + Copy>(value: T, prime: T, dimension: usize) -> Vec<T> {
    let mut val = value;

    let mut res = Vec::with_capacity(dimension);
    
    for _ in 0..=dimension {
        res.push(val % prime);
        val = val / prime;
    }

    res
}

pub fn point_to_hash<T: Integer + Copy>(point: Vec<T>, prime: T) -> T {
    let mut hash_q = T::zero();

    for coefficient in point.into_iter().rev() {
        hash_q = hash_q * prime + coefficient;
    }

    hash_q
}

pub fn enum_points(prime: u64, dimension: usize) -> impl Iterator<Item=Vec<u64>> {
    let mut current_range = 1u64;
    let highest_range = prime.pow(dimension as u32);
    let mut value = 1u64;

    std::iter::from_fn(move || {
        if value >= 2*current_range {
            current_range *= prime;
            value = current_range;
            if current_range > highest_range {
                return None;
            }
        }
        let oldval = value;
        value += 1;
        Some(point_from_hash(oldval, prime, dimension))
    })
}

pub fn normalize_coordinates<T: Integer + Copy>(point: &mut Vec<T>, prime: T) {
    let len_points = point.len();
    let mut last_coefficient = T::zero();

    for coefficient in 0..len_points {
        let val = point[coefficient] % prime;
        point[coefficient] = val;
        if val != T::zero() {
            last_coefficient = val;
        }
    }

    let m_inverse = mod_inverse(last_coefficient, prime);

    for coefficient in 0..len_points {
        point[coefficient] = (point[coefficient] * m_inverse) % prime;
    }
}

pub struct Powerset<T> where T: Clone {
    items: Vec<T>,
    state: u64,
    num_items: usize,
}

impl<T: Clone> Iterator for Powerset<T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Vec<T>> {
        if self.state == (2 << (self.num_items - 1)) {
            return None;
        }
        let curr_set: Vec<T> = self.items.iter()
                        .enumerate()
                        .filter(|&(t, _)| (self.state >> t) % 2 == 1)
                        .map(|(_, element)| element.clone())
                        .collect();
        
        self.state += 1;
        Some(curr_set)
    }
}

pub fn powerset<T>(s: Vec<T>) -> Powerset<T>
where T: Clone {
    let l = s.len();
    assert!(l <= 64);
    Powerset {items: s, state: 0, num_items: l}
}

pub fn lcm_all<T: Integer + Copy>(v: &[T]) -> T {
    let mut res = T::one();

    for x in v {
        res = res.lcm(x);
    }

    res
}

// TODO: Maybe use algorithm 4.79 from http://cacr.uwaterloo.ca/hac/about/chap4.pdf?
// This algorithm is linear time which is probably worse than it needs to be
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