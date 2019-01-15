// Copyright (c) 2015, Mikhail Vorotilov
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use super::super::FloatType;
use super::super::Roots;

/// Solves a depressed cubic equation x^3 + a1*x + a0 = 0.
///
/// In case more than one roots are present, they are returned in the increasing order.
///
/// # Examples
///
/// ```
/// use roots::find_roots_cubic_depressed;
///
/// let one_root = find_roots_cubic_depressed(0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^3 = 0' has one root 0
///
/// let three_roots = find_roots_cubic_depressed(-1f32, 0f32);
/// // Returns Roots::Three([-1f32, -0f32, 1f32]) as 'x^3 - x = 0' has roots -1, 0, and 1
/// ```
pub fn find_roots_cubic_depressed<F: FloatType>(a1: F, a0: F) -> Roots<F> {
    if a1 == F::zero() {
        Roots::One([-a0.cbrt()])
    } else if a0 == F::zero() {
        super::quadratic::find_roots_quadratic(F::one(), F::zero(), a1).add_new_root(F::zero())
    } else {
        let d = a0 * a0 / F::four() + a1 * a1 * a1 / F::twenty_seven();
        if d < F::zero() {
            // n*a0^2 + m*a1^3 < 0 => a1 < 0
            let a = (-F::four() * a1 / F::three()).sqrt();

            let phi = (-F::four() * a0 / (a * a * a)).acos() / F::three();
            Roots::One([a * phi.cos()])
                .add_new_root(a * (phi + F::two_third_pi()).cos())
                .add_new_root(a * (phi - F::two_third_pi()).cos())
        } else {
            let sqrt_d = d.sqrt();
            let a0_div_2 = a0 / F::two();
            let x1 = (sqrt_d - a0_div_2).cbrt() - (sqrt_d + a0_div_2).cbrt();
            if d == F::zero() {
                // one real root and one double root
                Roots::One([x1]).add_new_root(a0_div_2)
            } else {
                // one real root
                Roots::One([x1])
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_cubic_depressed() {
        assert_eq!(find_roots_cubic_depressed(0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_cubic_depressed(-1f64, 0f64), Roots::Three([-1f64, 0f64, 1f64]));

        match find_roots_cubic_depressed(-2f64, 2f64) {
            Roots::One(x) => {
                assert_float_array_eq!(1e-15, x, [-1.769292354238631415240409f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_cubic_depressed(-3f64, 2f64) {
            Roots::Two(x) => {
                assert_float_array_eq!(1e-15, x, [-2f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_cubic_depressed(-2f64, 1f64) {
            Roots::Three(x) => {
                assert_float_array_eq!(1e-15, x, [(-1f64 - 5f64.sqrt()) / 2f64, (-1f64 + 5f64.sqrt()) / 2f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }
    }
}
