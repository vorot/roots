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

/// Solves a normalized cubic equation x^3 + a2*x^2 + a1*x + a0 = 0.
///
/// In case more than one roots are present, they are returned in the increasing order.
///
/// # Examples
///
/// ```
/// use roots::find_roots_cubic_normalized;
///
/// let one_root = find_roots_cubic_normalized(0f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^3 = 0' has one root 0
///
/// let three_roots = find_roots_cubic_normalized(0f32, -1f32, 0f32);
/// // Returns Roots::Three([-1f32, -0f32, 1f32]) as 'x^3 - x = 0' has roots -1, 0, and 1
/// ```
pub fn find_roots_cubic_normalized<F: FloatType>(a2: F, a1: F, a0: F) -> Roots<F> {
    let q = (F::three() * a1 - a2 * a2) / F::nine();
    let r = (F::nine() * a2 * a1 - F::twenty_seven() * a0 - F::two() * a2 * a2 * a2) / (F::two() * F::twenty_seven());
    let q3 = q * q * q;
    let d = q3 + r * r;
    let a2_div_3 = a2 / F::three();

    if d < F::zero() {
        let phi_3 = (r / (-q3).sqrt()).acos() / F::three();
        let sqrt_q_2 = F::two() * (-q).sqrt();

        Roots::One([sqrt_q_2 * phi_3.cos() - a2_div_3])
            .add_new_root(sqrt_q_2 * (phi_3 - F::two_third_pi()).cos() - a2_div_3)
            .add_new_root(sqrt_q_2 * (phi_3 + F::two_third_pi()).cos() - a2_div_3)
    } else {
        let sqrt_d = d.sqrt();
        let s = (r + sqrt_d).cbrt();
        let t = (r - sqrt_d).cbrt();

        if s == t {
            if s + t == F::zero() {
                Roots::One([s + t - a2_div_3])
            } else {
                Roots::One([s + t - a2_div_3]).add_new_root(-(s + t) / F::two() - a2_div_3)
            }
        } else {
            Roots::One([s + t - a2_div_3])
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_cubic_normalized() {
        assert_eq!(find_roots_cubic_normalized(0f32, 0f32, 0f32), Roots::One([0f32]));

        match find_roots_cubic_normalized(0f64, -1f64, 0f64) {
            Roots::Three(x) => {
                assert_float_array_eq!(1e-15, x, [-1f64, 0f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_cubic_normalized(1f64, -2f64, 2f64) {
            Roots::One(x) => {
                assert_float_array_eq!(1e-15, x, [-2.269530842081142770853135f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_cubic_normalized(0f64, -3f64, 2f64) {
            Roots::Two(x) => {
                assert_float_array_eq!(1e-15, x, [-2f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_cubic_normalized(-2f64, -3f64, 2f64) {
            Roots::Three(x) => {
                assert_float_array_eq!(
                    1e-15,
                    x,
                    [
                        -1.342923082777170208054859f64,
                        0.5293165801288393926136314f64,
                        2.813606502648330815441228f64
                    ]
                );
            }
            _ => {
                assert!(false);
            }
        }
    }
}
