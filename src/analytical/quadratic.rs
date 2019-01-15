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

/// Solves a quadratic equation a2*x^2 + a1*x + a0 = 0.
///
/// In case two roots are present, the first returned root is less than the second one.
///
/// # Examples
///
/// ```
/// use roots::Roots;
/// use roots::find_roots_quadratic;
///
/// let no_roots = find_roots_quadratic(1f32, 0f32, 1f32);
/// // Returns Roots::No([]) as 'x^2 + 1 = 0' has no roots
///
/// let one_root = find_roots_quadratic(1f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^2 = 0' has one root 0
///
/// let two_roots = find_roots_quadratic(1f32, 0f32, -1f32);
/// // Returns Roots::Two([-1f32,1f32]) as 'x^2 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_quadratic<F: FloatType>(a2: F, a1: F, a0: F) -> Roots<F> {
    // Handle non-standard cases
    if a2 == F::zero() {
        // a2 = 0; a1*x+a0=0; solve linear equation
        super::linear::find_roots_linear(a1, a0)
    } else {
        // Rust lacks a simple way to convert an integer constant to generic type F
        let discriminant = a1 * a1 - F::four() * a2 * a0;
        if discriminant < F::zero() {
            Roots::No([])
        } else {
            let a2x2 = F::two() * a2;
            if discriminant == F::zero() {
                Roots::One([-a1 / a2x2])
            } else {
                // To improve precision, do not use the smallest divisor.
                // See https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
                let sq = discriminant.sqrt();

                let (same_sign, diff_sign) = if a1 < F::zero() {
                    (-a1 + sq, -a1 - sq)
                } else {
                    (-a1 - sq, -a1 + sq)
                };

                let (x1, x2) = if same_sign.abs() > a2x2.abs() {
                    let a0x2 = F::two() * a0;
                    if diff_sign.abs() > a2x2.abs() {
                        // 2*a2 is the smallest divisor, do not use it
                        (a0x2 / same_sign, a0x2 / diff_sign)
                    } else {
                        // diff_sign is the smallest divisor, do not use it
                        (a0x2 / same_sign, same_sign / a2x2)
                    }
                } else {
                    // 2*a2 is the greatest divisor, use it
                    (diff_sign / a2x2, same_sign / a2x2)
                };

                // Order roots
                if x1 < x2 {
                    Roots::Two([x1, x2])
                } else {
                    Roots::Two([x2, x1])
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_quadratic() {
        assert_eq!(find_roots_quadratic(0f32, 0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_quadratic(1f32, 0f32, 1f32), Roots::No([]));
        assert_eq!(find_roots_quadratic(1f64, 0f64, -1f64), Roots::Two([-1f64, 1f64]));
    }

    #[test]
    fn test_find_roots_quadratic_small_a2() {
        assert_eq!(
            find_roots_quadratic(1e-20f32, -1f32, -1e-30f32),
            Roots::Two([-1e-30f32, 1e20f32])
        );
        assert_eq!(
            find_roots_quadratic(-1e-20f32, 1f32, 1e-30f32),
            Roots::Two([-1e-30f32, 1e20f32])
        );
        assert_eq!(find_roots_quadratic(1e-20f32, -1f32, 1f32), Roots::Two([1f32, 1e20f32]));
        assert_eq!(find_roots_quadratic(-1e-20f32, 1f32, 1f32), Roots::Two([-1f32, 1e20f32]));
        assert_eq!(find_roots_quadratic(-1e-20f32, 1f32, -1f32), Roots::Two([1f32, 1e20f32]));
    }

    #[test]
    fn test_find_roots_quadratic_big_a1() {
        assert_eq!(find_roots_quadratic(1f32, -1e15f32, -1f32), Roots::Two([-1e-15f32, 1e15f32]));
        assert_eq!(find_roots_quadratic(-1f32, 1e15f32, 1f32), Roots::Two([-1e-15f32, 1e15f32]));
    }
}
