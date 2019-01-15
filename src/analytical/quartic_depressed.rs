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

/// Solves a depressed quartic equation x^4 + a2*x^2 + a1*x + a0 = 0.
///
/// Returned roots are ordered. Precision is about 1e-14 for f64.
///
/// # Examples
///
/// ```
/// use roots::find_roots_quartic_depressed;
///
/// let one_root = find_roots_quartic_depressed(1f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^4 = 0' has one root 0
///
/// let two_roots = find_roots_quartic_depressed(1f32, 0f32, -1f32);
/// // Returns Roots::Two([-1f32, 1f32]) as 'x^4 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_quartic_depressed<F: FloatType>(a2: F, a1: F, a0: F) -> Roots<F> {
    // Handle non-standard cases
    if a1 == F::zero() {
        // a1 = 0; x^4 + a2*x^2 + a0 = 0; solve biquadratic equation
        super::biquadratic::find_roots_biquadratic(F::one(), a2, a0)
    } else if a0 == F::zero() {
        // a0 = 0; x^4 + a2*x^2 + a1*x = 0; reduce to normalized cubic and add zero root
        super::cubic_normalized::find_roots_cubic_normalized(F::zero(), a2, a1).add_new_root(F::zero())
    } else {
        // Solve the auxiliary equation y^3 + (5/2)*a2*y^2 + (2*a2^2-a0)*y + (a2^3/2 - a2*a0/2 - a1^2/8) = 0
        let a2_pow_2 = a2 * a2;
        let a1_div_2 = a1 / F::two();
        let (b2, b1, b0) = (
            a2 * F::five() / F::two(),
            F::two() * a2_pow_2 - a0,
            (a2_pow_2 * a2 - a2 * a0 - a1_div_2 * a1_div_2) / F::two(),
        );

        // At least one root always exists. The last root is the maximal one.
        let y = *super::cubic_normalized::find_roots_cubic_normalized(b2, b1, b0)
            .as_ref()
            .iter()
            .last()
            .unwrap();

        let _a2_plus_2y = a2 + F::two() * y;
        if _a2_plus_2y > F::zero() {
            let sqrt_a2_plus_2y = _a2_plus_2y.sqrt();
            let q0a = a2 + y - a1_div_2 / sqrt_a2_plus_2y;
            let q0b = a2 + y + a1_div_2 / sqrt_a2_plus_2y;

            let mut roots = super::quadratic::find_roots_quadratic(F::one(), sqrt_a2_plus_2y, q0a);
            for x in super::quadratic::find_roots_quadratic(F::one(), -sqrt_a2_plus_2y, q0b)
                .as_ref()
                .iter()
            {
                roots = roots.add_new_root(*x);
            }
            roots
        } else {
            Roots::No([])
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_quartic_depressed() {
        assert_eq!(find_roots_quartic_depressed(0f32, 0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_quartic_depressed(1f32, 1f32, 1f32), Roots::No([]));

        // Thanks WolframAlpha for the test data
        match find_roots_quartic_depressed(1f64, 1f64, -1f64) {
            Roots::Two(x) => {
                assert_float_array_eq!(1e-15, x, [-1f64, 0.5698402909980532659114f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_quartic_depressed(-10f64, 5f64, 1f64) {
            Roots::Four(x) => {
                assert_float_array_eq!(
                    1e-15,
                    x,
                    [
                        -3.3754294311910698f64,
                        -0.1531811728532153f64,
                        0.67861075799846644f64,
                        2.84999984604581877f64
                    ]
                );
            }
            _ => {
                assert!(false);
            }
        }
    }
}
