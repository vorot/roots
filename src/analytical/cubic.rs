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

/// Solves a cubic equation a3*x^3 + a2*x^2 + a1*x + a0 = 0.
///
/// General formula (complex numbers) is implemented for three roots.
///
/// Note that very small values of a3 (comparing to other coefficients) will cause the loss of precision.
///
/// In case more than one roots are present, they are returned in the increasing order.
///
/// # Examples
///
/// ```
/// use roots::Roots;
/// use roots::find_roots_cubic;
///
/// let no_roots = find_roots_cubic(0f32, 1f32, 0f32, 1f32);
/// // Returns Roots::No([]) as 'x^2 + 1 = 0' has no roots
///
/// let one_root = find_roots_cubic(1f64, 0f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^3 = 0' has one root 0
///
/// let three_roots = find_roots_cubic(1f32, 0f32, -1f32, 0f32);
/// // Returns Roots::Three([-1f32, 0f32, 1f32]) as 'x^3 - x = 0' has roots -1, 0, and 1
///
/// let three_roots_less_precision = find_roots_cubic(
///            -0.000000000000000040410628481035f64,
///            0.0126298310280606f64,
///            -0.100896606408756f64,
///            0.0689539597036461f64);
/// // Returns Roots::Three([0.7583841816097057f64, 7.233267996296344f64, 312537357195212.9f64])
/// // while online math expects 0.7547108770537f64, 7.23404258961f64, 312537357195213f64
/// ```
pub fn find_roots_cubic<F: FloatType>(a3: F, a2: F, a1: F, a0: F) -> Roots<F> {
    // Handle non-standard cases
    if a3 == F::zero() {
        // a3 = 0; a2*x^2+a1*x+a0=0; solve quadratic equation
        super::quadratic::find_roots_quadratic(a2, a1, a0)
    } else if a2 == F::zero() {
        // a2 = 0; a3*x^3+a1*x+a0=0; solve depressed cubic equation
        super::cubic_depressed::find_roots_cubic_depressed(a1 / a3, a0 / a3)
    } else if a3 == F::one() {
        // solve normalized cubic expression
        super::cubic_normalized::find_roots_cubic_normalized(a2, a1, a0)
    } else {
        let _2 = F::from(2i16);
        let _3 = F::from(3i16);
        let _4 = F::from(4i16);
        let _9 = F::from(9i16);
        let _18 = F::from(18i16);
        let _27 = F::from(27i16);

        // standard case
        let d = _18 * a3 * a2 * a1 * a0 - _4 * a2 * a2 * a2 * a0 + a2 * a2 * a1 * a1
            - _4 * a3 * a1 * a1 * a1
            - _27 * a3 * a3 * a0 * a0;
        let d0 = a2 * a2 - _3 * a3 * a1;
        let d1 = _2 * a2 * a2 * a2 - _9 * a3 * a2 * a1 + _27 * a3 * a3 * a0;
        if d < F::zero() {
            // one real root
            let sqrt = (-_27 * a3 * a3 * d).sqrt();
            let c = F::cbrt(if d1 < F::zero() { d1 - sqrt } else { d1 + sqrt } / _2);
            let x = -(a2 + c + d0 / c) / (_3 * a3);
            Roots::One([x])
        } else if d == F::zero() {
            // multiple roots
            if d0 == F::zero() {
                // triple root
                Roots::One([-a2 / (a3 * _3)])
            } else {
                // single root and double root
                Roots::One([(_9 * a3 * a0 - a2 * a1) / (d0 * _2)])
                    .add_new_root((_4 * a3 * a2 * a1 - _9 * a3 * a3 * a0 - a2 * a2 * a2) / (a3 * d0))
            }
        } else {
            // three real roots
            let c3_img = F::sqrt(_27 * a3 * a3 * d) / _2;
            let c3_real = d1 / _2;
            let c3_module = F::sqrt(c3_img * c3_img + c3_real * c3_real);
            let c3_phase = _2 * F::atan(c3_img / (c3_real + c3_module));
            let c_module = F::cbrt(c3_module);
            let c_phase = c3_phase / _3;
            let c_real = c_module * F::cos(c_phase);
            let c_img = c_module * F::sin(c_phase);
            let x0_real = -(a2 + c_real + (d0 * c_real) / (c_module * c_module)) / (_3 * a3);

            let e_real = -F::one() / _2;
            let e_img = F::sqrt(_3) / _2;
            let c1_real = c_real * e_real - c_img * e_img;
            let c1_img = c_real * e_img + c_img * e_real;
            let x1_real = -(a2 + c1_real + (d0 * c1_real) / (c1_real * c1_real + c1_img * c1_img)) / (_3 * a3);

            let c2_real = c1_real * e_real - c1_img * e_img;
            let c2_img = c1_real * e_img + c1_img * e_real;
            let x2_real = -(a2 + c2_real + (d0 * c2_real) / (c2_real * c2_real + c2_img * c2_img)) / (_3 * a3);

            Roots::One([x0_real]).add_new_root(x1_real).add_new_root(x2_real)
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_cubic() {
        assert_eq!(find_roots_cubic(1f32, 0f32, 0f32, 0f32), Roots::One([0f32]));

        match find_roots_cubic(1f64, 0f64, -1f64, 0f64) {
            Roots::Three(x) => {
                assert_float_array_eq!(1e-15, x, [-1f64, 0f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }
    }

    #[test]
    fn test_find_roots_cubic_small_discriminant() {
        // Try to find roots of the cubic polynomial where the highest coefficient is very small
        // (as reported by Andrew Hunter in July 2019)
        match find_roots_cubic(
            -0.000000000000000040410628481035f64,
            0.0126298310280606f64,
            -0.100896606408756f64,
            0.0689539597036461f64,
        ) {
            Roots::Three(x) => {
                // (According to Wolfram Alpha, roots must be 0.7547108770537f64, 7.23404258961f64, 312537357195213f64)
                // Actual result differ a little due to the limited precision of calculations.
                assert_float_array_eq!(1e-8, x, [0.7583841816097057f64, 7.233267996296344f64, 312537357195212.9f64]);
            }
            _ => {
                assert!(false);
            }
        }
    }
}
