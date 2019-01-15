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

/// Solves a bi-quadratic equation a4*x^4 + a2*x^2 + a0 = 0.
///
/// Returned roots are arranged in the increasing order.
///
/// # Examples
///
/// ```
/// use roots::find_roots_biquadratic;
///
/// let no_roots = find_roots_biquadratic(1f32, 0f32, 1f32);
/// // Returns Roots::No([]) as 'x^4 + 1 = 0' has no roots
///
/// let one_root = find_roots_biquadratic(1f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^4 = 0' has one root 0
///
/// let two_roots = find_roots_biquadratic(1f32, 0f32, -1f32);
/// // Returns Roots::Two([-1f32, 1f32]) as 'x^4 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_biquadratic<F: FloatType>(a4: F, a2: F, a0: F) -> Roots<F> {
    // Handle non-standard cases
    if a4 == F::zero() {
        // a4 = 0; a2*x^2 + a0 = 0; solve quadratic equation
        super::quadratic::find_roots_quadratic(a2, F::zero(), a0)
    } else if a0 == F::zero() {
        // a0 = 0; a4*x^4 + a2*x^2 = 0; solve quadratic equation and add zero root
        super::quadratic::find_roots_quadratic(a4, F::zero(), a2).add_new_root(F::zero())
    } else {
        // solve the corresponding quadratic equation and order roots
        let mut roots = Roots::No([]);
        for x in super::quadratic::find_roots_quadratic(a4, a2, a0).as_ref().iter() {
            if *x > F::zero() {
                let sqrt_x = x.sqrt();
                roots = roots.add_new_root(-sqrt_x).add_new_root(sqrt_x);
            } else if *x == F::zero() {
                roots = roots.add_new_root(F::zero());
            }
        }
        roots
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_biquadratic() {
        assert_eq!(find_roots_biquadratic(0f32, 0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_biquadratic(1f32, 0f32, 1f32), Roots::No([]));
        assert_eq!(find_roots_biquadratic(1f64, 0f64, -1f64), Roots::Two([-1f64, 1f64]));
        assert_eq!(
            find_roots_biquadratic(1f64, -5f64, 4f64),
            Roots::Four([-2f64, -1f64, 1f64, 2f64])
        );
    }
}
