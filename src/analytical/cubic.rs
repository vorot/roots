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

use super::super::FloatWithConstants;
use super::super::Roots;

/// Solves a cubic equation a3*x^3 + a2*x^2 + a1*x + a0 = 0.
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
/// ```
pub fn find_roots_cubic<F:FloatWithConstants>(a3:F, a2:F, a1:F, a0:F) -> Roots<F> {
  // Handle non-standard cases
  if a3 == F::zero() {
    // a3 = 0; a2*x^2+a1*x+a0=0; solve quadratic equation
    super::quadratic::find_roots_quadratic(a2, a1, a0)
  } else if a2 == F::zero() {
    // a2 = 0; a3*x^3+a1*x+a0=0; solve depressed cubic equation
    super::cubic_depressed::find_roots_cubic_depressed(a1/a3, a0/a3)
  } else {
    // solve normalized cubic expression
    super::cubic_normalized::find_roots_cubic_normalized(a2/a3, a1/a3, a0/a3)
  }
}

#[test]
fn test_find_roots_cubic() {
  assert_eq!(find_roots_cubic(1f32, 0f32, 0f32, 0f32), Roots::One([0f32]));

  match find_roots_cubic(1f64, 0f64, -1f64, 0f64) {
    Roots::Three([x1, x2, x3]) => {
      assert_float_eq!(1e-15, x1, -1f64 );
      assert_float_eq!(1e-15, x2, 0f64 );
      assert_float_eq!(1e-15, x3, 1f64 );
    },
    _ => { assert!(false); }
  }
}
