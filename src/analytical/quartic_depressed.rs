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

use std::num::Float;
use std::cmp::Ordering;

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
/// // Returns (Some(0f64), None, None, None) as 'x^4 = 0' has one root 0
///
/// let two_roots = find_roots_quartic_depressed(1f32, 0f32, -1f32);
/// // Returns (Some(-1f32), Some(1f32), None, None) as 'x^4 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_quartic_depressed<F:Float>(a2:F, a1:F, a0:F) -> Vec<F> {
  // Handle non-standard cases
  let mut roots = if a1 == F::zero() {
    // a1 = 0; x^4 + a2*x^2 + a0 = 0; solve biquadratic equation
    super::biquadratic::find_roots_biquadratic(F::one(), a2, a0)
  } else if a0 == F::zero() {
    // a0 = 0; x^4 + a2*x^2 + a1*x = 0; reduce to normalized cubic and arrange results
    let mut tmp = vec![F::zero()];
    tmp.push_all(super::cubic_normalized::find_roots_cubic_normalized(F::zero(), a2, a1).as_slice());
    tmp
  } else {
    // Solve the auxiliary equation y^3 + (5/2)*a2*y^2 + (2*a2^2-a0)*y + (a2^3/2 - a2*a0/2 - a1^2/8) = 0
    let _2 = F::one() + F::one();
    let _3 = _2 + F::one();
    let _4 = _2 + _2;
    let _5 = _2 + _3;
    let _a2_2 = a2*a2;
    let _a2_div_2 = a2 / _2;

    // last root is the maximal one
    let y = super::cubic_normalized::find_roots_cubic_normalized(a2*_5/_2, _2*_a2_2-a0, (_a2_2*a2 - a2*a0 - a1*a1/_4)/_2).pop().unwrap();

    let _a2_plus_2y = a2 + _2*y;
    if _a2_plus_2y > F::zero() {
      let sqrt_a2_plus_2y = _a2_plus_2y.sqrt();
      let q1 = sqrt_a2_plus_2y;
      let q0a = a2 + y - a1/(_2*sqrt_a2_plus_2y);
      let q0b = a2 + y + a1/(_2*sqrt_a2_plus_2y);

      let mut tmp = super::quadratic::find_roots_quadratic(F::one(), q1, q0a);
      tmp.push_all(super::quadratic::find_roots_quadratic(F::one(), -q1, q0b).as_slice());
      tmp
    } else {
      vec![]
    }
  };

  roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
  roots.dedup();
  roots
}

#[test]
fn test_find_roots_biquadratic() {
  assert_eq!(find_roots_quartic_depressed(0f32, 0f32, 0f32), [0f32]);
  assert_eq!(find_roots_quartic_depressed(1f32, 1f32, 1f32), []);

  // Thanks WolframAlpha for the test data
  match find_roots_quartic_depressed(1f64, 1f64, -1f64).as_slice() {
    [x1, x2] => {
      assert_float_eq!(1e-15, x1, -1f64 );
      assert_float_eq!(1e-15, x2, 0.5698402909980532659114f64 );
    },
    _ => { assert!(false); }
  }

  match find_roots_quartic_depressed(-10f64, 5f64, 1f64).as_slice() {
    [x1, x2, x3, x4] => {
      assert_float_eq!(1e-15, x1, -3.3754294311910698f64 );
      assert_float_eq!(1e-15, x2, -0.1531811728532153f64 );
      assert_float_eq!(1e-15, x3, 0.67861075799846644f64 );
      assert_float_eq!(1e-15, x4, 2.84999984604581877f64 );
    },
    _ => { assert!(false); }
  }
}
