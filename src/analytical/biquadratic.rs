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
/// // Returns [] as 'x^4 + 1 = 0' has no roots
///
/// let one_root = find_roots_biquadratic(1f64, 0f64, 0f64);
/// // Returns [0f64] as 'x^4 = 0' has one root 0
///
/// let two_roots = find_roots_biquadratic(1f32, 0f32, -1f32);
/// // Returns [-1f32, 1f32] as 'x^4 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_biquadratic<F:Float>(a4:F, a2:F, a0:F) -> Vec<F> {
  // Handle non-standard cases
  if a4 == F::zero() {
    // a4 = 0; a2*x^2 + a0 = 0; solve quadratic equation
    super::quadratic::find_roots_quadratic(a2, F::zero(), a0)
  } else {
    // solve the corresponding quadratic equation and order roots
    // The function find_roots_quadratic returns an increasing tuple
    // Match evaluates branches consequently
    let mut roots = match super::quadratic::find_roots_quadratic(a4, a2, a0).as_slice() {
      [x1] if x1 >= F::zero() => vec![-x1.sqrt(), x1.sqrt()],
      [x1,x2] if x1 >= F::zero() => vec![-x2.sqrt(), -x1.sqrt(), x1.sqrt(), x2.sqrt()],
      [_,x2] if x2 >= F::zero() => vec![-x2.sqrt(), x2.sqrt()],
      _ => vec![],
    };

    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    roots.dedup();
    roots
  }
}

#[test]
fn test_find_roots_biquadratic() {
  assert_eq!(find_roots_biquadratic(0f32, 0f32, 0f32), [0f32]);
  assert_eq!(find_roots_biquadratic(1f32, 0f32, 1f32), []);
  assert_eq!(find_roots_biquadratic(1f64, 0f64, -1f64), [-1f64, 1f64]);
  assert_eq!(find_roots_biquadratic(1f64, -5f64, 4f64), [-2f64, -1f64, 1f64, 2f64]);
}
