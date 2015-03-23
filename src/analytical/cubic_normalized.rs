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
/// // Returns [0f64] as 'x^3 = 0' has one root 0
///
/// let three_roots = find_roots_cubic_normalized(0f32, -1f32, 0f32);
/// // Returns [-1f32, -0f32, 1f32] as 'x^3 - x = 0' has roots -1, 0, and 1
/// ```
pub fn find_roots_cubic_normalized<F:Float>(a2:F, a1:F, a0:F) -> Vec<F> {
  let _2 = F::one() + F::one();
  let _3 = _2 + F::one();
  let _9 = _3 * _3;
  let _27 = _9 * _3;
  let _54 = _27 * _2;

  let q = (_3*a1 - a2*a2) / _9;
  let r = (_9*a2*a1 - _27*a0 - _2*a2*a2*a2) / _54;
  let q3 = q*q*q;
  let d = q3 + r*r;
  let a2_3 = a2 / _3;

  let mut roots = if d < F::zero() {
    let phi_3 = (r/(-q3).sqrt()).acos() / _3;
    let sqrt_q_2 = _2*(-q).sqrt();
    let pi = F::acos(-F::one());
    let pi2_3 = pi*_2 / _3;
    let pi4_3 = pi2_3 *_2;

    vec![sqrt_q_2*phi_3.cos() - a2_3, sqrt_q_2*(phi_3 + pi2_3).cos() - a2_3, sqrt_q_2*(phi_3 + pi4_3).cos() - a2_3]
  } else {
    let sqrt_d = d.sqrt();
    let s = cbrt(r + sqrt_d);
    let t = cbrt(r - sqrt_d);

    if s == t {
      if s + t == F::zero() {
        vec![s + t - a2_3]
      } else {
        vec![s + t - a2_3, -(s + t) / _2 - a2_3]
      }
    } else {
      vec![s + t - a2_3]
    }
  };

  roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
  roots.dedup();
  roots
}

fn cbrt<F:Float>(x:F) -> F {
  // Rust lacks a simple way to convert an integer constant to generic type F
  let _1_3 = F::one() / (F::one() + F::one() + F::one());

  if x < F::zero() {
    -(-x).powf(_1_3)
  } else {
    x.powf(_1_3)
  }
}

#[test]
fn test_cbrt() {
  assert_eq!(cbrt(-8f64), -2f64);
}

#[test]
fn test_find_roots_cubic_normalized() {
  assert_eq!(find_roots_cubic_normalized(0f32, 0f32, 0f32), [0f32]);

  match find_roots_cubic_normalized(0f64, -1f64, 0f64).as_slice() {
    [x1, x2, x3] => {
      assert_float_eq!(1e-15, x1, -1f64 );
      assert_float_eq!(1e-15, x2, 0f64 );
      assert_float_eq!(1e-15, x3, 1f64 );
    },
    _ => { assert!(false); }
  }
}
