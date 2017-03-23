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

/// Solves a linear equation a1*x + a0 = 0.
///
/// # Examples
///
/// ```
/// use roots::Roots;
/// use roots::find_roots_linear;
///
/// // Returns Roots::No([]) as '0*x + 1 = 0' has no roots;
/// let no_root = find_roots_linear(0f32, 1f32);
/// assert_eq!(no_root, Roots::No([]));
///
/// // Returns Roots::Two([0f64]) as '1*x + 0 = 0' has the root 0
/// let root = find_roots_linear(1f64, 0f64);
/// assert_eq!(root, Roots::One([0f64]));
///
/// // Returns Roots::One([0f32]) as 0 is one of roots of '0*x + 0 = 0'
/// let zero_root = find_roots_linear(0f32, 0f32);
/// assert_eq!(zero_root, Roots::One([0f32]));
/// ```
pub fn find_roots_linear<F: FloatType>(a1: F, a0: F) -> Roots<F> {
    if a1 == F::zero() {
        if a0 == F::zero() {
            Roots::One([F::zero()])
        } else {
            Roots::No([])
        }
    } else {
        Roots::One([-a0 / a1])
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_linear() {
        assert_eq!(find_roots_linear(0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_linear(2f64, 1f64), Roots::One([-0.5f64]));
        assert_eq!(find_roots_linear(0f32, 1f32), Roots::No([]));
    }
}
