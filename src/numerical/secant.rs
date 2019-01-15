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
use super::Convergency;
use super::SearchError;

/// Find a root of the function f(x) = 0 using the secant method.
///
/// Pro
///
/// + Simple
/// + No need for initial bracketing
/// + No need for derivative function
///
/// Contra
///
/// - Impossible to predict which root will be found when many roots exist
/// - Unstable convergency for non-trivial functions
/// - Cannot continue when two consecutive iterations have the same value
///
/// # Failures
/// ## ZeroDerivative
/// Two consecutive points have the same value. Algorithm cannot continue.
/// ## NoConvergency
/// Algorithm cannot find a root within the given number of iterations.
/// # Examples
///
/// ```
/// use roots::SimpleConvergency;
/// use roots::find_root_secant;
///
/// let f = |x| { 1f64*x*x - 1f64 };
/// let mut convergency = SimpleConvergency { eps:1e-15f64, max_iter:30 };
///
/// let root1 = find_root_secant(10f64, 0f64, &f, &mut convergency);
/// // Returns approximately Ok(1);
///
/// let root2 = find_root_secant(-10f64, 0f64, &f, &mut 1e-15f64);
/// // Returns approximately Ok(-1);
/// ```
pub fn find_root_secant<F, Func>(first: F, second: F, f: Func, convergency: &mut Convergency<F>) -> (Result<F, SearchError>)
where
    F: FloatType,
    Func: Fn(F) -> F,
{
    let mut x1 = first;
    let mut y1 = f(x1);
    if convergency.is_root_found(y1) {
        return Ok(x1);
    }
    let mut x2 = second;
    let mut y2 = f(x2);
    if convergency.is_root_found(y2) {
        return Ok(x2);
    }

    let mut iter = 0;
    loop {
        if convergency.is_root_found(y1 - y2) {
            return Err(SearchError::ZeroDerivative);
        }
        let x = x2 - y2 * (x2 - x1) / (y2 - y1);
        if convergency.is_converged(x, x2) {
            return Ok(x);
        }
        let y = f(x);
        if convergency.is_root_found(y) {
            return Ok(x);
        }

        x1 = x2;
        y1 = y2;
        x2 = x;
        y2 = y;

        iter = iter + 1;
        if convergency.is_iteration_limit_reached(iter) {
            return Err(SearchError::NoConvergency);
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::*;
    use super::*;

    #[test]
    fn test_find_root_secant() {
        let f = |x| 1f64 * x * x - 1f64;
        let mut conv = debug_convergency::DebugConvergency::new(1e-15f64, 30);

        conv.reset();
        assert_float_eq!(1e-15f64, find_root_secant(10f64, 0f64, &f, &mut conv).ok().unwrap(), 1f64);
        assert_eq!(12, conv.get_iter_count());

        conv.reset();
        assert_float_eq!(1e-15f64, find_root_secant(-10f64, 0f64, &f, &mut conv).ok().unwrap(), -1f64);
        assert_eq!(12, conv.get_iter_count());

        conv.reset();
        assert_eq!(
            find_root_secant(10f64, -10f64, &f, &mut conv),
            Err(SearchError::ZeroDerivative)
        );
        assert_eq!(0, conv.get_iter_count());
    }
}
