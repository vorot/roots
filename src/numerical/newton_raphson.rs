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

/// Find a root of the function f(x) = 0 using the Newton-Raphson method.
///
/// Pro
///
/// + Simple
/// + Fast convergency for well-behaved functions
/// + No need for initial bracketing
///
/// Contra
///
/// - Needs derivative function
/// - Impossible to predict which root will be found when many roots exist
/// - Unstable convergency for non-trivial functions
/// - Cannot continue when derivative is zero
///
/// # Failures
/// ## ZeroDerivative
/// The stationary point of the function is encountered. Algorithm cannot continue.
/// ## NoConvergency
/// Algorithm cannot find a root within the given number of iterations.
/// # Examples
/// ```
/// use roots::SimpleConvergency;
/// use roots::find_root_newton_raphson;
///
/// let f = |x| { 1f64*x*x - 1f64 };
/// let d = |x| { 2f64*x };
/// let mut convergency = SimpleConvergency { eps:1e-15f64, max_iter:30 };
///
/// let root1 = find_root_newton_raphson(10f64, &f, &d, &mut convergency);
/// // Returns approximately Ok(1);
///
/// let root2 = find_root_newton_raphson(-10f64, &f, &d, &mut 1e-15f64);
/// // Returns approximately Ok(-1);
/// ```
pub fn find_root_newton_raphson<F, Func, Deriv>(
    start: F,
    f: Func,
    d: Deriv,
    convergency: &mut Convergency<F>,
) -> Result<F, SearchError>
where
    F: FloatType,
    Func: Fn(F) -> F,
    Deriv: Fn(F) -> F,
{
    let mut x = start;

    let mut iter = 0;
    loop {
        let f = f(x);
        let d = d(x);
        if convergency.is_root_found(f) {
            return Ok(x);
        }
        // Derivative is 0; try to correct the bad starting point
        if convergency.is_root_found(d) {
            if iter == 0 {
                x = x + F::one();
                iter = iter + 1;
                continue;
            } else {
                return Err(SearchError::ZeroDerivative);
            }
        }

        let x1 = x - f / d;
        if convergency.is_converged(x, x1) {
            return Ok(x1);
        }

        x = x1;
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
    fn test_find_root_newton_raphson() {
        let f = |x| 1f64 * x * x - 1f64;
        let d = |x| 2f64 * x;
        let mut conv = debug_convergency::DebugConvergency::new(1e-15f64, 30);

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_newton_raphson(10f64, &f, &d, &mut conv).ok().unwrap(),
            1f64
        );
        assert_eq!(8, conv.get_iter_count());

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_newton_raphson(-10f64, &f, &d, &mut conv).ok().unwrap(),
            -1f64
        );
        assert_eq!(8, conv.get_iter_count());
    }
}
