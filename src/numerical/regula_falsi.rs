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

/// Illinois modification to the classical method
#[derive(Debug, PartialEq)]
enum Edge {
    /// Value is close to X1, reduce the Y1 weight
    EdgeX1,
    /// Value is in the middle of the interval
    NoEdge,
    /// Value is close to X2, reduce the Y2 weight
    EdgeX2,
}

/// Find a root of the function f(x) = 0 using the Illinois modification of the regula falsi method.
///
/// Pro
///
/// + Simple
/// + Robust
/// + No need for derivative function
///
/// Contra
///
/// - Slow
/// - Needs initial bracketing
///
/// # Failures
/// ## NoBracketing
/// Initial values do not bracket the root.
/// ## NoConvergency
/// Algorithm cannot find a root within the given number of iterations.
/// # Examples
///
/// ```
/// use roots::SimpleConvergency;
/// use roots::find_root_regula_falsi;
///
/// let f = |x| { 1f64*x*x - 1f64 };
/// let mut convergency = SimpleConvergency { eps:1e-15f64, max_iter:30 };
///
/// let root1 = find_root_regula_falsi(10f64, 0f64, &f, &mut convergency);
/// // Returns approximately Ok(1);
///
/// let root2 = find_root_regula_falsi(-10f64, 0f64, &f, &mut 1e-15f64);
/// // Returns approximately Ok(-1);
/// ```
pub fn find_root_regula_falsi<F, Func>(a: F, b: F, f: Func, convergency: &mut Convergency<F>) -> Result<F, SearchError>
where
    F: FloatType,
    Func: Fn(F) -> F,
{
    let (mut x1, mut x2) = if a > b { (b, a) } else { (a, b) };
    let mut y1 = f(x1);
    if convergency.is_root_found(y1) {
        return Ok(x1);
    }
    let mut y2 = f(x2);
    if convergency.is_root_found(y2) {
        return Ok(x2);
    }
    if y1 * y2 > F::zero() {
        return Err(SearchError::NoBracketing);
    }
    let mut edge = Edge::NoEdge;
    let mut iter = 0;
    loop {
        let x = (x1 * y2 - x2 * y1) / (y2 - y1);
        if convergency.is_converged(x1, x2) {
            return Ok(x);
        }
        let y = f(x);
        if convergency.is_root_found(y) {
            return Ok(x);
        }

        if y * y1 > F::zero() {
            x1 = x;
            y1 = y;
            if edge == Edge::EdgeX1 {
                y2 = y2 / F::two();
            }
            edge = Edge::EdgeX1;
        } else if y * y2 > F::zero() {
            x2 = x;
            y2 = y;
            if edge == Edge::EdgeX2 {
                y1 = y1 / F::two();
            }
            edge = Edge::EdgeX2;
        } else {
            return Ok(x);
        }

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
    fn test_find_root_regula_falsi() {
        let f = |x| 1f64 * x * x - 1f64;
        let mut conv = debug_convergency::DebugConvergency::new(1e-15f64, 30);

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_regula_falsi(10f64, 0f64, &f, &mut conv).ok().unwrap(),
            1f64
        );
        assert_eq!(11, conv.get_iter_count());

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_regula_falsi(-10f64, 0f64, &f, &mut conv).ok().unwrap(),
            -1f64
        );
        assert_eq!(11, conv.get_iter_count());

        conv.reset();
        assert_eq!(
            find_root_regula_falsi(10f64, 20f64, &f, &mut conv),
            Err(SearchError::NoBracketing)
        );
        let result = find_root_regula_falsi(10f64, 20f64, &f, &mut conv);
        assert_eq!(result.unwrap_err().description(), "Initial values do not bracket zero");
        assert_eq!(0, conv.get_iter_count());
    }
}
