// Copyright (c) 2019, Mikhail Vorotilov
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

use super::super::find_roots_quadratic;
use super::super::FloatType;
use super::Convergency;
use super::Interval;
use super::Sample;
use super::SearchError;

/// Definition of the quadratic equation a*x^2 + b*x + c
#[derive(Debug, PartialEq)]
pub struct Parabola<F>
where
    F: FloatType,
{
    a: F,
    b: F,
    c: F,
}

impl<F> Parabola<F>
where
    F: FloatType,
{
    /// Restore coefficients of the quadratic equation by 3 points
    fn from_three_points(p1: &Sample<F>, p2: &Sample<F>, p3: &Sample<F>) -> Self {
        let denom = (p1.x - p2.x) * (p1.x - p3.x) * (p2.x - p3.x);
        let a = (p3.x * (p2.y - p1.y) + p2.x * (p1.y - p3.y) + p1.x * (p3.y - p2.y)) / denom;
        let b = (p1.x * p1.x * (p2.y - p3.y) + p3.x * p3.x * (p1.y - p2.y) + p2.x * p2.x * (p3.y - p1.y)) / denom;
        let c = (p2.x * p2.x * (p3.x * p1.y - p1.x * p3.y)
            + p2.x * (p1.x * p1.x * p3.y - p3.x * p3.x * p1.y)
            + p1.x * p3.x * (p3.x - p1.x) * p2.y)
            / denom;

        Parabola { a: a, b: b, c: c }
    }
}

/// Find a root of the function f(x) = 0 using inverse quadratic approximation.
///
/// Pro
///
/// + Faster than linear approximation
/// + No need for derivative function
///
/// Contra
///
/// - sqrt is calculated on every step
/// - only works for polynomial-like functions
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
/// use roots::find_root_inverse_quadratic;
///
/// let f = |x| { 1f64*x*x - 1f64 };
/// let mut convergency = SimpleConvergency { eps:1e-15f64, max_iter:30 };
///
/// let root1 = find_root_inverse_quadratic(10f64, 0f64, &f, &mut convergency);
/// // Returns approximately Ok(1);
///
/// let root2 = find_root_inverse_quadratic(-10f64, 0f64, &f, &mut 1e-15f64);
/// // Returns approximately Ok(-1);
/// ```
pub fn find_root_inverse_quadratic<F, Func>(a: F, b: F, f: Func, convergency: &mut Convergency<F>) -> Result<F, SearchError>
where
    F: FloatType,
    Func: Fn(F) -> F,
{
    let (x1, x2) = if a > b { (b, a) } else { (a, b) };
    let sample1 = Sample { x: x1, y: f(x1) };
    if convergency.is_root_found(sample1.y) {
        return Ok(sample1.x);
    }
    let sample2 = Sample { x: x2, y: f(x2) };
    if convergency.is_root_found(sample2.y) {
        return Ok(sample2.x);
    }
    if !sample1.is_bracketed_with(&sample2) {
        return Err(SearchError::NoBracketing);
    }

    // Initially, find x3 using the regula falsi method
    let mut interval = Interval {
        begin: sample1,
        end: sample2,
    };
    let mut x3 = interval.middle();
    if interval.is_converged(convergency) {
        return Ok(x3);
    }
    let mut sample3 = Sample { x: x3, y: f(x3) };
    if convergency.is_root_found(sample3.y) {
        return Ok(sample3.x);
    }

    // Iterate quadratically
    let mut iter = 0;
    loop {
        println! {"interval:{:?}, middle:{:?}",interval,sample3};
        let parabola = Parabola::from_three_points(&interval.begin, &interval.end, &sample3);
        println! {"parabola:{:?}",parabola};

        // Find the new approximation quadratically
        x3 = if let Some(root) = find_roots_quadratic(parabola.a, parabola.b, parabola.c)
            .as_ref()
            .iter()
            .find(|x| interval.contains_x(x))
        {
            println! {"root:{:?}",root};
            *root
        } else {
            // no roots inside interval, fallback to linear approximation
            println! {"fallback middle:{:?}",interval.middle()};
            interval.middle()
        };

        // calculate the approximated value
        sample3 = Sample { x: x3, y: f(x3) };

        if convergency.is_root_found(sample3.y) {
            return Ok(sample3.x);
        }

        // Narrow down the search interval while keeping the root bracketed
        if sample3.is_bracketed_with(&interval.begin) {
            interval.end = Sample {
                x: sample3.x,
                y: sample3.y,
            };
        } else {
            interval.begin = Sample {
                x: sample3.x,
                y: sample3.y,
            };
        }

        if interval.is_converged(convergency) {
            return Ok(interval.middle());
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
    fn test_find_root_inverse_quadratic() {
        let f = |x| 1f64 * x * x - 1f64;
        let mut conv = debug_convergency::DebugConvergency::new(1e-15f64, 30);

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_inverse_quadratic(10f64, 0f64, &f, &mut conv).ok().unwrap(),
            1f64
        );
        assert_eq!(0, conv.get_iter_count());

        conv.reset();
        assert_float_eq!(
            1e-15f64,
            find_root_inverse_quadratic(-10f64, 0f64, &f, &mut conv).ok().unwrap(),
            -1f64
        );
        assert_eq!(0, conv.get_iter_count());

        conv.reset();
        assert_eq!(
            find_root_inverse_quadratic(10f64, 20f64, &f, &mut conv),
            Err(SearchError::NoBracketing)
        );
        let result = find_root_inverse_quadratic(10f64, 20f64, &f, &mut conv);
        assert_eq!(result.unwrap_err().description(), "Initial values do not bracket zero");
        assert_eq!(0, conv.get_iter_count());
    }

    #[test]
    fn test_from_three_points() {
        assert_eq!(
            Parabola {
                a: 1f64,
                b: 0f64,
                c: -1f64
            },
            Parabola::from_three_points(
                &Sample { x: -10f64, y: 99f64 },
                &Sample { x: -2f64, y: 3f64 },
                &Sample { x: 0f64, y: -1f64 }
            )
        );
        assert_eq!(
            Parabola {
                a: 1f64,
                b: 0f64,
                c: -1f64
            },
            Parabola::from_three_points(
                &Sample { x: 10f64, y: 99f64 },
                &Sample { x: 2f64, y: 3f64 },
                &Sample { x: 0f64, y: -1f64 }
            )
        );
        assert_eq!(
            Parabola {
                a: 1f64,
                b: 0f64,
                c: -1f64
            },
            Parabola::from_three_points(
                &Sample { x: -3f64, y: 8f64 },
                &Sample { x: 2f64, y: 3f64 },
                &Sample { x: 0f64, y: -1f64 }
            )
        );
    }
}
