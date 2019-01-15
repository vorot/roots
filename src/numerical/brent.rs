// Copyright 2015 Mikhail Vorotilov. See the COPYRIGHT
// file at the top-level directory of this distribution and at
// http://rust-lang.org/COPYRIGHT.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::super::FloatType;
use super::Convergency;
use super::SearchError;

/// Arrange two points so that the greatest value is first
fn arrange<F: FloatType>(a: F, ya: F, b: F, yb: F) -> (F, F, F, F) {
    if ya.abs() > yb.abs() {
        (a, ya, b, yb)
    } else {
        (b, yb, a, ya)
    }
}

/// Find a root of the function f(x) = 0 using the Brent method.
///
/// Pro
///
/// + Fast
/// + Robust
/// + No need for derivative function
///
/// Contra
///
/// - Complicated
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
/// use roots::find_root_brent;
///
/// let f = |x| { 1f64*x*x - 1f64 };
/// let mut convergency = SimpleConvergency { eps:1e-15f64, max_iter:30 };
///
/// let root1 = find_root_brent(10f64, 0f64, &f, &mut convergency);
/// // Returns approximately Ok(1);
///
/// let root2 = find_root_brent(-10f64, 0f64, &f, &mut 1e-15f64);
/// // Returns approximately Ok(-1);
/// ```
pub fn find_root_brent<F, Func>(a: F, b: F, f: Func, convergency: &mut Convergency<F>) -> Result<F, SearchError>
where
    F: FloatType,
    Func: Fn(F) -> F,
{
    let (mut a, mut ya, mut b, mut yb) = arrange(a, f(a), b, f(b));
    if ya * yb > F::zero() {
        return Err(SearchError::NoBracketing);
    }

    let (mut c, mut yc, mut d) = (a, ya, a);
    let mut flag = true;

    let mut iter = 0;
    loop {
        if convergency.is_root_found(ya) {
            return Ok(a);
        }
        if convergency.is_root_found(yb) {
            return Ok(b);
        }
        if convergency.is_converged(a, b) {
            return Ok(c);
        }
        let mut s = if (ya != yc) && (yb != yc) {
            a * yb * yc / ((ya - yb) * (ya - yc)) + b * ya * yc / ((yb - ya) * (yb - yc)) + c * ya * yb / ((yc - ya) * (yc - yb))
        } else {
            b - yb * (b - a) / (yb - ya)
        };

        let cond1 = (s - b) * (s - (F::three() * a + b) / F::four()) > F::zero();
        let cond2 = flag && (s - b).abs() >= (b - c).abs() / F::two();
        let cond3 = !flag && (s - b).abs() >= (c - d).abs() / F::two();
        let cond4 = flag && convergency.is_converged(b, c);
        let cond5 = !flag && convergency.is_converged(c, d);

        if cond1 || cond2 || cond3 || cond4 || cond5 {
            s = (a + b) / F::two();
            flag = true;
        } else {
            flag = false;
        }

        let ys = f(s);
        d = c;
        c = b;
        yc = yb;
        if ya * ys < F::zero() {
            // Root bracketed between a ans s
            match arrange(a, f(a), s, ys) {
                (_a, _ya, _b, _yb) => {
                    a = _a;
                    ya = _ya;
                    b = _b;
                    yb = _yb;
                }
            }
        } else {
            // Root bracketed between s ans b
            match arrange(s, ys, b, f(b)) {
                (_a, _ya, _b, _yb) => {
                    a = _a;
                    ya = _ya;
                    b = _b;
                    yb = _yb;
                }
            }
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
    fn test_find_root_brent() {
        let f = |x| 1f64 * x * x - 1f64;
        let mut conv = debug_convergency::DebugConvergency::new(1e-15f64, 30);

        conv.reset();
        assert_float_eq!(1e-15f64, find_root_brent(10f64, 0f64, &f, &mut conv).ok().unwrap(), 1f64);
        assert_eq!(10, conv.get_iter_count());

        conv.reset();
        assert_float_eq!(1e-15f64, find_root_brent(-10f64, 0f64, &f, &mut conv).ok().unwrap(), -1f64);
        assert_eq!(10, conv.get_iter_count());

        conv.reset();
        assert_eq!(find_root_brent(10f64, 20f64, &f, &mut conv), Err(SearchError::NoBracketing));
        assert_eq!(0, conv.get_iter_count());
    }

    #[test]
    fn test_find_root_brent_simple() {
        let f = |x| 1f64 * x * x - 1f64;

        assert_float_eq!(1e-15f64, find_root_brent(10f64, 0f64, &f, &mut 1e-15f64).ok().unwrap(), 1f64);

        assert_float_eq!(
            1e-15f64,
            find_root_brent(-10f64, 0f64, &f, &mut 1e-15f64).ok().unwrap(),
            -1f64
        );

        assert_eq!(
            find_root_brent(10f64, 20f64, &f, &mut 1e-15f64),
            Err(SearchError::NoBracketing)
        );
    }
}
