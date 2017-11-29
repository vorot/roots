// Copyright (c) 2017, Mikhail Vorotilov
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
use super::SearchError;
use super::Convergency;
use super::super::find_roots_linear;
use super::super::find_roots_quadratic;
use super::super::find_roots_cubic;
use super::super::find_roots_quartic;
use std::fmt::Debug;

trait Polynom<F> {
    fn value(self, x: F) -> F;
    fn derivative_polynom(self) -> Vec<F>;
}

impl<'a, F> Polynom<F> for &'a [F]
    where F: FloatType + Debug
{
    fn value(self, x: F) -> F {
        let mut result = F::zero();
        let mut xn = F::one();

        // Sum starting with a0
        for a in self.iter().rev() {
            result = result + *a * xn;
            xn = xn * x;
        }

        // The highest coefficient of the normalized polynom is 1
        result + xn
    }

    fn derivative_polynom(self) -> Vec<F> {
        let mut result = Vec::from(self);
        result.truncate(self.len() - 1);
        let n: F = F::from(self.len() as i16);
        let mut ni = F::one();

        for x in result.iter_mut().rev() {
            println!("{:?} {:?}", ni, *x);
            *x = (*x * ni) / n;
            ni = ni + F::one();
        }

        result
    }
}

/// Find all roots of the normalized polynom
/// 1*x^n + a[n-1]*x^(n-1) + a[n-2]*x^(n-2) + ... + a[0] = 0.
///
/// # Failures
/// ## ZeroDerivative
/// Two consecutive points have the same value. Algorithm cannot continue.
/// ## NoConvergency
/// Algorithm cannot find a root within the given number of iterations.
/// # Examples
///
/// ```
/// use roots::find_roots_sturm;
///
/// let polynom = &[1f64,1f64,1f64,1f64,1f64,1f64];
///
/// let roots = find_roots_sturm(polynom, &mut 1e-6);
/// // Returns vector of roots;
/// ```
pub fn find_roots_sturm<F>(a: &[F], convergency: &mut Convergency<F>) -> Result<Vec<F>, SearchError>
    where F: FloatType
{
    match a.len() {
        0 => Err(SearchError::ZeroDerivative),
        1 => Ok(find_roots_linear(F::one(), a[0]).as_ref().to_vec()),
        2 => Ok(find_roots_quadratic(F::one(), a[0], a[1]).as_ref().to_vec()),
        3 => Ok(find_roots_cubic(F::one(), a[0], a[1], a[2]).as_ref().to_vec()),
        4 => Ok(find_roots_quartic(F::one(), a[0], a[1], a[2], a[3]).as_ref().to_vec()),
        _ => Err(SearchError::NoConvergency),
    }

}

#[cfg(test)]
mod test {
    use super::*;
    use super::super::*;

    #[test]
    fn test_find_roots_sturm() {
        let polynom = &[-2f64, 1f64];
        let roots = find_roots_sturm(polynom, &mut 1e-6f64).ok().unwrap();
        assert_float_array_eq!(1e-15, roots, [1f64]);
    }

    #[test]
    fn test_polynom_value() {
        let polynom = [1f64, -2f64, 1f64];
        assert_eq!(1f64, polynom.value(0f64));
        assert_eq!(1f64, polynom.value(1f64));
        assert_eq!(3f64, polynom.value(-1f64));
    }

    #[test]
    fn test_derivative_polynom() {
        let polynom = [1f64, -2f64, 1f64];
        let derivative = polynom.derivative_polynom();
        assert_float_array_eq!(1e-15, derivative, [2f64 / 3f64, -2f64 / 3f64]);
    }
}
