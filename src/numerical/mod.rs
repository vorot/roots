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

use super::FloatType;
use std::error::Error;
use std::fmt;

/// Pair of the independent variable x and the function value y=F(x)
#[derive(Debug, PartialEq)]
pub struct Sample<F>
where
    F: FloatType,
{
    /// Value of the independent variable (X-axis)
    x: F,
    /// Value of the dependent variable (Y-axis)
    y: F,
}

impl<F> Sample<F>
where
    F: FloatType,
{
    fn is_bracketed_with(&self, other: &Self) -> bool {
        self.y * other.y <= F::zero()
    }
}

/// Interval between two samples, including these samples
#[derive(Debug, PartialEq)]
pub struct Interval<F>
where
    F: FloatType,
{
    /// First sample
    begin: Sample<F>,
    /// Last sample
    end: Sample<F>,
}

impl<F> Interval<F>
where
    F: FloatType,
{
    fn is_bracketed(&self) -> bool {
        self.begin.is_bracketed_with(&self.end)
    }
    fn is_converged(&self, convergency: &mut Convergency<F>) -> bool {
        convergency.is_converged(self.begin.x, self.end.x)
    }
    /// Check if the given X is inside the interval
    fn contains_x(&self, x: &F) -> bool {
        *x <= self.end.x && *x >= self.begin.x
    }
    /// Returns a point somewhere in middle of the interval for narrowing this interval down.
    /// Rules are as follows:
    /// * If the interval is bracketed, use the secant to find the middle point.
    /// ** The middle point may not be too close to either range of the interval.
    /// * If the interval is not bracketed (why would one use an unbracketed interval?), bisect it.
    fn middle(&self) -> F {
        if self.is_bracketed() && self.begin.y != self.end.y {
            let mut shift = -self.begin.y * (self.end.x - self.begin.x) / (self.end.y - self.begin.y);
            if shift < (self.end.x - self.begin.x) / F::twenty_seven() {
                shift = (self.end.x - self.begin.x) / F::twenty_seven();
            }
            if shift > (self.end.x - self.begin.x) * (F::twenty_seven() - F::one()) / F::twenty_seven() {
                shift = (self.end.x - self.begin.x) * (F::twenty_seven() - F::one()) / F::twenty_seven();
            }
            self.begin.x + shift
        } else {
            (self.begin.x + self.end.x) / F::two()
        }
    }
}

/// Possible errors
#[derive(Debug, PartialEq)]
pub enum SearchError {
    /// The algorithm could not converge within the given number of iterations
    NoConvergency,
    /// Initial values do not bracket zero
    NoBracketing,
    /// The algorithm cannot continue from the point where the derivative is zero
    ZeroDerivative,
}

impl fmt::Display for SearchError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SearchError::NoConvergency => write!(f, "Convergency Error"),
            SearchError::NoBracketing => write!(f, "Bracketing Error"),
            SearchError::ZeroDerivative => write!(f, "Zero Derivative Error"),
        }
    }
}
impl Error for SearchError {
    fn description(&self) -> &str {
        match self {
            SearchError::NoConvergency => "The algorithm could not converge within the given number of iterations",
            SearchError::NoBracketing => "Initial values do not bracket zero",
            SearchError::ZeroDerivative => "The algorithm cannot continue from the point where the derivative is zero",
        }
    }
}

/// The way to check if the algorithm has finished by either finding a root
/// or reaching the iteration limit.
pub trait Convergency<F: FloatType> {
    /// Return true if the given Y value is close enough to the zero
    fn is_root_found(&mut self, y: F) -> bool;
    /// Return true if given x values are close enough to each other
    fn is_converged(&mut self, x1: F, x2: F) -> bool;
    /// Return true if no more iterations desired
    fn is_iteration_limit_reached(&mut self, iter: usize) -> bool;
}

impl<F: FloatType> Convergency<F> for F {
    /// Return true if the given Y value is close enough to the zero
    fn is_root_found(&mut self, y: F) -> bool {
        y.abs() < self.abs()
    }
    /// Return true if given x values are close enough to each other
    fn is_converged(&mut self, x1: F, x2: F) -> bool {
        (x1 - x2).abs() < self.abs()
    }
    /// Return true if no more iterations desired
    fn is_iteration_limit_reached(&mut self, iter: usize) -> bool {
        iter >= 30
    }
}

pub mod brent;
pub mod eigen;
pub mod inverse_quadratic;
pub mod newton_raphson;
pub mod polynom;
pub mod regula_falsi;
pub mod secant;

pub mod debug_convergency;
pub mod simple_convergency;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn sample_bracketed() {
        let sample1 = Sample { x: 0f64, y: 0f64 };
        let sample2 = Sample { x: 1f64, y: 1f64 };
        let sample3 = Sample { x: 1f64, y: -1f64 };
        let sample4 = Sample { x: -1f64, y: 0f64 };
        let sample5 = Sample { x: -1f64, y: 1f64 };
        assert_eq!(true, sample1.is_bracketed_with(&sample2));
        assert_eq!(true, sample1.is_bracketed_with(&sample3));
        assert_eq!(true, sample1.is_bracketed_with(&sample4));
        assert_eq!(true, sample1.is_bracketed_with(&sample5));
        assert_eq!(true, sample2.is_bracketed_with(&sample3));
        assert_eq!(true, sample2.is_bracketed_with(&sample4));
        assert_eq!(false, sample2.is_bracketed_with(&sample5));
        assert_eq!(true, sample3.is_bracketed_with(&sample4));
        assert_eq!(true, sample3.is_bracketed_with(&sample5));
        assert_eq!(true, sample4.is_bracketed_with(&sample5));
    }

    #[test]
    fn root_interval_bracketed() {
        let sut1 = Interval {
            begin: Sample { x: 0f64, y: 0f64 },
            end: Sample { x: 0f64, y: 0f64 },
        };
        let sut2 = Interval {
            begin: Sample { x: 0f32, y: 0f32 },
            end: Sample { x: 1f32, y: 0f32 },
        };
        let sut3 = Interval {
            begin: Sample { x: 0f64, y: 0f64 },
            end: Sample { x: 0f64, y: 1f64 },
        };
        let sut4 = Interval {
            begin: Sample { x: -1f64, y: 0f64 },
            end: Sample { x: 0f64, y: 0f64 },
        };
        let sut5 = Interval {
            begin: Sample { x: -1f64, y: 0f64 },
            end: Sample { x: 0f64, y: 1f64 },
        };
        let sut6 = Interval {
            begin: Sample { x: -1f32, y: -1f32 },
            end: Sample { x: 0f32, y: 1f32 },
        };
        let sut7 = Interval {
            begin: Sample { x: 0f64, y: 1f64 },
            end: Sample { x: 1f64, y: -1f64 },
        };
        assert_eq!(true, sut1.is_bracketed());
        assert_eq!(true, sut2.is_bracketed());
        assert_eq!(true, sut3.is_bracketed());
        assert_eq!(true, sut4.is_bracketed());
        assert_eq!(true, sut5.is_bracketed());
        assert_eq!(true, sut6.is_bracketed());
        assert_eq!(true, sut7.is_bracketed());
    }

    #[test]
    fn root_interval_not_bracketed() {
        let sut1 = Interval {
            begin: Sample { x: 0f64, y: 1f64 },
            end: Sample { x: 1f64, y: 1f64 },
        };
        let sut2 = Interval {
            begin: Sample { x: -1f64, y: -1f64 },
            end: Sample { x: 1f64, y: -1f64 },
        };
        assert_eq!(false, sut1.is_bracketed());
        assert_eq!(false, sut2.is_bracketed());
    }

    #[test]
    fn root_interval_middle() {
        let sut1 = Interval {
            begin: Sample { x: 0f64, y: 1f64 },
            end: Sample { x: 2f64, y: -3f64 },
        };
        let sut2 = Interval {
            begin: Sample { x: -1f64, y: 0f64 },
            end: Sample { x: 1f64, y: 0f64 },
        };
        assert_eq!(0.5f64, sut1.middle());
        assert_eq!(0f64, sut2.middle());
    }

}
