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

use super::super::find_roots_cubic;
use super::super::find_roots_linear;
use super::super::find_roots_quadratic;
use super::super::find_roots_quartic;
use super::super::FloatType;
use super::Convergency;
use super::Interval;
use super::Sample;
use super::SearchError;

#[derive(Debug, PartialEq)]
struct ValueAndDerivative<F>
where
    F: FloatType,
{
    value: Sample<F>,
    derivative: F,
}

trait Polynom<F>
where
    F: FloatType,
{
    fn value(&self, x: &F) -> F;
    fn value_and_derivative(&self, x: &F) -> ValueAndDerivative<F>;
    fn find_root(&self, bracketed_start: &mut Interval<F>, convergency: &mut Convergency<F>) -> Result<F, SearchError>;
    fn derivative_polynom(&self) -> Vec<F>;
}

impl<F> Polynom<F> for [F]
where
    F: FloatType,
{
    fn value(&self, x: &F) -> F {
        let mut result = F::zero();
        let mut xn = F::one();

        // Sum starting with a0
        for a in self.iter().rev() {
            result = result + *a * xn;
            xn = xn * *x;
        }

        // The highest coefficient of the normalized polynom is 1
        result + xn
    }

    fn value_and_derivative(&self, x: &F) -> ValueAndDerivative<F> {
        let mut xn = F::one(); // x^n for SUM(A(n)*x^(n))
        let mut value = F::zero();

        let mut xn1 = F::zero(); // x^n-1 for SUM(n*A(n-1)*x^(n-1))
        let mut derivative = F::zero();
        let mut n = F::zero();

        // Sum starting with a0
        for a in self.iter().rev() {
            value = value + *a * xn;
            derivative = derivative + *a * n * xn1;
            xn1 = xn;
            xn = xn * *x;
            n = n + F::one();
        }

        // The highest coefficient of the normalized polynom is 1
        ValueAndDerivative {
            value: Sample { x: *x, y: value + xn },
            derivative: derivative + n * xn1,
        }
    }

    fn find_root(&self, bracketed_start: &mut Interval<F>, convergency: &mut Convergency<F>) -> Result<F, SearchError> {
        if bracketed_start.is_bracketed() {
            let interval = bracketed_start;
            let mut iter = 0;
            loop {
                if convergency.is_root_found(interval.begin.y) {
                    break Ok(interval.begin.x);
                } else if convergency.is_root_found(interval.end.y) {
                    break Ok(interval.end.x);
                } else if interval.is_converged(convergency) {
                    break Ok(interval.middle());
                } else {
                    let middle = self.value_and_derivative(&interval.middle());
                    let next_sample = if middle.derivative != F::zero() {
                        let newton_raphson = middle.value.x - middle.value.y / middle.derivative;
                        if newton_raphson >= interval.begin.x && newton_raphson <= interval.end.x {
                            let newton_raphson_value = self.value(&newton_raphson);
                            if newton_raphson_value.abs() < middle.value.y.abs() {
                                Sample {
                                    x: newton_raphson,
                                    y: newton_raphson_value,
                                }
                            } else {
                                middle.value
                            }
                        } else {
                            middle.value
                        }
                    } else {
                        middle.value
                    };
                    if interval.begin.is_bracketed_with(&next_sample) {
                        interval.end = Sample {
                            x: next_sample.x,
                            y: next_sample.y,
                        };
                    } else {
                        interval.begin = Sample {
                            x: next_sample.x,
                            y: next_sample.y,
                        };
                    }
                }
                iter = iter + 1;
                if convergency.is_iteration_limit_reached(iter) {
                    break Err(SearchError::NoConvergency);
                }
            }
        } else {
            Err(SearchError::NoBracketing)
        }
    }

    fn derivative_polynom(&self) -> Vec<F> {
        let mut result = Vec::from(self);
        result.truncate(self.len() - 1);
        let n: F = F::from(self.len() as i16);
        let mut ni = F::one();

        for x in result.iter_mut().rev() {
            *x = (*x * ni) / n;
            ni = ni + F::one();
        }

        result
    }
}

/// Interval for searching roots
enum SearchInterval<F>
where
    F: FloatType,
{
    /// [-infinity .. +infinity]
    Whole,
    /// [-infinity .. x]
    First(Sample<F>),
    /// [x .. +infinity ]
    Last(Sample<F>),
    /// [x1 .. x2 ]
    Middle(Interval<F>),
}

enum BracketingDirection {
    TowardsPositive,
    TowardsNegative,
}

fn initial_bracket<F>(
    initial_sample: &Sample<F>,
    direction: &BracketingDirection,
    polynom: &[F],
    derivative_polynom: &[F],
    convergency: &mut Convergency<F>,
) -> Result<Interval<F>, SearchError>
where
    F: FloatType,
{
    let mut iter = 0;
    let towards_positive = match direction {
        &BracketingDirection::TowardsPositive => true,
        &BracketingDirection::TowardsNegative => false,
    };
    let mut step = if towards_positive { F::one() } else { -F::one() };
    let initial_copy = Sample {
        x: initial_sample.x,
        y: initial_sample.y,
    };
    let mut next_x = initial_sample.x + step;
    let result = loop {
        let mut next_y = polynom.value(&next_x);
        let mut next_sample = Sample { x: next_x, y: next_y };
        if next_sample.is_bracketed_with(&initial_sample) {
            break Ok(if towards_positive {
                Interval {
                    begin: initial_copy,
                    end: next_sample,
                }
            } else {
                Interval {
                    begin: next_sample,
                    end: initial_copy,
                }
            });
        } else {
            let derivative = derivative_polynom.value(&next_x);
            if derivative > F::zero() {
                next_x = next_x - next_y / derivative;
                next_y = polynom.value(&next_x);
                next_sample = Sample { x: next_x, y: next_y };
                if next_sample.is_bracketed_with(&initial_sample) {
                    break Ok(if towards_positive {
                        Interval {
                            begin: initial_copy,
                            end: next_sample,
                        }
                    } else {
                        Interval {
                            begin: next_sample,
                            end: initial_copy,
                        }
                    });
                }
            };
            step = step * F::two();
            next_x = next_x + step;
            iter = iter + 1;
            if convergency.is_iteration_limit_reached(iter) {
                break Err(SearchError::NoConvergency);
            };
        }
    };
    result
}

fn narrow_down<F>(
    initial_interval: &SearchInterval<F>,
    polynom: &[F],
    derivative_polynom: &[F],
    convergency: &mut Convergency<F>,
) -> Result<Interval<F>, SearchError>
where
    F: FloatType,
{
    match initial_interval {
        &SearchInterval::Whole => {
            let zero_sample = Sample {
                x: F::zero(),
                y: polynom.value(&F::zero()),
            };
            let zero_interval = if zero_sample.y > F::zero() {
                SearchInterval::First(zero_sample)
            } else {
                SearchInterval::Last(zero_sample)
            };
            narrow_down(&zero_interval, polynom, derivative_polynom, convergency)
        }
        &SearchInterval::First(ref end) => initial_bracket(
            &end,
            &BracketingDirection::TowardsNegative,
            polynom,
            derivative_polynom,
            convergency,
        ),
        &SearchInterval::Last(ref begin) => initial_bracket(
            &begin,
            &BracketingDirection::TowardsPositive,
            polynom,
            derivative_polynom,
            convergency,
        ),
        &SearchInterval::Middle(ref interval) => {
            if interval.is_bracketed() {
                let middle_x = if interval.begin.y == interval.end.y {
                    (interval.begin.x + interval.end.x) / F::two()
                } else {
                    interval.begin.x - interval.begin.y * (interval.end.x - interval.begin.x) / (interval.end.y - interval.begin.y)
                };
                let mut middle_sample = Sample {
                    x: middle_x,
                    y: polynom.value(&middle_x),
                };
                let derivative = derivative_polynom.value(&middle_x);
                if derivative != F::zero() {
                    let closer_x = middle_sample.x - middle_sample.y / derivative;
                    if closer_x >= interval.begin.x && closer_x <= interval.end.x {
                        middle_sample = Sample {
                            x: closer_x,
                            y: polynom.value(&closer_x),
                        };
                    }
                }
                if interval.begin.is_bracketed_with(&middle_sample) {
                    Ok(Interval {
                        begin: Sample {
                            x: interval.begin.x,
                            y: interval.begin.y,
                        },
                        end: middle_sample,
                    })
                } else {
                    Ok(Interval {
                        begin: middle_sample,
                        end: Sample {
                            x: interval.end.x,
                            y: interval.end.y,
                        },
                    })
                }
            } else {
                Err(SearchError::NoBracketing)
            }
        }
    }
}

fn find_root_intervals<F>(
    polynom: &[F],
    derivative_polynom: &[F],
    convergency: &mut Convergency<F>,
) -> Result<Vec<SearchInterval<F>>, SearchError>
where
    F: FloatType,
{
    let mut result = Vec::new();
    let derivative_roots = find_roots_sturm(&derivative_polynom, convergency);
    let symmetric_polynom = polynom.len() % 2 == 0;
    let mut expect_positive = !symmetric_polynom;
    let mut previous_interval: SearchInterval<F> = SearchInterval::Whole;
    // Iterate through all roots of the derivative polynom
    for derivative_root in derivative_roots.iter().filter_map(|s| match s {
        &Ok(ref x) => Some(x),
        &Err(_) => None,
    }) {
        let value = polynom.value(derivative_root);
        if (expect_positive && value >= F::zero()) || (!expect_positive && value < F::zero()) {
            // Transition found
            let interval_to_add = match &previous_interval {
                &SearchInterval::Whole => SearchInterval::First(Sample {
                    x: *derivative_root,
                    y: value,
                }),
                &SearchInterval::First(ref previous_end) => SearchInterval::Middle(Interval {
                    begin: Sample {
                        x: previous_end.x,
                        y: previous_end.y,
                    },
                    end: Sample {
                        x: *derivative_root,
                        y: value,
                    },
                }),
                _ => panic!("Unexpected type of the previous root interval!"),
            };
            result.push(interval_to_add);
            expect_positive = !expect_positive;
        }
        previous_interval = SearchInterval::First(Sample {
            x: *derivative_root,
            y: value,
        });
    }
    // All roots are checked, now the final step
    match previous_interval {
        SearchInterval::Whole => {
            if !symmetric_polynom {
                result.push(SearchInterval::Whole);
            }
            Ok(result)
        }
        SearchInterval::First(sample) => {
            if sample.x < F::zero() {
                result.push(SearchInterval::Last(sample));
            }
            Ok(result)
        }
        _ => Err(SearchError::NoBracketing),
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
pub fn find_roots_sturm<F>(a: &[F], convergency: &mut Convergency<F>) -> Vec<Result<F, SearchError>>
where
    F: FloatType,
{
    match a.len() {
        0 => Vec::new(),
        1 => find_roots_linear(F::one(), a[0]).as_ref().iter().map(|s| Ok(*s)).collect(),
        2 => find_roots_quadratic(F::one(), a[0], a[1])
            .as_ref()
            .iter()
            .map(|s| Ok(*s))
            .collect(),
        3 => find_roots_cubic(F::one(), a[0], a[1], a[2])
            .as_ref()
            .iter()
            .map(|s| Ok(*s))
            .collect(),
        4 => find_roots_quartic(F::one(), a[0], a[1], a[2], a[3])
            .as_ref()
            .iter()
            .map(|s| Ok(*s))
            .collect(),
        _ => {
            let mut result = Vec::new();
            let derivative_polynom = a.derivative_polynom();
            match find_root_intervals(a, &derivative_polynom, convergency) {
                Ok(root_intervals) => {
                    for root_interval in &root_intervals {
                        if let Ok(mut narrowed) = narrow_down(&root_interval, a, &derivative_polynom, convergency) {
                            result.push(a.find_root(&mut narrowed, convergency));
                        }
                    }
                }
                Err(error) => {
                    result.push(Err(error));
                }
            }
            result
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::*;
    use super::*;

    #[test]
    fn test_find_roots_sturm() {
        let polynom = &[-2f64, 1f64];
        let roots = find_roots_sturm(polynom, &mut 1e-6f64);
        assert_eq!(roots, [Ok(1f64)]);
    }

    #[test]
    fn test_polynom_value() {
        let polynom = [1f64, -2f64, 1f64];
        assert_eq!(1f64, polynom.value(&0f64));
        assert_eq!(1f64, polynom.value(&1f64));
        assert_eq!(3f64, polynom.value(&-1f64));
    }

    #[test]
    fn test_polynom_value_and_derivative() {
        let polynom = [1f64, -2f64, 1f64];
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: 0f64, y: 1f64 },
                derivative: -2f64
            },
            polynom.value_and_derivative(&0f64)
        );
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: 1f64, y: 1f64 },
                derivative: 3f64
            },
            polynom.value_and_derivative(&1f64)
        );
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: -1f64, y: 3f64 },
                derivative: -1f64
            },
            polynom.value_and_derivative(&-1f64)
        );
    }

    #[test]
    fn test_derivative_polynom_3() {
        // x^3 + 1*x^2 - 2*x^1 + 1*x^0 => 3*x^2 + 2*x^1 - 2*x^0 => x^2 + (2/3)*x^1 - (2/3)*x^0
        let polynom = [1f64, -2f64, 1f64];
        let derivative = polynom.derivative_polynom();
        assert_float_array_eq!(1e-15, derivative, [2f64 / 3f64, -2f64 / 3f64]);
    }

    #[test]
    fn test_derivative_polynom_5() {
        // x^5 - 2*x^4 - 3*x^3 + 4*x^2 + 0*x^1 + 0*x^0 => 5*x^4 - 8*x^3 - 9*x^2 + 8*x^1 + 0*x^0 => x^4 - (8/5)*x^3 - (9/5)*x^2 + (8/5)*x^1 + 0*x^0
        let polynom = [-2f64, -3f64, 4f64, 0f64, 0f64];
        let derivative = polynom.derivative_polynom();
        assert_float_array_eq!(1e-15, derivative, [-8f64 / 5f64, -9f64 / 5f64, 8f64 / 5f64, 0f64]);
    }
}
