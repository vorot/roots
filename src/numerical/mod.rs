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
/// Possible errors
#[derive(Debug,PartialEq)]
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
        match self{
            SearchError::NoConvergency=>write!(f, "Convergency Error"),
            SearchError::NoBracketing=>write!(f, "Bracketing Error"),
            SearchError::ZeroDerivative=>write!(f, "Zero Derivative Error"),
        }
    }
}
impl Error for SearchError {
    fn description(&self) -> &str {
        match self{
            SearchError::NoConvergency=>"The algorithm could not converge within the given number of iterations",
            SearchError::NoBracketing=>"Initial values do not bracket zero",
            SearchError::ZeroDerivative=>"The algorithm cannot continue from the point where the derivative is zero",
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

pub mod newton_raphson;
pub mod brent;
pub mod secant;
pub mod regula_falsi;
pub mod polynom;

pub mod simple_convergency;
pub mod debug_convergency;
