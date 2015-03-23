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

use std::cell::Cell;
use std::num::Float;
use std::fmt::Display;
use super::Convergency;

/// Convergency provider for debugging.
pub struct DebugConvergency<F:Float>{
  /// Precision for both X and Y axes
  eps: F,
  /// Maximum number of iterations
  max_iter: usize,
  /// Last iteration
  iter: Cell<usize>,
}

impl<F:Float> DebugConvergency<F> {
  pub fn new(eps: F, max_iter: usize) -> DebugConvergency<F> {
    DebugConvergency { eps:eps, max_iter:max_iter,iter:Cell::new(0) }
  }

  pub fn reset(self: &DebugConvergency<F>) {
    self.iter.set(0);
  }

  pub fn get_iter_count(self: &DebugConvergency<F>) -> usize {
    self.iter.get()
  }
}

impl<F:Float+Display> Convergency<F> for DebugConvergency<F> {
  /// Prints the value being checked
  fn is_root_found(&self, y:F) -> bool { println!("#{} check root {}", self.iter.get(), y); y.abs() < self.eps.abs() }
  /// Prints values being checked
  fn is_converged(&self, x1:F, x2:F) -> bool { println!("#{} check convergency {}-{}", self.iter.get(), x1, x2); (x1-x2).abs() < self.eps.abs() }
  /// Updates internal iteration counter
  fn is_iteration_limit_reached(&self, iter:usize) -> bool { println!("#{} check iteration limit {}", self.iter.get(), iter); self.iter.set(iter); iter >= self.max_iter }
}
