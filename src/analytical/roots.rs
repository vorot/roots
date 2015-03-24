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

use super::super::FloatWithConstants;

#[derive(Debug,PartialEq)]
pub enum Roots<F:FloatWithConstants> {
  No([F;0]),
  One([F;1]),
  Two([F;2]),
  Three([F;3]),
  Four([F;4]),
}

impl<F:FloatWithConstants> AsRef<[F]> for Roots<F> {
  fn as_ref(&self) -> &[F] {
    match self {
      &Roots::No(ref x) => x,
      &Roots::One(ref x) => x,
      &Roots::Two(ref x) => x,
      &Roots::Three(ref x) => x,
      &Roots::Four(ref x) => x,
    }
  }
}

impl<F:FloatWithConstants> Roots<F> {
  pub fn add_sorted(self, x:F) -> Self {
    match self {
      Roots::No([]) => Roots::One([x]),
      Roots::One([x1]) if x==x1 => self,
      Roots::One([x1]) if x<x1 => Roots::Two([x,x1]),
      Roots::One([x1]) => Roots::Two([x1,x]),
      Roots::Two([x1,x2]) if (x==x1) || (x==x2) => self,
      Roots::Two([x1,x2]) if x<x1 => Roots::Three([x,x1,x2]),
      Roots::Two([x1,x2]) if x<x2 => Roots::Three([x1,x,x2]),
      Roots::Two([x1,x2]) => Roots::Three([x1,x2,x]),
      Roots::Three([x1,x2,x3]) if (x==x1) || (x==x2) || (x==x3) => self,
      Roots::Three([x1,x2,x3]) if x<x1 => Roots::Four([x,x1,x2,x3]),
      Roots::Three([x1,x2,x3]) if x<x2 => Roots::Four([x1,x,x2,x3]),
      Roots::Three([x1,x2,x3]) if x<x3 => Roots::Four([x1,x2,x,x3]),
      Roots::Three([x1,x2,x3]) => Roots::Four([x1,x2,x3,x]),
      Roots::Four([x1,x2,x3,x4]) if (x==x1) || (x==x2) || (x==x3) || (x==x4) => self,
      Roots::Four(_) => panic!("Full"),
    }
  }
}
