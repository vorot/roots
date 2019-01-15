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

/// Sorted and unique list of roots of an equation.
#[derive(Debug, PartialEq)]
pub enum Roots<F: FloatType> {
    /// Equation has no roots
    No([F; 0]),
    /// Equation has one root (or all roots of the equation are the same)
    One([F; 1]),
    /// Equation has two roots
    Two([F; 2]),
    /// Equation has three roots
    Three([F; 3]),
    /// Equation has four roots
    Four([F; 4]),
}

impl<F: FloatType> AsRef<[F]> for Roots<F> {
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

impl<F: FloatType> Roots<F> {
    fn check_new_root(&self, new_root: F) -> (bool, usize) {
        let mut pos = 0;
        let mut exists = false;

        for x in self.as_ref().iter() {
            if *x == new_root {
                exists = true;
                break;
            }
            if *x > new_root {
                break;
            }
            pos = pos + 1;
        }

        (exists, pos)
    }

    /// Add a new root to existing ones keeping the list of roots ordered and unique.
    pub fn add_new_root(self, new_root: F) -> Self {
        match self {
            Roots::No(_) => Roots::One([new_root]),
            _ => {
                let (exists, pos) = self.check_new_root(new_root);

                if exists {
                    self
                } else {
                    let old_roots = self.as_ref();
                    match (old_roots.len(), pos) {
                        (1, 0) => Roots::Two([new_root, old_roots[0]]),
                        (1, 1) => Roots::Two([old_roots[0], new_root]),
                        (2, 0) => Roots::Three([new_root, old_roots[0], old_roots[1]]),
                        (2, 1) => Roots::Three([old_roots[0], new_root, old_roots[1]]),
                        (2, 2) => Roots::Three([old_roots[0], old_roots[1], new_root]),
                        (3, 0) => Roots::Four([new_root, old_roots[0], old_roots[1], old_roots[2]]),
                        (3, 1) => Roots::Four([old_roots[0], new_root, old_roots[1], old_roots[2]]),
                        (3, 2) => Roots::Four([old_roots[0], old_roots[1], new_root, old_roots[2]]),
                        (3, 3) => Roots::Four([old_roots[0], old_roots[1], old_roots[2], new_root]),
                        _ => panic!("Cannot add root"),
                    }
                }
            }
        }
    }
}

#[test]
fn test_roots() {
    let mut roots = Roots::One([1f32]);
    assert_eq!(roots, Roots::One([1f32]));

    roots = roots.add_new_root(1f32);
    assert_eq!(roots, Roots::One([1f32]));

    roots = roots.add_new_root(0f32);
    assert_eq!(roots, Roots::Two([0f32, 1f32]));

    roots = roots.add_new_root(0f32);
    assert_eq!(roots, Roots::Two([0f32, 1f32]));

    roots = roots.add_new_root(3f32);
    assert_eq!(roots, Roots::Three([0f32, 1f32, 3f32]));

    roots = roots.add_new_root(2f32);
    assert_eq!(roots, Roots::Four([0f32, 1f32, 2f32, 3f32]));
}
