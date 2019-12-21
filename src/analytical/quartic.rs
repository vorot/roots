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
use super::super::Roots;

/// Solves a quartic equation a4*x^4 + a4*x^3 + a2*x^2 + a1*x + a0 = 0.
/// pp, rr, and dd are already computed while searching for multiple roots
fn find_roots_via_depressed_quartic<F: FloatType>(a4: F, a3: F, a2: F, a1: F, a0: F, pp: F, rr: F, dd: F) -> Roots<F> {
    // Depressed quartic
    // https://en.wikipedia.org/wiki/Quartic_function#Converting_to_a_depressed_quartic

    let _2 = F::from(2i16);
    let _3 = F::from(3i16);
    let _4 = F::from(4i16);
    let _6 = F::from(6i16);
    let _8 = F::from(8i16);
    let _12 = F::from(12i16);
    let _16 = F::from(16i16);
    let _256 = F::from(256i16);

    // a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0 => y^4 + p*y^2 + q*y + r.
    let a4_pow_2 = a4 * a4;
    let a4_pow_3 = a4_pow_2 * a4;
    let a4_pow_4 = a4_pow_2 * a4_pow_2;
    // Re-use pre-calculated values
    let p = pp / (_8 * dbg!(a4_pow_2));
    let q = rr / (_8 * dbg!(a4_pow_3));
    let r = (dd + _16 * a4_pow_2 * (_12 * a0 * a4 - _3 * a1 * a3 + a2 * a2)) / (_256 * dbg!(a4_pow_4));

    let mut roots = Roots::No([]);
    for y in super::quartic_depressed::find_roots_quartic_depressed(dbg!(p), dbg!(q), dbg!(r))
        .as_ref()
        .iter()
    {
        roots = roots.add_new_root(dbg!(*y) - a3 / (_4 * a4));
    }
    roots
}

/// Solves a quartic equation a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0.
///
/// Returned roots are ordered.
/// Precision is about 5e-15 for f64, 5e-7 for f32.
///
/// # Examples
///
/// ```
/// use roots::find_roots_quartic;
///
/// let one_root = find_roots_quartic(1f64, 0f64, 0f64, 0f64, 0f64);
/// // Returns Roots::One([0f64]) as 'x^4 = 0' has one root 0
///
/// let two_roots = find_roots_quartic(1f32, 0f32, 0f32, 0f32, -1f32);
/// // Returns Roots::Two([-1f32, 1f32]) as 'x^4 - 1 = 0' has roots -1 and 1
/// ```
pub fn find_roots_quartic<F: FloatType>(a4: F, a3: F, a2: F, a1: F, a0: F) -> Roots<F> {
    // Handle non-standard cases
    if a4 == F::zero() {
        // a4 = 0; a3*x^3 + a2*x^2 + a1*x + a0 = 0; solve cubic equation
        super::cubic::find_roots_cubic(a3, a2, a1, a0)
    } else if a0 == F::zero() {
        // a0 = 0; x^4 + a2*x^2 + a1*x = 0; reduce to cubic and arrange results
        super::cubic::find_roots_cubic(a4, a3, a2, a1).add_new_root(F::zero())
    } else if a1 == F::zero() && a3 == F::zero() {
        // a1 = 0, a3 =0; a4*x^4 + a2*x^2 + a0 = 0; solve bi-quadratic equation
        super::biquadratic::find_roots_biquadratic(a4, a2, a0)
    } else {
        let _3 = F::three();
        let _4 = F::four();
        let _6 = F::four() + F::two();
        let _8 = F::four() + F::four();
        let _9 = _6 + _3;
        let _10 = _6 + _4;
        let _12 = _8 + _4;
        let _16 = _8 + _8;
        let _18 = _16 + F::two();
        let _27 = F::twenty_seven();
        let _64 = _8 * _8;
        let _72 = _64 + _8;
        let _80 = _64 + _16;
        let _128 = _64 + _64;
        let _144 = _80 + _64;
        let _192 = _128 + _64;
        let _256 = _192 + _64;
        // Discriminant
        // https://en.wikipedia.org/wiki/Quartic_function#Nature_of_the_roots
        let discriminant = dbg!(a4 * a0 * a4 * (_256 * a4 * a0 * a0 - _192 * dbg!(a3) * a1 * a0 + _144 * a2 * dbg!(a1) * a1))
            + dbg!(a4 * a0 * a2 * a2 * (_16 * a2 * a2 - _80 * a3 * a1 - _128 * a4 * a0))
            + dbg!(
                a3 * a3
                    * (a4 * a0 * (_144 * a2 * a0 - _6 * a1 * a1)
                        + (a0 * (_18 * a3 * a2 * a1 - _27 * a3 * a3 * a0 - _4 * a2 * a2 * a2)
                            + a1 * a1 * (a2 * a2 - _4 * a3 * a1)))
            )
            + dbg!(a4 * a1 * a1 * (_18 * a3 * a2 * a1 - _27 * a4 * a1 * a1 - _4 * a2 * a2 * a2));
        let pp = _8 * a4 * a2 - _3 * a3 * a3;
        let rr = a3 * a3 * a3 + _8 * a4 * a4 * a1 - _4 * a4 * a3 * a2;
        let delta0 = a2 * a2 - _3 * a3 * a1 + _12 * a4 * a0;
        let dd = _64 * a4 * a4 * a4 * a0 - _16 * a4 * a4 * a2 * a2 + _16 * a4 * a3 * a3 * a2
            - _16 * a4 * a4 * a3 * a1
            - _3 * a3 * a3 * a3 * a3;

        // Handle special cases
        let double_root = dbg!(discriminant) == F::zero();
        if dbg!(double_root) {
            let triple_root = double_root && dbg!(delta0) == F::zero();
            let quadruple_root = triple_root && dbg!(dd) == F::zero();
            let no_roots = dd == F::zero() && dbg!(pp) > F::zero() && dbg!(rr) == F::zero();
            if dbg!(quadruple_root) {
                // Wiki: all four roots are equal
                Roots::One([-a3 / (_4 * a4)])
            } else if dbg!(triple_root) {
                // Wiki: At least three roots are equal to each other
                // x0 is the unique root of the remainder of the Euclidean division of the quartic by its second derivative
                //
                // Solved by SymPy (ra is the desired reminder)
                // a, b, c, d, e = symbols('a,b,c,d,e')
                // f=a*x**4+b*x**3+c*x**2+d*x+e     // Quartic polynom
                // g=6*a*x**2+3*b*x+c               // Second derivative
                // q, r = div(f, g)                 // SymPy only finds the highest power
                // simplify(f-(q*g+r)) == 0         // Verify the first division
                // qa, ra = div(r/a,g/a)            // Workaround to get the second division
                // simplify(f-((q+qa)*g+ra*a)) == 0 // Verify the second division
                // solve(ra,x)
                // ----- yields
                // (−72*a^2*e+10*a*c^2−3*b^2*c)/(9*(8*a^2*d−4*a*b*c+b^3))
                let x0 = (-_72 * a4 * a4 * a0 + _10 * a4 * a2 * a2 - _3 * a3 * a3 * a2)
                    / (_9 * (_8 * a4 * a4 * a1 - _4 * a4 * a3 * a2 + a3 * a3 * a3));
                let roots = dbg!(Roots::One([x0]));
                roots.add_new_root(dbg!(-(a3 / a4 + _3 * x0)))
            } else if dbg!(no_roots) {
                // Wiki: two complex conjugate double roots
                Roots::No([])
            } else {
                find_roots_via_depressed_quartic(a4, a3, a2, a1, a0, pp, dbg!(rr), dd)
            }
        } else {
            let no_roots = dbg!(pp) > F::zero() || dbg!(dd) > F::zero();
            if dbg!(no_roots) {
                // Wiki: two pairs of non-real complex conjugate roots
                Roots::No([])
            } else {
                find_roots_via_depressed_quartic(a4, a3, a2, a1, a0, pp, dbg!(rr), dd)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::super::*;

    #[test]
    fn test_find_roots_quartic() {
        assert_eq!(find_roots_quartic(1f32, 0f32, 0f32, 0f32, 0f32), Roots::One([0f32]));
        assert_eq!(find_roots_quartic(1f64, 0f64, 0f64, 0f64, -1f64), Roots::Two([-1f64, 1f64]));
        assert_eq!(
            find_roots_quartic(1f64, -10f64, 35f64, -50f64, 24f64),
            Roots::Four([1f64, 2f64, 3f64, 4f64])
        );

        match find_roots_quartic(3f64, 5f64, -5f64, -5f64, 2f64) {
            Roots::Four(x) => {
                assert_float_array_eq!(2e-15f64, x, [-2f64, -1f64, 0.33333333333333333f64, 1f64]);
            }
            _ => {
                assert!(false);
            }
        }

        match find_roots_quartic(3f32, 5f32, -5f32, -5f32, 2f32) {
            Roots::Four(x) => {
                assert_float_array_eq!(5e-7, x, [-2f32, -1f32, 0.33333333333333333f32, 1f32]);
            }
            _ => {
                assert!(false);
            }
        }
    }

    #[test]
    fn test_find_roots_quartic_tim_luecke() {
        assert_eq!(
            find_roots_quartic(-14.0625f64, -3.75f64, 29.75f64, 4.0f64, -16.0f64),
            Roots::Two([-1.1016116464173349f64, 0.9682783130840016f64])
        );
        assert_eq!(
            find_roots_quartic(-14.0625f32, -3.75f32, 29.75f32, 4.0f32, -16.0f32),
            Roots::Two([-1.1016117f32, 0.96827835f32])
        );
    }

    #[test]
    fn test_find_roots_quartic_triple_root() {
        // (x+3)(3x-1)^3 == 27 x^4 + 54 x^3 - 72 x^2 + 26 x - 3
        assert_eq!(
            find_roots_quartic(27f64, 54f64, -72f64, 26f64, -3f64),
            Roots::Two([-3.0f64, 0.3333333333333333f64])
        );
        assert_eq!(
            find_roots_quartic(27f32, 54f32, -72f32, 26f32, -3f32),
            Roots::Two([-3.0f32, 0.33333333f32])
        );
    }

    #[test]
    fn test_find_roots_quartic_quadruple_root() {
        // (7x+2)^4 == 2401 x^4 + 2744 x^3 + 1176 x^2 + 224 x + 16
        assert_eq!(
            find_roots_quartic(2401f64, 2744f64, 1176f64, 224f64, 16f64),
            Roots::One([-0.2857142857142857f64])
        );
        // 32-bit floating point is less accurate
        assert_eq!(
            find_roots_quartic(2401f32, 2744f32, 1176f32, 224f32, 16f32),
            Roots::One([-0.2857143f32])
        );
    }
}
