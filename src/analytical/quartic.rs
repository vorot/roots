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
        let discriminant = 
            _256 * a4 * a4 * a4 * a0 * a0 * a0
          - _192 * a4 * a4 * a3 * a1 * a0 * a0
          - _128 * a4 * a4 * a2 * a2 * a0 * a0
          + _144 * a4 * a4 * a2 * a1 * a1 * a0
          -  _27 * a4 * a4 * a1 * a1 * a1 * a1
          + _144 * a4 * a3 * a3 * a2 * a0 * a0
          -   _6 * a4 * a3 * a3 * a1 * a1 * a0
          -  _80 * a4 * a3 * a2 * a2 * a1 * a0
          +  _18 * a4 * a3 * a2 * a1 * a1 * a1
          +  _16 * a4 * a2 * a2 * a2 * a2 * a0
          -   _4 * a4 * a2 * a2 * a2 * a1 * a1
          -  _27 * a3 * a3 * a3 * a3 * a0 * a0
          +  _18 * a3 * a3 * a3 * a2 * a1 * a0
          -   _4 * a3 * a3 * a3 * a1 * a1 * a1
          -   _4 * a3 * a3 * a2 * a2 * a2 * a0
          +        a3 * a3 * a2 * a2 * a1 * a1
          ;
        let pp = _8 * a4 * a2 - _3 * a3 * a3;
        let rr = a3 * a3 * a3 + _8 * a4 * a4 * a1 - _4 * a4 * a3 * a2;
        let delta0 = a2 * a2 - _3 * a3 * a1 + _12 * a4 * a0;
        let dd = 
             _64 * a4 * a4 * a4 * a0
           - _16 * a4 * a4 * a2 * a2
           + _16 * a4 * a3 * a3 * a2
           - _16 * a4 * a4 * a3 * a1
           -  _3 * a3 * a3 * a3 * a3;
        println!("quartic ∆:{:?},P:{:?},R:{:?},∆0:{:?},D:{:?}",discriminant,pp,rr,delta0,dd);

        // Handle special cases
        if discriminant == F::zero() && delta0 == F::zero() && dd == F::zero() {
            // Wiki: all four roots are equal
            Roots::One([-a3/(_4*a4)])
        }
        else if discriminant == F::zero() && delta0 == F::zero() {
            println!("quartic discriminant is zero, multiple roots path");
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
            let x0 = (_72 * a4 * a4 * a0 + _10 * a4 * a2 * a2 - _3 * a1 * a1 * a2) / (_9 * (_8*a4*a4*a1-_4*a4*a3*a2+a3*a3*a3));
            let roots = Roots::One([x0]);
            roots.add_new_root(-(a3/a4+_3*x0))
        }
        else if discriminant == F::zero() && dd == F::zero() && pp > F::zero() && rr == F::zero() {
            // Wiki: two complex conjugate double roots
            println!("no real roots path");
            Roots::No([])
        }
        else {
            // Depressed quartic
            // https://en.wikipedia.org/wiki/Quartic_function#Converting_to_a_depressed_quartic

            // a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0 => x^4 + a*x^3 + b*x^2 + c*x + d = 0.
            let (a, b, c, d) = (a3 / a4, a2 / a4, a1 / a4, a0 / a4);
            // x^4 + a*x^3 + b*x^2 + c*x + d = 0 => y^4 + p*y^2 + q*y + r.
            let a_pow_2 = a * a;
            let a_pow_3 = a_pow_2 * a;
            let a_pow_4 = a_pow_2 * a_pow_2;
            let subst = -a3 / (F::four() * a4);
            let p = b - F::three() * a_pow_2 / _8;
            let q = a_pow_3 / _8 - a * b / F::two() + c;
            let r = d - F::three() * a_pow_4 / _256 - c * a / F::four() + a_pow_2 * b / _16;

            println!("Quartic depressed p:{:?},q:{:?},r:{:?}",p,q,r);
            let mut roots = Roots::No([]);
            for x in super::quartic_depressed::find_roots_quartic_depressed(p, q, r)
                .as_ref()
                .iter()
            {
                println!("Quartic depressed root {:?}",*x);
                roots = roots.add_new_root(*x + subst);
            }
            roots
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
        assert_eq!(find_roots_quartic(-14.0625f64, -3.75f64, 29.75f64, 4.0f64, -16.0f64), Roots::Two([-1.1016117f64,0.96827835f64]));
        assert_eq!(find_roots_quartic(-14.0625f32, -3.75f32, 29.75f32, 4.0f32, -16.0f32), Roots::Two([-1.1016117f32,0.96827835f32]));
    }
}
