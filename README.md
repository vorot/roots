#roots [![Build Status](https://travis-ci.org/vorot/roots.svg)](https://travis-ci.org/vorot/roots)[![Crates.io](https://img.shields.io/crates/v/roots.svg)](https://crates.io/crates/roots) #

Library of well known algorithms for numerical root finding.

## Features
  - Iterative approximation:
   - [Newton-Raphson](http://en.wikipedia.org/wiki/Newton%27s_method) method
   - [Secant](http://en.wikipedia.org/wiki/Secant_method) method
   - [Regula falsi](http://en.wikipedia.org/wiki/False_position_method) method (with Illinois modification)
   - [Brent-Dekker](http://en.wikipedia.org/wiki/Brent%27s_method) method
  - Solving polynomial equations
   - [Linear](http://en.wikipedia.org/wiki/Linear_equation) equation (editors' choice)
   - [Quadratic](http://en.wikipedia.org/wiki/Quadratic_equation) equation
   - [Cubic](http://en.wikipedia.org/wiki/Cubic_function) equation
   - [Quartic](http://en.wikipedia.org/wiki/Quartic_function) equation

## Usage

```rust
extern crate roots;
use roots::Roots;
use roots::find_roots_cubic;
use roots::find_root_brent;
use roots::find_root_secant;
use roots::SimpleConvergency;

// Find the root of a complex function in the area determined by a simpler polynom
fn find_solution(enormous_function: &Fn(f64)->f64, root_area_polynom:(f64,f64,f64,f64)) -> Option<f64> {
  // de-structure polynom coefficients
  match root_area_polynom {
    (a3,a2,a1,a0) => {
      // Set convergency conditions: required precision - 1e-8, max 30 iterations
      let conv = SimpleConvergency{eps:1e-8f64; max_iter:30};
      // Find root area by solving the polynom
      match find_roots_cubic(a3,a2,a1,a0) {
        // Try to find the root by one of iterative methods
        Roots::Three(roots) => { find_root_brent(roots[0],roots[2],enormous_function,conv).ok() },
        Roots::Two(roots) => { find_root_brent(roots[0],roots[1],enormous_function,conv).ok() },
        Roots::One(roots) => { find_root_secant(roots[0]-1f64,roots[0]+1f64,enormous_function,conv).ok() },
        _ => None,
      }
    },
    _ => None,
  }
}
```
