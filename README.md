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
use roots::find_roots_cubic;
use roots::find_root_brent;
use roots::find_root_secant;
use roots::SimpleConvergency;

// Find the root of a complex function in the area determined by a simpler polynom
fn find_solution(enormous_function: &Fn(f64)->f64, root_area_polynom:(f64,f64,f64,f64)) -> Option<f64> {
  // destructure polynom coefficients
  match root_area_polynom {
    (a3,a2,a1,a0) => {
      // Set convergency conditions: required precision - 1e-8, max 30 iterations
      let conv = SimpleConvergency{eps:1e-8f64; max_iter:30};
      // Find root area by solving the polynom
      match find_roots_cubic(a3,a2,a1,a0).as_slice() {
        // Try to find the root by one of iterative methods
        [x1,x2,x3] => { find_root_brent(x1,x3,enormous_function,conv).ok() },
        [x1,x2] => { find_root_brent(x1,x2,enormous_function,conv).ok() },
        [x1] => { find_root_secant(x1-1f64,x1+1f64,enormous_function,conv).ok() },
        _ => None,
      }
    },
    _ => None,
  }
}
```
