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

extern crate test;
use self::test::Bencher;

use super::*;

fn x2_min_1(x:f64) -> f64 { x*x-1f64 }
fn x2_min_1_d(x:f64) -> f64 { 2f64*x }

#[bench]
fn bench_10_numerical_newton_raphson_x2_min_1(b: &mut Bencher) {
  let conv = SimpleConvergency {eps:1e-15, max_iter:30};

  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_root_newton_raphson(10f64, &x2_min_1, &x2_min_1_d, &conv).ok().unwrap();
    }
  } );
}
#[bench]
fn bench_10_numerical_regula_falsi_x2_min_1(b: &mut Bencher) {
  let conv = SimpleConvergency {eps:1e-15, max_iter:30};

  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_root_regula_falsi(10f64, 0f64, &x2_min_1, &conv).ok().unwrap();
    }
  } );
}

#[bench]
fn bench_10_numerical_secant_x2_min_1(b: &mut Bencher) {
  let conv = SimpleConvergency {eps:1e-15, max_iter:30};

  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_root_secant(0f64, 10f64, &x2_min_1, &conv).ok().unwrap();
    }
  } );
}

#[bench]
fn bench_10_numerical_brent_x2_min_1(b: &mut Bencher) {
  let conv = SimpleConvergency {eps:1e-15, max_iter:30};

  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_root_brent(0f64, 10f64, &x2_min_1, &conv).ok().unwrap();
    }
  } );
}

#[bench]
fn bench_10_linear_x_plus_1(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_linear(1f64, 1f64);
    }
  } );
}

#[bench]
fn bench_10_quadratic_x2_plus_x_min_1(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_quadratic(1f64, 1f64, -1f64);
    }
  } );
}

#[bench]
fn bench_10_cubic_depressed_x3_plus_x_min_1(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_cubic_depressed(1f64, -1f64);
    }
  } );
}

#[bench]
fn bench_10_cubic_normalized_x3_min_x2_plus_x_min_1(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_cubic_normalized(-1f64, 1f64, -1f64);
    }
  } );
}

#[bench]
fn bench_10_cubic_2x3_min_x2_plus_x_min_1(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_cubic(2f64, -1f64, 1f64, -1f64);
    }
  } );
}

#[bench]
fn bench_10_quartic_3x4_plus_5x3_min_5x2_min_5x_plus_2(b: &mut Bencher) {
  b.iter( || {
    for _x in 0..test::black_box(10) {
      let _y = find_roots_quartic(3f64, 5f64, -5f64, -5f64, 2f64);
    }
  } );
}
