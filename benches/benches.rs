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

#[macro_use]
extern crate bencher;
extern crate roots;
use self::bencher::Bencher;
use roots::find_root_brent;
use roots::find_root_newton_raphson;
use roots::find_root_regula_falsi;
use roots::find_root_secant;
use roots::find_roots_biquadratic;
use roots::find_roots_quadratic;
use roots::find_roots_quartic;

fn x2_min_1(x: f64) -> f64 {
    x * x - 1f64
}

fn x4_min_1(x: f64) -> f64 {
    x * x * x * x - 1f64
}

fn x2_min_1_derivative(x: f64) -> f64 {
    2f64 * x
}

fn x4_min_1_derivative(x: f64) -> f64 {
    4f64 * x * x * x
}

fn secant_x2_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_secant(0f64, 10f64, &x2_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn secant_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_secant(0f64, 10f64, &x4_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn regula_falsi_x2_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_regula_falsi(0f64, 10f64, &x2_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn regula_falsi_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_regula_falsi(0f64, 10f64, &x4_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn brent_x2_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_brent(0f64, 10f64, &x2_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn brent_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_brent(0f64, 10f64, &x4_min_1, &mut 1e-15f64).ok().unwrap();
        }
    });
}

fn newton_raphson_x2_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_newton_raphson(0.5f64, &x2_min_1, &x2_min_1_derivative, &mut 1e-15f64)
                .ok()
                .unwrap();
        }
    });
}

fn newton_raphson_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_root_newton_raphson(0.5f64, &x4_min_1, &x4_min_1_derivative, &mut 1e-15f64)
                .ok()
                .unwrap();
        }
    });
}

fn quadratic_x2_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_roots_quadratic(1f64, 0f64, -1f64);
        }
    });
}

fn biquadratic_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_roots_biquadratic(1f64, 0f64, -1f64);
        }
    });
}

fn quartic_x4_min_1_x1000(b: &mut Bencher) {
    b.iter(|| {
        for _x in 0..1000 {
            let _y = find_roots_quartic(1f64, 0f64, 0f64, 0f64, -1f64);
        }
    });
}

benchmark_group!(
    benches,
    quadratic_x2_min_1_x1000,
    biquadratic_x4_min_1_x1000,
    quartic_x4_min_1_x1000,
    secant_x2_min_1_x1000,
    secant_x4_min_1_x1000,
    regula_falsi_x2_min_1_x1000,
    regula_falsi_x4_min_1_x1000,
    brent_x2_min_1_x1000,
    brent_x4_min_1_x1000,
    newton_raphson_x2_min_1_x1000,
    newton_raphson_x4_min_1_x1000
);
benchmark_main!(benches);
