/*
eigen.rs 0.2

This piece of code is transpiled EigenvalueDecomposition.java from Jama framework.
I had to do this, because I haven't found any buildable rust opensource project
to calculate real eigen values in Rust. There are many packages which are usually bound to
openBLAS, but I can't build them with gnu toolchain, and that's the only toolchain, which
allows gdb (and thus IDE) debug nowadays.

Quality code is far from perfect, hopefully someone will appreciate my one day of
manual code conversion nightmare and mention me in the source code.

Stepan Yakovenko,
https://github.com/stiv-yakovenko

*/
use std::fmt;
use std::ops::Index;
use std::collections::VecDeque;
use std::ops::IndexMut;
use std::cmp;
extern crate num_traits;

use eigen::num_traits::Float;
//use self::num_traits;

pub struct Matrix {
    data:VecDeque<f64>,
    n:usize
}
impl Matrix {
    pub fn new(n:usize)->Matrix{
        let mut data = VecDeque::new();
        data.resize(n * n, 0.);
        Matrix{data:data,n:n}
    }
}
impl fmt::Debug for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{").ok();
        for r in 0..self.n {
            for c in 0..self.n {
                write!(f, "{:.3?} ", self[[c,r]]).ok();
            }
            writeln!(f, "").ok();
        }
        write!(f, "}}")
    }
}
impl Index<[usize; 2]> for Matrix{
    type Output = f64;
    fn index(&self, idx: [usize; 2]) -> &f64 {
        &self.data[idx[0]+self.n*idx[1]]
    }
}
impl IndexMut<[usize; 2]> for Matrix{
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut f64 {
        self.data.get_mut(idx[0]+self.n*idx[1]).unwrap()
    }
}
#[allow(dead_code)]
fn cdiv(xr: f64, xi: f64, yr: f64, yi: f64) -> (f64, f64) {
    let r: f64;
    let d: f64;
    if yr.abs() > yi.abs() {
        r = yi / yr;
        d = yr + r * yi;
        ((xr + r * xi) / d, (xi - r * xr) / d)
    } else {
        r = yr / yi;
        d = yi + r * yr;
        ((r * xr + xi) / d, (r * xi - xr) / d)
    }    
}

#[allow(dead_code)]
pub fn hqr2(n_in: usize, h: &mut Matrix, v: &mut Matrix, d: &mut Vec<f64>, e: &mut Vec<f64>) {
    //  This is derived from the Algol procedure hqr2,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    // Initialize
    let nn = n_in;
    let mut n = nn as i16 - 1;
    let low = 0;
    let high = nn - 1;
    let eps = (2.0).powf(-52.0);
    let mut exshift = 0.0;
    let mut p = 0.;
    let mut q = 0.;
    let mut r = 0.;
    let mut s = 0.;
    let mut z = 0.;
    let mut t;
    let mut w;
    let mut x;
    let mut y;
    // Store roots isolated by balanc and compute matrix norm
    let mut norm = 0.0;
    let mut i = 0 as usize;
    while i < nn {
        if i < low || i > high {
            d[i] = h[[i, i]];
            e[i] = 0.0;
        }
        let mut j = cmp::max(i as i16 - 1, 0) as usize;
        while j < nn {
            norm = norm + (h[[i, j]]).abs();
            j = j + 1;
        }
        i = i + 1;
    }
    // Outer loop over eigenvalue index
    let mut iter = 0;
    while n >= low as i16 {
        // Look for single small sub-diagonal element
        let mut l = n;
        while l > low as i16 {
            s = (h[[l as usize - 1, l as usize - 1]]).abs() + (h[[l as usize, l as usize]]).abs();
            if s == 0.0 {
                s = norm;
            }
            if (h[[l as usize, l as usize - 1]]).abs() < eps * s {
                break;
            }
            l = l - 1;
        }
        // Check for convergence
        // One root found
        if l == n {
            h[[n as usize, n as usize]] = h[[n as usize, n as usize]] + exshift;
            d[n as usize] = h[[n as usize, n as usize]];
            e[n as usize] = 0.0;
            n = n - 1;
            iter = 0;
            // Two roots found
        } else if l == n - 1 {
            w = h[[n as usize, n as usize - 1]] * h[[n as usize - 1, n as usize]];
            p = (h[[n as usize - 1, n as usize - 1]] - h[[n as usize, n as usize]]) / 2.0;
            q = p * p + w;
            z = (q).abs().sqrt();
            h[[n as usize, n as usize]] = h[[n as usize, n as usize]] + exshift;
            h[[n as usize - 1, n as usize - 1]] = h[[n as usize - 1, n as usize - 1]] + exshift;
            x = h[[n as usize, n as usize]];
            // Real pair
            if q >= 0. {
                if p >= 0. {
                    z = p + z;
                } else {
                    z = p - z;
                }
                d[n as usize - 1] = x + z;
                d[n as usize] = d[n as usize - 1];
                if z != 0.0 {
                    d[n as usize] = x - w / z;
                }
                e[n as usize - 1] = 0.0;
                e[n as usize] = 0.0;
                x = h[[n as usize, n as usize - 1]];
                s = (x).abs() + (z).abs();
                p = x / s;
                q = z / s;
                r = (p * p + q * q).sqrt();
                p = p / r;
                q = q / r;
                // Row modification
                let mut j = n - 1;
                while j < nn as i16 {
                    z = h[[n as usize - 1, j as usize]];
                    h[[n as usize - 1, j as usize]] = q * z + p * h[[n as usize, j as usize]];
                    h[[n as usize, j as usize]] = q * h[[n as usize, j as usize]] - p * z;
                    j = j + 1;
                }
                // Column modification
                let mut i = 0;
                while i <= n {
                    z = h[[i as usize, n as usize - 1]];
                    h[[i as usize, n as usize - 1]] = q * z + p * h[[i as usize, n as usize]];
                    h[[i as usize, n as usize]] = q * h[[i as usize, n as usize]] - p * z;
                    i = i + 1;
                }
                // Accumulate transformations
                let mut i = low;
                while i <= high {
                    z = v[[i as usize, n as usize - 1]];
                    v[[i as usize, n as usize - 1]] = q * z + p * v[[i as usize, n as usize]];
                    v[[i as usize, n as usize]] = q * v[[i as usize, n as usize]] - p * z;
                    i = i + 1;
                }
                // Complex pair
            } else {
                d[n as usize - 1] = x + p;
                d[n as usize] = x + p;
                e[n as usize - 1] = z;
                e[n as usize] = -z;
            }
            n = n - 2;
            iter = 0;
            // No convergence yet
        } else {
            // Form shift
            x = h[[n as usize, n as usize]];
            y = 0.0;
            w = 0.0;
            if l < n {
                y = h[[n as usize - 1, n as usize - 1]];
                w = h[[n as usize, n as usize - 1]] * h[[n as usize - 1, n as usize]];
            }
            // Wilkinson's original ad hoc shift
            if iter == 10 {
                exshift += x;
                let mut i = low;
                while i <= n as usize {
                    h[[i, i]] -= x;
                    i = i + 1;
                }
                s = (h[[n as usize, n as usize - 1]]).abs() + (h[[n as usize - 1, n as usize - 2]]).abs();
                y = 0.75 * s;
                x = y;
                w = -0.4375 * s * s;
            }
            // MATLAB's new ad hoc shift
            if iter == 30 {
                s = (y - x) / 2.0;
                s = s * s + w;
                if s > 0. {
                    s = s.sqrt();
                    if y < x {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    let mut i = low;
                    while i <= n as usize {
                        h[[i, i]] -= s;
                        i = i + 1;
                    }
                    exshift += s;
                    x = 0.964;
                    y = x;
                    w = y;
                }
            }
            iter = iter + 1;   // (Could check iteration count here.)
            // Look for two consecutive small sub-diagonal elements
            let mut m = n - 2;
            while m >= l {
                z = h[[m as usize, m as usize]];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[[m as usize + 1, m as usize]] + h[[m as usize, m as usize + 1]];
                q = h[[m as usize + 1, m as usize + 1]] - z - r - s;
                r = h[[m as usize + 2, m as usize + 1]];
                s = (p).abs() + (q).abs() + (r).abs();
                p = p / s;
                q = q / s;
                r = r / s;
                if m == l {
                    break;
                }
                if h[[m as usize, m as usize - 1]].abs() * (q).abs() + (r).abs() <
                    eps * ((p).abs() * ((h[[m as usize - 1, m as usize - 1]]).abs() + (z).abs() +
                        (h[[m as usize + 1, m as usize + 1]]).abs()
                    )) {
                    break;
                }
                m = m - 1;
            }
            let mut i = m + 2;
            while i <= n {
                h[[i as usize, i as usize - 2]] = 0.0;
                if i > m + 2 {
                    h[[i as usize, i as usize - 3]] = 0.0;
                }
                i = i + 1;
            }
            // Double QR step involving rows l:n and columns m:n
            let mut k = m;
            while k <= n - 1 {
                let notlast = if k != n - 1 { true } else { false };
                if k != m {
                    p = h[[k as usize, k as usize - 1]];
                    q = h[[k as usize + 1, k as usize - 1]];
                    r = if notlast { h[[k as usize + 2, k as usize - 1]] } else { 0.0 };
                    x = (p).abs() + (q).abs() + (r).abs();
                    if x == 0.0 {
                        k=k+1;
                        continue;
                    }
                    p = p / x;
                    q = q / x;
                    r = r / x;
                }
                s = (p * p + q * q + r * r).sqrt();
                if p < 0. {
                    s = -s;
                }
                if s != 0. {
                    if k != m {
                        h[[k as usize, k as usize - 1]] = -s * x;
                    } else if l != m {
                        h[[k as usize, k as usize - 1]] = -h[[k as usize, k as usize - 1]];
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;
                    // Row modification
                    let mut j = k;
                    while j < nn as i16 {
                        p = h[[k as usize, j as usize]] + q * h[[k as usize + 1, j as usize]];
                        if notlast {
                            p = p + r * h[[k as usize + 2, j as usize]];
                            h[[k as usize + 2, j as usize]] = h[[k as usize + 2, j as usize]] - p * z;
                        }
                        h[[k as usize, j as usize]] = h[[k as usize, j as usize]] - p * x;
                        h[[k as usize + 1, j as usize]] = h[[k as usize + 1, j as usize]] - p * y;
                        j = j + 1;
                    }
                    // Column modification
                    let mut i = 0;
                    while i <= cmp::min(n as usize, k as usize + 3) {
                        p = x * h[[i, k as usize]] + y * h[[i as usize, k as usize + 1]];
                        if notlast {
                            p = p + z * h[[i, k as usize + 2]];
                            h[[i, k as usize + 2]] = h[[i, k as usize + 2]] - p * r;
                        }
                        h[[i, k as usize]] = h[[i, k as usize]] - p;
                        h[[i, k as usize + 1]] = h[[i, k as usize + 1]] - p * q;
                        i = i + 1;
                    }
                    // Accumulate transformations
                    let mut i = low;
                    while i <= high {
                        p = x * v[[i, k as usize]] + y * v[[i, k as usize + 1]];
                        if notlast {
                            p = p + z * v[[i as usize, k as usize + 2]];
                            v[[i as usize, k as usize + 2]] = v[[i as usize, k as usize + 2]] - p * r;
                        }
                        v[[i, k as usize]] = v[[i, k as usize]] - p;
                        v[[i, k as usize + 1]] = v[[i, k as usize + 1]] - p * q;
                        i = i + 1;
                    }
                }  // (s != 0)
                k = k + 1;
            }  // k loop
        }  // check convergence
    }  // while n >= low
    // Backsubstitute to find vectors of upper triangular form
    if norm == 0.0 {
        return;
    }
    n = nn as i16 - 1;
    while n >= 0 {
        p = d[n as usize];
        q = e[n as usize];
        // Real vector
        if q == 0. {
            let mut l = n;
            h[[n as usize, n as usize]] = 1.0;
            let mut i = n as i16 - 1;
            while i >= 0 {
                w = h[[i as usize, i as usize]] - p;
                r = 0.0;
                let mut j = l;
                while j <= n {
                    r = r + h[[i as usize, j as usize]] * h[[j as usize, n as usize]];
                    j = j + 1;
                }
                if e[i as usize] < 0.0 {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if e[i as usize] == 0.0 {
                        if w != 0.0 {
                            h[[i as usize, n as usize]] = -r / w;
                        } else {
                            h[[i as usize, n as usize]] = -r / (eps * norm);
                        }
                        // Solve real equations
                    } else {
                        x = h[[i as usize, i as usize + 1]];
                        y = h[[i as usize + 1, i as usize]];
                        q = (d[i as usize] - p) * (d[i as usize] - p) + e[i as usize] * e[i as usize];
                        t = (x * s - z * r) / q;
                        h[[i as usize, n as usize]] = t;
                        if (x).abs() > (z).abs() {
                            h[[i as usize + 1, n as usize]] = (-r - w * t) / x;
                        } else {
                            h[[i as usize + 1, n as usize]] = (-s - y * t) / z;
                        }
                    }
                    // Overflow control
                    t = h[[i as usize, n as usize]];
                    if (eps * t).abs() * t > 1. {
                        let mut j = i;
                        while j <= n as i16 {
                            h[[j as usize, n as usize]] = h[[j as usize, n as usize]] / t;
                            j = j + 1;
                        }
                    }
                }
                i = i - 1;
            }
            // Complex vector
        } else if q < 0. {
            let mut l = n - 1;
            // Last vector component imaginary so matrix is triangular
            if (h[[n as usize, n as usize - 1]]).abs() > (h[[n as usize - 1, n as usize]]).abs() {
                h[[n as usize - 1, n as usize - 1]] = q / h[[n as usize, n as usize - 1]];
                h[[n as usize - 1, n as usize]] = -(h[[n as usize, n as usize]] - p) / h[[n as usize, n as usize - 1]];
            } else {
                let (cdivr, cdivi) = cdiv(0.0, -h[[n as usize - 1, n as usize]], h[[n as usize - 1, n as usize - 1]] - p, q);
                h[[n as usize - 1, n as usize - 1]] = cdivr;
                h[[n as usize - 1, n as usize]] = cdivi;
            }
            h[[n as usize, n as usize - 1]] = 0.0;
            h[[n as usize, n as usize]] = 1.0;
            let mut i = n - 2;
            while i >= 0 {
                let mut ra = 0.;
                let mut sa = 0.;
                let mut vr;
                let vi;
                let mut j = l;
                while j <= n {
                    ra = ra + h[[i as usize, j as usize]] * h[[j as usize, n as usize - 1]];
                    sa = sa + h[[i as usize, j as usize]] * h[[j as usize, n as usize]];
                    j = j + 1;
                }
                w = h[[i as usize, i as usize]] - p;
                if e[i as usize] < 0.0 {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if e[i as usize] == 0. {
                        let (cdivr, cdivi) = cdiv(-ra, -sa, w, q);
                        h[[i as usize, n as usize - 1]] = cdivr;
                        h[[i as usize, n as usize]] = cdivi;
                    } else {
                        // Solve complex equations
                        x = h[[i as usize, i as usize + 1]];
                        y = h[[i as usize + 1, i as usize]];
                        vr = (d[i as usize] - p) * (d[i as usize] - p) + e[i as usize] * e[i as usize] - q * q;
                        vi = (d[i as usize] - p) * 2.0 * q;
                        if vr == 0.0 && vi == 0.0 {
                            vr = eps * norm * ((w).abs() + (q).abs() +
                                (x).abs() + (y).abs() + (z)).abs();
                        }
                        let (cdivr, cdivi) = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                        h[[i as usize, n as usize - 1]] = cdivr;
                        h[[i as usize, n as usize]] = cdivi;
                        if (x).abs() > ((z).abs() + (q).abs()) {
                            h[[i as usize + 1, n as usize - 1]] = (-ra - w * h[[i as usize, n as usize - 1]] + q * h[[i as usize, n as usize]]) / x;
                            h[[i as usize + 1, n as usize]] = (-sa - w * h[[i as usize, n as usize]] - q * h[[i as usize, n as usize - 1]]) / x;
                        } else {
                            let (cdivr, cdivi) = cdiv(-r - y * h[[i as usize, n as usize - 1]], -s - y * h[[i as usize, n as usize]], z, q);
                            h[[i as usize + 1, n as usize - 1]] = cdivr;
                            h[[i as usize + 1, n as usize]] = cdivi;
                        }
                    }
                    // Overflow control
                    t = (h[[i as usize, n as usize - 1]]).abs().max(h[[i as usize, n as usize]].abs());
                    if (eps * t) * t > 1. {
                        let mut j = i;
                        while j <= n {
                            j = j + 1;
                            h[[j as usize, n as usize - 1]] = h[[j as usize, n as usize - 1]] / t;
                            h[[j as usize, n as usize]] = h[[j as usize, n as usize]] / t;
                        }
                    }
                }
                i = i - 1;
            }
        }
        n = n - 1;
    }
    // Vectors of isolated roots
    let mut i = 0;
    while i < nn {
        if i < low || i > high {
            let mut j = i;
            while j < nn {
                v[[i, j]] = h[[i, j]];
                j = j + 1;
            }
        }
        i = i + 1;
    }
    // Back transformation to get eigenvectors of original matrix
    let mut j = nn as i16 - 1;
    while j >= low as i16 {
        let mut i = low;
        while i <= high {
            z = 0.0;
            let mut k = low;
            while k <= cmp::min(j as usize, high) {
                z = z + v[[i, k]] * h[[k, j as usize]];
                k = k + 1;
            }
            v[[i, j as usize]] = z;
            i = i + 1;
        }
        j = j - 1;
    }
}

//  This is derived from the Algol procedures orthes and ortran,
//  by Martin and Wilkinson, Handbook for Auto. Comp.,
//  Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutines in EISPACK.
#[allow(dead_code)]
pub fn orthes(m: &mut Matrix, h_mat: &mut Matrix, v_mat: &mut Matrix) {
    let low = 0;
    let n = m.n;
    let high = n - 1;
    let mut m = low + 1;
    let mut ort = vec!(0.; n);
       while m < high - 1 {
        // Scale column.
        let mut scale = 0.0;
        let mut i = m;
        //for (int        i = m;        i < = high;        i + +)
        while i <= high {
            scale = scale + (h_mat[[i, m - 1]]).abs();
            i = i + 1;
        }
        if scale != 0.0 {
            // Compute Householder transformation.
            let mut h = 0.0;
            let mut i = high;
            while i >= m {
                ort[i] = h_mat[[i, m - 1]] / scale;
                h += ort[i] * ort[i];
                i = i - 1;
            }
            let mut g = h.sqrt();
            if ort[m] > 0. {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)
            let mut j = m;
            while j < n {
                let mut f = 0.0;
                let mut i = high;
                while i >= m {
                    f += ort[i] * h_mat[[i, j]];
                    i = i - 1;
                }
                f = f / h;
                let mut i = m;
                while
                    i <= high {
                    h_mat[[i, j]] -= f * ort[i];
                    i = i + 1;
                }
                j = j + 1;
            }
            let mut i = 0;
            while i <= high {
                let mut f = 0.0;
                let mut j = high;
                while j >= m {
                    f += ort[j] * h_mat[[i, j]];
                    j = j - 1;
                }
                f = f / h;
                let mut j = m;
                while j <= high {
                    h_mat[[i, j]] -= f * ort[j];
                    j = j + 1;
                }
                i = i + 1;
            }
            ort[m] = scale * ort[m];
            h_mat[[m, m - 1]] = scale * g;
        }
        m = m + 1;
    }
    // Accumulate transformations (Algol's ortran).
    for i in 0..n {
        for j in 0..n {
            v_mat[[i, j]] = if i == j { 1.0 } else { 0.0 };
        }
    }
    let mut m = high - 1;
    while m >= low + 1 {
        if h_mat[[m, m - 1]] != 0.0 {
            let mut i = m + 1;
            while i <= high {
                ort[i] = h_mat[[i, m - 1]];
                i = i + 1;
            }
            let mut j = m;
            while j <= high {
                let mut g = 0.0;
                let mut i = m;
                while i <= high {
                    g += ort[i] * v_mat[[i, j]];
                    i = i + 1;
                }
                // Double division avoids possible underflow
                g = (g / ort[m]) / h_mat[[m, m - 1]];
                let mut i = m;
                while i <= high {
                    v_mat[[i, j]] += g * ort[i];
                    i = i + 1;
                }
                j = j + 1;
            }
        }
        m = m - 1;
    }
}

#[allow(dead_code)]
fn calc_eigen(m: &mut Matrix) -> Vec<(f64, f64)>{
    let n = m.n;
    let mut h_mat = Matrix::new(n);
    let mut v_mat = Matrix::new(n);
    //let mut ort = vec!(0.;n);
    let mut d = vec!(0.; n);
    let mut e = vec!(0.; n);
    for i in 0..n {
        for j in 0..n {
            h_mat[[i, j]] = m[[i, j]];
        }
    }
    orthes(m, &mut h_mat, &mut v_mat);
    hqr2(n, &mut h_mat, &mut v_mat, &mut d, &mut e);
    let mut r= vec!((0.,0.); n);
    for i in 0..n {
        r[i]=(d[i],e[i])
    }
    r
}

#[allow(dead_code)]
// 6 x x x  + 9 x x + 3  x + 1 == 0 is vec![1.,3.,9.,6.] 
pub fn solve_poly(c: Vec<f64>) -> VecDeque<f64> {
    let n = c.len();
    let mut m = Matrix::new(n);
    for i in 0..(n-1) {
        m[[i+1,i]]=1.;
    }
    for i in 0..(n) {
        m[[i,n-1]]=-c[i];
    }
    //println!("{:?}",m);
    let ei = calc_eigen(&mut m);
    let mut r = VecDeque::new();
    for c in ei {
        if c.1*c.1==0. {
            r.push_back(c.0);
        }
    }
    r
}

/*fn main() {
    let c =vec![1.,3.,9.,6.];
    let r =solve_poly(c);
    println!("{:?}",r);
    //    let m = Matrix::new(10);
    //    *m.get(5, 5) = 7.;
}
*/