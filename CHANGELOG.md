# Change Log
All notable changes to this project will be documented in this file.
This project is going to adhere to [Semantic Versioning](http://semver.org/)
after version 1.0.0.

## [Unreleased]

## [0.0.8] - 2022-12-21
* Remove debug print - thanks to J-F-Liu, Ralith
* Fixed discriminant check in quartic equations - thanks to jingnanshi
* Allowed mutable closures for secant etc. - thanks to vron

## [0.0.7] - 2021-06-17
* SearchError becomes public - thanks to JP-Ellis
* Unnormalized cubic equations are solved using the general formula rather than trigonometrically - thanks to Logicalshift

## [0.0.6] - 2019-12-22
* Fixed cubic equations with very small a3 - thanks to Andrew Hunter
* Improved quartic equations with multiple roots (for f64; f32 is still a problem) - thanks to Tim Lueke
* Removed warnings of rustc 1.40.0
* Switched benchmarks from Bencher to Criterion

## [0.0.5] - 2019-01-20
* Trait Error implemented for SearchError - thanks to phillyfan1138
* Find roots of higher-degree polynomials using eigenvalues - thanks to stiv-yakovenko
* Find roots of higher-degree polynomials using Sturm's theorem recursively (experimental)
* Inverse quadratic approximation

## [0.0.4] - 2017-09-05
* Reduced the performance overhead by using generics - thanks to aepsil0n
* Handle special cases of quadratic equations - thanks to stiv-yakovenko

## [0.0.3] - 2017-03-28
* New version of the compiler
* Benchmarks
* Improved the speed of the Brent-Dekker method
* Reduced the convergency boilerplate

## [0.0.2] - 2015-06-08
* Fight against the compiler

## [0.0.1] - 2015-03-24
* Initial version

[Unreleased]: https://github.com/vorot/roots/compare/v0.0.8...HEAD
[0.0.8]: https://github.com/vorot/roots/compare/v0.0.7...v0.0.8
[0.0.7]: https://github.com/vorot/roots/compare/v0.0.6...v0.0.7
[0.0.6]: https://github.com/vorot/roots/compare/v0.0.5...v0.0.6
[0.0.5]: https://github.com/vorot/roots/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/vorot/roots/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/vorot/roots/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/vorot/roots/compare/v0.0.1...v0.0.2
