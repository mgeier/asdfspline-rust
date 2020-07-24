/*!
Rust implementation of the spline type that's used in the Audio Scene
Description Format (ASDF), see
<https://AudioSceneDescriptionFormat.readthedocs.io/>.

# Requirements

* Rust compiler, Cargo (<https://rustup.rs/>)

The required Rust packages (a.k.a. "crates") are listed in the file
`Cargo.toml`.

# Tests

```text
cargo test --all
```

There are further tests (using Python) in the `python/` directory.

# API Documentation

Run `cargo doc --all` in the main directory to create the documentation.
The generated HTML documentation can be accessed via
[target/doc/asdfspline/index.html](index.html) and
[target/doc/asdfspline_ffi/index.html](../asdfspline_ffi/index.html).

# Updating `README.md`

Using [cargo-readme](https://github.com/livioribeiro/cargo-readme) (`cargo install cargo-readme`):

```text
cargo readme -o README.md
```
*/
//use std::fmt::Debug;
//use std::ops::{Add, Div, DivAssign, Mul, Sub};
use std::ops::DivAssign;

//use num_traits::{Float, FromPrimitive, NumAssign};

pub mod asdfspline;
pub mod centripetalkochanekbartelsspline;
pub mod cubichermitespline;
pub mod monotonecubicspline;
pub mod piecewisecubiccurve;
pub mod shapepreservingcubicspline;
pub mod utilities;

pub use crate::asdfspline::AsdfPosSpline;
pub use crate::monotonecubicspline::MonotoneCubicSpline;
pub use crate::piecewisecubiccurve::PiecewiseCubicCurve;

/// A trait that is automatically implemented for all types that can be used as scalars,
/// e.g. time values.
pub trait Scalar: alga::general::RealField {}

impl<T: alga::general::RealField> Scalar for T {}

/// A trait that is automatically implemented for all types that can be used as positions,
/// polynomial coefficients, tangent vectors etc.
pub trait Vector<S: Scalar>: alga::linear::VectorSpace<Field = S> {}

impl<S: Scalar, T: alga::linear::VectorSpace<Field = S>> Vector<S> for T {}

pub trait VectorWithNorm<S: Scalar>: alga::linear::NormedSpace<RealField = S, ComplexField = S> + DivAssign<S> {}

impl<S: Scalar, T: alga::linear::NormedSpace<RealField = S, ComplexField = S> + DivAssign<S>> VectorWithNorm<S> for T {}
