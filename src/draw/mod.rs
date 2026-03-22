/*!
This module contains the [`geometry`] and [`intersection`] submodules.

- [`geometry`] holds the `draw` implementations for all
[`Primitive`](crate::primitive::Primitive) and
[`Composite`](crate::composite::Composite) types of this crate.
- [`drawable`] TODO
- [`intersection`] offers `draw` implementations for
[`Intersection`](crate::composite::Intersection)s.
 */

pub mod geometry;
pub use geometry::*;

pub mod intersection;
pub use intersection::*;

pub mod drawable;
pub use drawable::*;
