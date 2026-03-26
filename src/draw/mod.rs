/*!
This module contains the [`geometry`] and [`intersection`] submodules.

- [`geometry`] holds the `draw` implementations for all
[`Primitive`](crate::primitive::Primitive) and
[`Composite`](crate::composite::Composite) types of this crate.
- [`drawable`] has wrappers around [`Geometry`](crate::geometry::Geometry),
[`GeometryRef`](crate::geometry::GeometryRef) and
[`GeometryCow`](crate::geometry::GeometryCow) which contain the enum and a
[`Style`] struct so they can be drawn onto a [`cairo::Context`].
- [`intersection`] offers `draw` implementations for
[`Intersection`](crate::composite::Intersection)s.
 */

pub mod geometry;
pub use geometry::*;

pub mod intersection;
pub use intersection::*;

pub mod drawable;
pub use drawable::*;
