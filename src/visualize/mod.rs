/*!
This module contains the [`imp`] and [`intersection`] submodules. The former
holds the `draw` implementations for all [`Primitive`](crate::primitive::Primitive)
and [`Composite`](crate::composite::Composite) types of this crate, while the
latter offers `draw` implementations for [`Intersection`](crate::composite::Intersection)s.
See their respective module documentations for more.
 */

pub mod imp;
pub use imp::*;

pub mod intersection;
pub use intersection::*;
