/*!
[`ArcSegment`]: crate::segment::ArcSegment
[`LineSegment`]: crate::segment::LineSegment
[`Segment`]: crate::segment::Segment
[`Polysegment`]: crate::polysegment::Polysegment
[`Contour`]: crate::contour::Contour
[`Shape`]: crate::shape::Shape
[`Geometry`]: crate::geometry::Geometry
[`intersections`]: crate::geometry::GeometryRef::intersections
[`Primitive`]: crate::primitive::Primitive
[`Composite`]: crate::composite::Composite
[`Transformation`]: crate::Transformation
[`DEFAULT_EPSILON`]: crate::DEFAULT_EPSILON
[`DEFAULT_MAX_RELATIVE`]: crate::DEFAULT_MAX_RELATIVE
[crate_index]: crate
[draw]: crate::draw
[`Context`]: cairo::Context
[cairo]: cairo
[approxim]: approx
[serde]: serde
[`relative_eq`]: approx::relative_eq

A Rust library for 2D geometry: geometric objects, algorithms and visualization.
 */
#![cfg_attr(feature = "doc-images",
cfg_attr(all(),
doc = ::embed_doc_image::embed_image!("type_overview.svg", "docs/img/type_overview.svg"),
doc = ::embed_doc_image::embed_image!("shape.svg", "docs/img/shape.svg"),
doc = ::embed_doc_image::embed_image!("intersection_segments.svg", "docs/img/intersection_segments.svg"),
doc = ::embed_doc_image::embed_image!("intersection_composites.svg", "docs/img/intersection_composites.svg"),
))]
#![cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images not enabled**. Compile docs with `cargo doc --features 'doc-images'` and Rust version >= 1.54."
)]
#![doc = include_str!("../docs/main.md")]
// #![deny(missing_docs)]

use bounding_box::BoundingBox;

// Base
pub mod error;
pub mod prelude;

// Primitives
pub mod line;
pub mod primitive;
pub mod segment;

// Composites
pub mod composite;
pub mod contour;
pub mod polysegment;
pub mod shape;

// Container for primitives and composites
pub mod geometry;

#[cfg(feature = "cairo")]
pub mod draw;

/**
A reasonable default value for the absolute tolerance.

Various methods within this crate (such as the `contains_` or `intersections_`
methods from the [`Primitive`] and [`Composite`] traits) perform floating point
comparisons using the [`relative_eq`](approx::relative_eq) macro from the
[`approxim`](approx), which requires specifying an absolute tolerance `epsilon`
and a relative tolerance `max_relative`. If those tolerances aren't explicitly
provided via a [`ToleranceContext`], this constant is used for `epsilon` (
and [`DEFAULT_MAX_RELATIVE`] for `max_relative`). These defaults have been
chosen to provide "intuitive" behaviour of the aforementioned methods. See
[`ToleranceContext`] for more information.

The value of this constant is the square root of the machine precision
(`f64::EPSILON.sqrt()`).
*/
pub const DEFAULT_EPSILON: f64 = 0.000000014901161193847656_f64;

/**
A reasonable default value for the relative tolerance.

Various methods within this crate (such as the `contains_` or `intersections_`
methods from the [`Primitive`] and [`Composite`] traits) perform floating point
comparisons using the [`relative_eq`](approx::relative_eq) macro from the
[`approxim`](approx), which requires specifying an absolute tolerance `epsilon`
and a relative tolerance `max_relative`. If those tolerances aren't explicitly
provided via a [`ToleranceContext`], this constant is used for `max_relative` (
and [`DEFAULT_EPSILON`] for `epsilon`). These defaults have been
chosen to provide "intuitive" behaviour of the aforementioned methods. See
[`ToleranceContext`] for more information.
 */
pub const DEFAULT_MAX_RELATIVE: f64 = 1e-8;

/**
A tolerance context wrapper around a geometric object that allows specifying
custom tolerances for geometric operations such as intersection calculation.
Tolerance contexts are usually created with the `with_tolerance` method.

Various methods within this crate (such as the `covers_` or `intersections_`
methods from the [`Primitive`](crate::primitive::Primitive) and [`Composite`](crate::composite::Composite) traits) perform floating point comparisons.
Floating point numbers have finite precision and use a binary representation.
Consequently, many decimal numbers (such as `0.1`) cannot be represented exactly
as an `f64`, and arithmetic involving them can produce results that differ
slightly from the mathematically exact result:
For example, `assert_eq!(0.1 + 0.2, 0.3)` will panic!

Therefore, the aforementioned methods use the
[`relative_eq`](approx::relative_eq) macro from the [`approxim`](approx) crate
when comparing floats, which requires specifying an absolute tolerance `epsilon`
and a relative tolerance `max_relative`. These default to [`DEFAULT_EPSILON`]
and [`DEFAULT_MAX_RELATIVE`], but can be overwritten by using a
[`ToleranceContext`], which can be created by the `with_tolerance` methods of
the geometric objects:

```
use planar_geo::prelude::*;

let line_segment = LineSegment::new([0.0, 0.0], [1.0, 0.0]).expect("start and end point are different");

// From intuition, the point [0.5, 0.5] is not covered by the line segment:
assert!(!line_segment.covers(&[0.5, 0.5]));

// However, if the tolerance becomes sufficiently large, it is covered by the line segment:
assert!(line_segment.with_tolerance(1.0, DEFAULT_MAX_RELATIVE).covers(&[0.5, 0.5]));
```

The presented
[`LineSegment::with_tolerance`](crate::segment::LineSegment::with_tolerance)
method creates a [`ToleranceContext<LineSegment>`] which can be used for one-off
calculations as shown above, but can also be reused across multiple operations:

```
use planar_geo::prelude::*;

let line_segment = LineSegment::new([0.0, 0.0], [1.0, 0.0]).expect("start and end point are different");

let tol_context = line_segment.with_tolerance(1.0, DEFAULT_MAX_RELATIVE);
assert!(tol_context.covers(&[0.5, 0.5]));
assert!(tol_context.covers(&[0.5, -0.5]));
assert!(!tol_context.covers(&[0.5, 1.5])); // This point is still not covered
```

The default values for [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_RELATIVE`]
have been chosen to provide intuitive behaviour without making the tolerances
unnecessarily large. The following example demonstrates why simply using a zero
tolerance can lead to strange results:

```
use planar_geo::prelude::*;

// ArcSegment::same_circle checks if center and radius of two ArcSegments are
// identical (i.e., if they lay on the same underlying circle).

let as1 = ArcSegment::new([0.1, 0.1], 0.1, 0.0, 1.0).expect("valid inputs");
let as2 = ArcSegment::from_start_middle_stop([0.2, 0.1], [0.1, 0.2], [0.0, 0.1]).expect("valid inputs");

// Intuitively, both arc segments are on the same circle with center [0.1, 0.1]
// and radius 0.1. With the default tolerances, this is indeed the case.
assert!(as1.same_circle(&as2));

// Using zero tolerances runs into the problem that 0.1 is not exactly
// representable in f64. Since the center of as2 needs to be calculated from
// three points, this leads to slight differences:
assert!(!as1.with_tolerance(0.0, 0.0).same_circle(&as2));
```

Unless there is a specific reason for customizing the tolerances, it is
recommended to simply use the default tolerances (which happens automatically
unless `with_tolerance` is inserted between the geometric object and the
method).

As stated in the [`approxim`](approx) documentation: "Floating point is hard!".
The following links, taken directly from the [`approxim`](approx) crate
documentation, provide more informatio regarding the behaviour of floating point
numbers, particularly when comparing them:
- [Comparing Floating Point Inners, 2012 Edition](https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/)
- [The Floating Point Guide - Comparison](https://floating-point-gui.de/errors/comparison/)
- [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
*/
pub struct ToleranceContext<'p, T> {
    /// A reference to the geometric object to which this context applies.
    ///
    /// Usually, the context is created by calling `with_tolerance` on this
    /// object.
    pub inner: &'p T,
    /// The absolute tolerance used within the context.
    pub epsilon: f64,
    /// The relative tolerance used within the context.
    pub max_relative: f64,
}

pub(crate) mod private {
    /// Sealed trait for `Composite`, `Primitive` and `WithTolerance` to prevent
    /// implementation of the traits outside the crate.
    pub trait Sealed {}
}

pub trait WithTolerance: Sized + private::Sealed {
    /// Wraps `self` in a [`ToleranceContext`] using the specified `epsilon` and
    /// `max_relative` tolerances.
    ///
    /// The [`ToleranceContext`] applies these tolerances to floating-point
    /// comparisons performed by geometric operations, such as finding
    /// intersections. See [`ToleranceContext`] for details and examples.
    fn with_tolerance<'a>(&'a self, epsilon: f64, max_relative: f64) -> ToleranceContext<'a, Self> {
        ToleranceContext {
            inner: self,
            epsilon,
            max_relative,
        }
    }
}

/**
Affine transformations for geometric types.

All geometric types within this crate as well as the basic "point"
`[f64;2]` and [`bounding_box::BoundingBox`] types implement this trait to allow
for easy affine transformations. All examples in the docstrings of the
ndividual trait methods are for the point type, because all other geometric
types are based on it. Hence, their implementation basically just delegates to
`impl Transformation for [f64; 2]` for all their points.
 */
pub trait Transformation {
    /**
    Translates `self` by the given `shift`.

    ```
    use planar_geo::prelude::Transformation;

    let mut pt = [1.0, 1.0];
    pt.translate([1.0, -1.0]);
    assert_eq!(pt, [2.0, 0.0]);
    ```
     */
    fn translate(&mut self, shift: [f64; 2]);

    /**
    Rotates `self` around the `center` by the given `angle` (in rad).

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::Transformation;

    let mut pt = [1.0, 1.0];
    pt.rotate([1.0, -1.0], PI);
    approx::assert_abs_diff_eq!(pt, [1.0, -3.0], epsilon = 1e-15);
    ```
     */
    fn rotate(&mut self, center: [f64; 2], angle: f64);

    /**
    Scales `self` by `factor` with respect to the origin `[0.0, 0.0]`.

    ```
    use planar_geo::prelude::Transformation;

    let mut pt = [1.0, 1.0];
    pt.scale(2.0);
    assert_eq!(pt, [2.0, 2.0]);
    ```
     */
    fn scale(&mut self, factor: f64);

    /**
    Mirrors `self` about a line defined by two points.

    ```
    use planar_geo::prelude::Transformation;

    let mut pt = [1.0, 1.0];
    pt.line_reflection([0.0, 0.0], [0.0, 2.0]); // Vertical line
    approx::assert_abs_diff_eq!(pt, [-1.0, 1.0], epsilon = 1e-15);
    pt.line_reflection([0.0, -1.0], [2.0, -1.0]); // Horizontal line
    approx::assert_abs_diff_eq!(pt, [-1.0, -3.0], epsilon = 1e-15);
    ```
     */
    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> ();

    /**
    Mirrors `self` about a `point`. This operation is equivalent to a rotation
    around the point with the angle PI.

    ```
    use planar_geo::prelude::Transformation;

    let mut pt = [1.0, 1.0];
    pt.point_reflection([0.0, 0.0]);
    approx::assert_abs_diff_eq!(pt, [-1.0, -1.0], epsilon = 1e-15);
    pt.point_reflection([0.0, 2.0]);
    approx::assert_abs_diff_eq!(pt, [1.0, 5.0], epsilon = 1e-15);
    ```
     */
    fn point_reflection(&mut self, point: [f64; 2]) -> () {
        self.rotate(point, std::f64::consts::PI);
    }
}

impl Transformation for [f64; 2] {
    fn translate(&mut self, shift: [f64; 2]) {
        *self = [self[0] + shift[0], self[1] + shift[1]];
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let t = Rotation2::new(angle);
        let pt = [self[0] - center[0], self[1] - center[1]];
        let mut p = t * pt;
        p.translate([center[0], center[1]]);
        *self = p;
    }

    fn scale(&mut self, factor: f64) {
        *self = [self[0] * factor, self[1] * factor];
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        // Treatment of special case: vertical line
        if start[0] == stop[0] {
            *self = [-self[0] + 2.0 * start[0], self[1]];

        // Treatment of special case: Horizontal line
        } else if start[1] == stop[1] {
            *self = [self[0], -self[1] + 2.0 * start[1]];
        // All other cases
        } else {
            // Solve the line equation
            let m = (stop[1] - start[1]) / (stop[0] - start[0]);
            let c = (stop[0] * start[1] - start[0] * stop[1]) / (stop[0] - start[0]);

            let d = (self[0] + (self[1] - c) * m) / (1.0 + m.powi(2));
            *self = [2.0 * d - self[0], 2.0 * d * m - self[1] + 2.0 * c];
        }
    }
}

impl Transformation for BoundingBox {
    fn translate(&mut self, shift: [f64; 2]) {
        BoundingBox::translate(self, shift);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let mut ll = [self.xmin(), self.ymin()];
        let mut ul = [self.xmin(), self.xmax()];
        let mut lr = [self.xmax(), self.ymin()];
        let mut ur = [self.xmax(), self.xmax()];
        ll.rotate(center, angle);
        ul.rotate(center, angle);
        lr.rotate(center, angle);
        ur.rotate(center, angle);
        *self = BoundingBox::from_points([ll, ul, lr, ur].into_iter())
            .expect("supplied more than zero points")
    }

    fn scale(&mut self, factor: f64) {
        BoundingBox::scale(self, factor);
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        let mut ll = [self.xmin(), self.ymin()];
        let mut ul = [self.xmin(), self.xmax()];
        let mut lr = [self.xmax(), self.ymin()];
        let mut ur = [self.xmax(), self.xmax()];
        ll.line_reflection(start, stop);
        ul.line_reflection(start, stop);
        lr.line_reflection(start, stop);
        ur.line_reflection(start, stop);
        *self = BoundingBox::from_points([ll, ul, lr, ur].into_iter())
            .expect("supplied more than zero points")
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CentroidData {
    pub(crate) area: f64,
    pub(crate) x: f64,
    pub(crate) y: f64,
}

impl CentroidData {
    /**
    Calculates the common centroid of `self` and `other`.
     */
    pub(crate) fn union(&self, other: &Self) -> Self {
        let area = self.area + other.area;
        if area == 0.0 {
            let x = 0.0;
            let y = 0.0;
            return Self { area, x, y };
        } else {
            let x = (self.area * self.x + other.area * other.x) / area;
            let y = (self.area * self.y + other.area * other.y) / area;
            return Self { area, x, y };
        }
    }

    /**
    Calculates the common centroid of `self` and `other`, but subtract `other`
    from `self` (i.e. the original contour of `other` is subtracted from `self`)
     */
    pub(crate) fn subtract(&self, other: &Self) -> Self {
        let area = self.area - other.area;
        if area == 0.0 {
            let x = 0.0;
            let y = 0.0;
            return Self { area, x, y };
        } else {
            let x = (self.area * self.x - other.area * other.x) / area;
            let y = (self.area * self.y - other.area * other.y) / area;
            return Self { area, x, y };
        }
    }
}

impl From<CentroidData> for [f64; 2] {
    fn from(value: CentroidData) -> Self {
        return [value.x, value.y];
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Rotation2 {
    sin: f64,
    cos: f64,
}

impl Rotation2 {
    pub(crate) fn new(angle: f64) -> Self {
        return Self {
            sin: angle.sin(),
            cos: angle.cos(),
        };
    }
}

impl std::ops::Mul<[f64; 2]> for Rotation2 {
    type Output = [f64; 2];

    fn mul(self, rhs: [f64; 2]) -> Self::Output {
        return [
            rhs[0] * self.cos - rhs[1] * self.sin,
            rhs[0] * self.sin + rhs[1] * self.cos,
        ];
    }
}
