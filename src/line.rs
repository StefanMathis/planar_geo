/*!
This module contains the [`Line`] struct, which represents a two-dimensional
mathematical line - an infinitely long object with no width or curvature.

In the context of this crate, a line is treated as a "primitive" geometric type.
Unlike the "segment" types ([`Segment`](crate::segment::Segment) and its
variants [`ArcSegment`](crate::segment::ArcSegment) and
[`LineSegment`](crate::segment::LineSegment)), it is not used in defining more
complex "composite" types (such as the
[`Polysegment`](crate::polysegment::Polysegment)). Its main purpose is to
serve as a tool for calculations, for example for intersection calculation.

See the docstring of [`Line`] for more information.
 */

use std::f64::{INFINITY, NEG_INFINITY};

use approx::relative_eq;
use bounding_box::{BoundingBox, ToBoundingBox};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE, Rotation2, Transformation,
    primitive::{Primitive, PrimitiveIntersections},
};

/**
A mathematical Line which has an infinite length, but no width.

A line can be represented by the equation `a*x + b*y + c = 0`. The three fields
[`Line::a`], [`Line::b`] and [`Line::c`] correspond to these coefficients.

The main purpose of this type is being a tool for calculations. For example,
the intersection of two [`LineSegment`](crate::segment::LineSegment)s (straight
connections between two points of finite length) can be calculated by first
calculating the intersection of the corresponding infinite lines and then
checking whether the found intersection point is actually covered by both
segments.

Obviously, a [`Line`] object can be directly created by providing its three
coefficients. Additionally, it is also possible to derive a [`Line`] from a
point it goes through and its angle ([`Line::from_point_angle`]) and from a
two-point representation ([`Line::from_two_points`]).
Because the [`Line`] is closely related to the
[`LineSegment`](crate::segment::LineSegment), a [`From`] implementation exists.

The [crate docstring](crate::line) describes the relationship between the
[`Line`] and the other geometric types provided by this crate.

# Serialization and deserialization

This struct can be serialized and deserialized from its three fields
[`Line::a`], [`Line::b`] and [`Line::c`] using the [`serde`] crate if the
`serde` feature is enabled.
 */
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct Line {
    /// Coefficient `a` in the line formula `a*x + b*y + c = 0`.
    pub a: f64,
    /// Coefficient `b` in the line formula `a*x + b*y + c = 0`.
    pub b: f64,
    /// Coefficient `c` in the line formula `a*x + b*y + c = 0`.
    pub c: f64,
}

impl Line {
    /**
    Creates a [`Line`] from its three coefficients. This is an alias for using
    the literal struct construction syntax `Line { a, b, c }`.
     */
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        return Self { a, b, c };
    }

    /**
    Creates a [`Line`] from a point it goes through and its angle (relative to
    the `x` / horizontal axis).

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    // A line with a 45° angle going through [2.0, 0.0].
    let line = Line::from_point_angle([2.0, 0.0], 0.25 * PI);

    assert!(line.covers_point([1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(line.covers_point([2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(line.covers_point([3.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn from_point_angle(pt: [f64; 2], angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        let a = -sin;
        let b = cos;
        let c = -(a * pt[0] + b * pt[1]);
        Self { a, b, c }
    }

    /**
    Creates a [`Line`] from two points. This constructor fails if the points are
    identical.

    # Examples

    ```
    use planar_geo::prelude::*;

    // A line with a 45° angle going through [2.0, 0.0].
    let line = Line::from_two_points([2.0, 0.0], [3.0, 1.0]).expect("points not identical");

    assert!(line.covers_point([1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(line.covers_point([2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(line.covers_point([3.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn from_two_points(pt1: [f64; 2], pt2: [f64; 2]) -> Result<Self, crate::error::Error> {
        Ok(crate::segment::LineSegment::new(pt1, pt2)?.into())
    }

    /**
    Returns `true` if the two given lines are parallel within the tolerance band
    defined by the absolute tolerance `epsilon` and the relative tolerance
    `max_relative`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let line_1 = Line::from_point_angle([2.0, 0.0], 1.0);
    let line_2 = Line::from_point_angle([2.0, 1.0], 1.0);
    let line_3 = Line::from_point_angle([2.0, 0.0], -1.0);
    assert!(line_1.parallel(&line_2, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!line_1.parallel(&line_3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!line_2.parallel(&line_3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn parallel(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        let zn = det(self.a, self.b, other.a, other.b);
        return relative_eq!(zn, 0.0, epsilon = epsilon, max_relative = max_relative);
    }

    /**
    Returns `true` if the two given lines are identical within the tolerance
    band defined by the absolute tolerance `epsilon` and the relative tolerance
    `max_relative`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let line_1 = Line::from_point_angle([2.0, 0.0], 0.0);
    let line_2 = Line::from_point_angle([-3.0, 0.0], 0.0);
    let line_3 = Line::from_point_angle([2.0, 1.0], 0.0);
    assert!(line_1.identical(&line_2, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!line_1.identical(&line_3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!line_2.identical(&line_3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn identical(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        return relative_eq!(
            det(self.a, self.b, other.a, other.b),
            0.0,
            epsilon = epsilon,
            max_relative = max_relative
        ) && relative_eq!(
            det(self.a, self.c, other.a, other.c),
            0.0,
            epsilon = epsilon,
            max_relative = max_relative
        ) && relative_eq!(
            det(self.b, self.c, other.b, other.c),
            0.0,
            epsilon = epsilon,
            max_relative = max_relative
        );
    }
}

impl PartialEq for Line {
    fn eq(&self, other: &Self) -> bool {
        return self.identical(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
    }
}

impl From<crate::segment::LineSegment> for Line {
    fn from(l: crate::segment::LineSegment) -> Self {
        return Line::from(&l);
    }
}

impl From<[[f64; 2]; 2]> for Line {
    fn from(pts: [[f64; 2]; 2]) -> Self {
        let [pt1, pt2] = pts;
        let a = pt2[1] - pt1[1];
        let b = pt1[0] - pt2[0];
        let c = pt2[0] * pt1[1] - pt1[0] * pt2[1];
        return Self { a, b, c };
    }
}

impl From<&crate::segment::LineSegment> for Line {
    fn from(l: &crate::segment::LineSegment) -> Self {
        return [l.start(), l.stop()].into();
    }
}

impl ToBoundingBox for Line {
    fn bounding_box(&self) -> BoundingBox {
        if self.a == 0.0 {
            // Horizontal line
            let y = self.c / self.b;
            return BoundingBox::new(NEG_INFINITY, INFINITY, y, y);
        } else if self.b == 0.0 {
            // Vertical line
            let x = self.c / self.a;
            return BoundingBox::new(x, x, NEG_INFINITY, INFINITY);
        } else {
            return BoundingBox::new(NEG_INFINITY, INFINITY, NEG_INFINITY, INFINITY);
        }
    }
}

fn det(a: f64, b: f64, c: f64, d: f64) -> f64 {
    return a * d - b * c;
}

impl Transformation for Line {
    fn translate(&mut self, shift: [f64; 2]) {
        self.c = self.c - (self.a * shift[0] + self.b * shift[1]);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let t = Rotation2::new(angle);
        let ab_r = t * [self.a, self.b];
        self.c = self.c + self.a * center[0] + self.b * center[1]
            - ab_r[0] * center[0]
            - ab_r[1] * center[1];
        self.a = ab_r[0];
        self.b = ab_r[1];
    }

    fn scale(&mut self, factor: f64) {
        self.a /= factor;
        self.b /= factor;
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        let [x1, y1] = start;
        let [x2, y2] = stop;

        // Direction of reflection line
        let dx = x2 - x1;
        let dy = y2 - y1;

        // Normal of reflection line (not unit yet)
        let mut mx = dy;
        let mut my = -dx;

        // Normalize m
        let len = (mx * mx + my * my).sqrt();
        mx /= len;
        my /= len;

        // Extract original line components
        let (a, b, c) = (self.a, self.b, self.c);

        // Step 1: translate so reflection line passes through origin
        // The point (x1, y1) goes to (0, 0),
        // so c becomes:
        let c0 = c + a * x1 + b * y1;

        // Step 2: reflect the normal vector
        let dot = a * mx + b * my;
        let a2 = a - 2.0 * dot * mx;
        let b2 = b - 2.0 * dot * my;

        // Step 3: translate back
        let c2 = c0 - (a2 * x1 + b2 * y1);

        self.a = a2;
        self.b = b2;
        self.c = c2;
    }
}

impl crate::primitive::private::Sealed for Line {}

impl Primitive for Line {
    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_relative: f64) -> bool {
        // Solution from https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        return relative_eq!(
            (self.a * point[0] + self.b * point[1] + self.c).abs()
                / (self.a.powi(2) + self.b.powi(2)).sqrt(),
            0.0,
            epsilon = epsilon,
            max_relative = max_relative
        );
    }

    fn covers_arc_segment(
        &self,
        _arc_segment: &crate::prelude::ArcSegment,
        _epsilon: f64,
        _max_relative: f64,
    ) -> bool {
        return false;
    }

    fn covers_line_segment(
        &self,
        line_segment: &crate::segment::LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> bool {
        let other = Self::from(line_segment);
        return self.covers_line(&other, epsilon, max_relative);
    }

    fn covers_line(&self, line: &Line, epsilon: f64, max_relative: f64) -> bool {
        return self.identical(line, epsilon, max_relative);
    }

    fn intersections_line(
        &self,
        line: &Line,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        fn det(a: f64, b: f64, c: f64, d: f64) -> f64 {
            return a * d - b * c;
        }

        /*
        Calculate the intersection between two lines. If the lines are parallel, return None instead.
        https://cp-algorithms.com/geometry/lines-intersection.html
        */
        let zn = det(self.a, self.b, line.a, line.b);
        if approx::relative_eq!(zn, 0.0, epsilon = epsilon, max_relative = max_relative) {
            return PrimitiveIntersections::Zero;
        } else {
            let x = -det(self.c, self.b, line.c, line.b) / zn;
            let y = -det(self.a, self.c, line.a, line.c) / zn;
            return PrimitiveIntersections::One([x, y]);
        }
    }

    fn intersections_line_segment(
        &self,
        line_segment: &crate::segment::LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        let other_line = Line::from(line_segment);
        match self.intersections_line(&other_line, epsilon, max_relative) {
            PrimitiveIntersections::Zero => PrimitiveIntersections::Zero,
            PrimitiveIntersections::One(pt) => {
                return line_segment.intersections_point(pt, epsilon, max_relative);
            }
            PrimitiveIntersections::Two(_) => {
                unreachable!("line-line can have either zero or one intersection")
            }
        }
    }

    fn intersections_arc_segment(
        &self,
        arc_segment: &crate::segment::ArcSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        arc_segment.intersections_line_circle(self.a, self.b, self.c, epsilon, max_relative)
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        return other.intersections_line(self, epsilon, max_relative);
    }
}
