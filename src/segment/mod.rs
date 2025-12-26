/*!
A module containing the [`Segment`] enum and various helper structs.

The [`Segment`] enum contains the segment variants available in this crate, such
as the [`ArcSegment`] or the [`LineSegment`]. A segment is a directed line
connecting a start to an end point. All types in this crate are two-dimensional
(planar), meaning that start and end points consist of two values (`[f64; 2]`),
where the first element is interpreted as the `x` value and the second element
is interpreted as the `y` value in the underlying cartesian coordinate system.

Multiple connected segments form a
[`SegmentChain`](crate::segment_chain::SegmentChain), which in turn is
used to define [`Contour`](crate::contour::Contour)s and by extension
[`Shape`](crate::shape::Shape)s, the [`Composite`](crate::composite::Composite)
types of this crate. However, segments can also used on their own to e.g.
calculate properties or finding intersections. They implement the [`Primitive`]
trait, which provides various methods shared with other primitive geometric
types such as [`Line`](crate::line::Line)s and points (`[f64; 2]`).
*/

use bounding_box::BoundingBox;
pub mod arc_segment;
pub mod line_segment;
pub use arc_segment::*;
pub use line_segment::*;

#[cfg(feature = "serde")]
use serde::Serialize;

#[cfg(feature = "serde")]
use deserialize_untagged_verbose_error::DeserializeUntaggedVerboseError;

use crate::{
    Transformation,
    primitive::{Primitive, PrimitiveIntersections},
};

/**
A directed connection between a start and an end point, which can take different
paths depending on the variant. All methods of [`Segment`] delegate to the
corresponding function of the underlying variant.

*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/example_segments.svg \"Different segment types\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/example_segments.svg)"
)]
/**

See the [module-level](crate::segment) docstring for more.

# Serialization and deserialization

When the `serde` feature is enabled, segments can be serialized and deserialized
using the [`serde`] crate. The enum is treated as untagged, meaning that the
serialized representation of a segment is identical to that of the segment
variant it contains ([`LineSegment`] or [`ArcSegment`]). See [`serde`]s
documentation for more.
 */
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, DeserializeUntaggedVerboseError))]
#[cfg_attr(feature = "serde", serde(untagged))]
pub enum Segment {
    /// A straight line segment between a start and an end point. See the
    /// docstring of [`LineSegment`].
    LineSegment(LineSegment),
    /// An arc segment between a start and an end point. See the docstring of
    /// [`ArcSegment`].
    ArcSegment(ArcSegment),
}

impl Segment {
    /**
    Returns the start point of the underlying segment variant.
     */
    pub fn start(&self) -> [f64; 2] {
        match self {
            Segment::LineSegment(line_segment) => line_segment.start(),
            Segment::ArcSegment(arc_segment) => arc_segment.start(),
        }
    }

    /**
    Returns the end / stop point of the underlying segment variant.
     */
    pub fn end(&self) -> [f64; 2] {
        match self {
            Segment::LineSegment(line_segment) => line_segment.stop(),
            Segment::ArcSegment(arc_segment) => arc_segment.stop(),
        }
    }

    /**
    Returns the end / stop point of the underlying segment variant. This is an
    alias for [`Segment::end`].
     */
    pub fn stop(&self) -> [f64; 2] {
        return self.end();
    }

    /**
    Returns the number of points of the underlying variant.
     */
    pub fn number_points(&self) -> usize {
        match self {
            Segment::LineSegment(line_segment) => line_segment.number_points(),
            Segment::ArcSegment(arc_segment) => arc_segment.number_points(),
        }
    }

    /**
    Returns the length of the underlying line / arc segment.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let ls = LineSegment::new([0.0, 0.0], [0.0, 2.0], 0.0, 0).unwrap();
    assert_eq!(ls.length(), 2.0);
    let s = Segment::from(ls);
    assert_eq!(s.length(), 2.0);

    let arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.0,
            PI,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
    approx::assert_abs_diff_eq!(arc.length(), 2.0 * PI);
    let s = Segment::from(arc);
    approx::assert_abs_diff_eq!(s.length(), 2.0 * PI);
    ```
     */
    pub fn length(&self) -> f64 {
        match self {
            Segment::LineSegment(line_segment) => line_segment.length(),
            Segment::ArcSegment(arc_segment) => arc_segment.length(),
        }
    }

    /**
    Reverses the underlying segment, i.e. switching the segment direction.
     */
    pub fn reverse(&mut self) {
        match self {
            Segment::LineSegment(line_segment) => line_segment.reverse(),
            Segment::ArcSegment(arc_segment) => arc_segment.reverse(),
        }
    }

    /**
    Returns a point on the segment defined by its normalized position on it.

    For example, `normalized = 0` returns the start point, `normalized = 1`
    returns the end point and `normalized = 0.5` returns the middle point of the
    segment. The input `normalized` is clamped to [0, 1].

    # Examples

    ```
    use planar_geo::prelude::*;
    let ls = LineSegment::new([0.0, 0.0], [0.0, 2.0], 0.0, 0).unwrap();
    let s = Segment::from(ls);

    assert_eq!(s.segment_point(0.0), s.start());
    assert_eq!(s.segment_point(-10.0), s.start());

    assert_eq!(s.segment_point(1.0), s.stop());
    assert_eq!(s.segment_point(10.0), s.stop());

    assert_eq!(s.segment_point(0.5), [0.0, 1.0]);
    ```
    */
    pub fn segment_point(&self, normalized: f64) -> [f64; 2] {
        match self {
            Segment::LineSegment(line_segment) => line_segment.segment_point(normalized),
            Segment::ArcSegment(arc_segment) => arc_segment.segment_point(normalized),
        }
    }

    /**
    Returns the points of a polygon chain which approximates `self`. The number
    of points is defined by the [`SegmentPolygonizer`] (see its docstring). The
    points are regularily distributed over the segment, which means that two
    subsequent points always have the same euclidian distance from each other.

    # Examples

    ```
    use std::f64::consts::FRAC_PI_2;
    use planar_geo::prelude::*;

    // Approximate an arc segment from 0 to 90 degrees by four straight
    // segments (five points). The returned points are regularily distributed
    // over the arc.
    let arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            1.0,
            0.0,
            FRAC_PI_2,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
    let s = Segment::from(arc);
    let mut iter = s.polygonize(SegmentPolygonizer::NumberSegments(4));

    assert_eq!(iter.next(), Some(s.segment_point(0.0)));
    assert_eq!(iter.next(), Some(s.segment_point(0.25)));
    assert_eq!(iter.next(), Some(s.segment_point(0.5)));
    assert_eq!(iter.next(), Some(s.segment_point(0.75)));
    assert_eq!(iter.next(), Some(s.segment_point(1.0)));
    assert!(iter.next().is_none());
    ```
     */
    pub fn polygonize<'a>(&'a self, polygonizer: SegmentPolygonizer) -> PolygonPointsIterator<'a> {
        match self {
            Segment::LineSegment(line_segment) => line_segment.polygonize(polygonizer),
            Segment::ArcSegment(arc_segment) => arc_segment.polygonize(polygonizer),
        }
    }

    /**
    Returns the centroid (center of mass) of the segment.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.0,
            PI,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
    let s = Segment::from(arc);
    approx::assert_abs_diff_eq!(s.centroid(), [0.0, 0.8488263631567751]);
    ```
     */
    pub fn centroid(&self) -> [f64; 2] {
        return crate::CentroidData::from(self).into();
    }
}

impl Transformation for Segment {
    fn translate(&mut self, shift: [f64; 2]) {
        match self {
            Segment::LineSegment(obj) => obj.translate(shift),
            Segment::ArcSegment(obj) => obj.translate(shift),
        }
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        match self {
            Segment::LineSegment(obj) => obj.rotate(center, angle),
            Segment::ArcSegment(obj) => obj.rotate(center, angle),
        }
    }

    fn scale(&mut self, factor: f64) {
        match self {
            Segment::LineSegment(obj) => obj.scale(factor),
            Segment::ArcSegment(obj) => obj.scale(factor),
        }
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        match self {
            Segment::LineSegment(obj) => obj.line_reflection(start, stop),
            Segment::ArcSegment(obj) => obj.line_reflection(start, stop),
        }
    }
}

impl crate::primitive::private::Sealed for Segment {}

impl Primitive for Segment {
    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        match self {
            Segment::LineSegment(s) => s.contains_point(point, epsilon, max_ulps),
            Segment::ArcSegment(s) => s.contains_point(point, epsilon, max_ulps),
        }
    }

    fn intersections_line(
        &self,
        line: &crate::line::Line,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        match self {
            Segment::LineSegment(s) => s.intersections_line(line, epsilon, max_ulps),
            Segment::ArcSegment(s) => s.intersections_line(line, epsilon, max_ulps),
        }
    }

    fn intersections_line_segment(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        match self {
            Segment::LineSegment(s) => {
                s.intersections_line_segment(line_segment, epsilon, max_ulps)
            }
            Segment::ArcSegment(s) => s.intersections_line_segment(line_segment, epsilon, max_ulps),
        }
    }

    fn intersections_arc_segment(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        match self {
            Segment::LineSegment(s) => s.intersections_arc_segment(arc_segment, epsilon, max_ulps),
            Segment::ArcSegment(s) => s.intersections_arc_segment(arc_segment, epsilon, max_ulps),
        }
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        other.intersections_segment(self, epsilon, max_ulps)
    }
}

impl From<&Segment> for BoundingBox {
    fn from(value: &Segment) -> BoundingBox {
        match value {
            Segment::LineSegment(obj) => obj.into(),
            Segment::ArcSegment(obj) => obj.into(),
        }
    }
}

impl From<&Segment> for crate::CentroidData {
    fn from(value: &Segment) -> Self {
        match value {
            Segment::LineSegment(line_segment) => line_segment.into(),
            Segment::ArcSegment(arc_segment) => arc_segment.into(),
        }
    }
}

impl From<LineSegment> for Segment {
    fn from(value: LineSegment) -> Self {
        return Self::LineSegment(value);
    }
}

impl From<ArcSegment> for Segment {
    fn from(value: ArcSegment) -> Self {
        return Self::ArcSegment(value);
    }
}

/**
This enum defines how many points should be in the polygonized representation
of a segment (created by [`Segment::polygonize`]). Depending on the selected
variant and its parametrization, a different number of points is created:
```
use std::f64::consts::PI;
use planar_geo::prelude::*;

let arc = ArcSegment::from_center_radius_start_offset_angle(
    [0.0, 0.0],
    1.0,
    1.25 * PI,
    0.75 * PI,
    DEFAULT_EPSILON,
    DEFAULT_MAX_ULPS,
)
.unwrap();
let s = Segment::from(arc);

// One segment (red line in the image below this code snippet)
let iter = s.polygonize(SegmentPolygonizer::MaximumSegmentLength(10.0));
assert_eq!(iter.count(), 2);

// Two segments (green line in the image below this code snippet)
let iter = s.polygonize(SegmentPolygonizer::MaximumAngle(0.5 * PI));
assert_eq!(iter.count(), 3);

// Three segments (blue line in the image below this code snippet)
let iter = s.polygonize(SegmentPolygonizer::NumberSegments(3));
assert_eq!(iter.count(), 4);
```
*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/polygonized_arc.svg \"Line style comparison\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/polygonized_arc.svg)"
)]
#[derive(Debug, Clone, Copy)]
pub enum SegmentPolygonizer {
    /**
    Number of line segments of the polygon chains given explicitly
    -> Number of points is equal to this plus one (end point of the last
    segment). If the number of line segments is specified to be 0, only the
    start point of the input is returned.
     */
    NumberSegments(usize),
    /**
    The given number is the upper limit for the allowed length of the line
    segments. The algorithm then calculates the fewest number of line segments
    possible where this constraint is fulfilled and returns the respective
    points. If the given value is zero or negative, only the start point is
    returned.
     */
    MaximumSegmentLength(f64),
    /**
    If applied to an [`ArcSegment`], this constraint limits the maximum angle
    which can be convered by a single line segment. The algorithm then
    calculates the fewest number of line segments possible where this constraint
    is fulfilled and returns the respective points. If the given value is
    zero or negative, only the start point is returned.
    If applied to a [`LineSegment`], only the start point is returned.
     */
    MaximumAngle(f64),
}

impl Default for SegmentPolygonizer {
    fn default() -> Self {
        return SegmentPolygonizer::NumberSegments(1);
    }
}

/**
The "borrowed" version of [`Segment`] which as the same variants as its brother,
but with references instead of owned instances of the underlying segment types.
 */
#[derive(Clone, Debug)]
pub(crate) enum SegmentRef<'a> {
    /// A reference to a [`LineSegment`].
    LineSegment(&'a LineSegment),
    /// A reference to an [`ArcSegment`].
    ArcSegment(&'a ArcSegment),
}

/**
An iterator created by calling [`Segment::polygonize`]. See the docstring of
[`SegmentPolygonizer`] for more.
 */
#[derive(Clone, Debug)]
pub struct PolygonPointsIterator<'a> {
    index: usize,
    num_segs: usize,
    segment: SegmentRef<'a>,
}

impl<'a> PolygonPointsIterator<'a> {
    pub(super) fn start_at_second_point(&mut self) {
        self.index = 1;
    }
    pub(super) fn skip_last_vertex(&mut self) {
        if self.num_segs > 0 {
            self.num_segs -= 1;
        }
    }
}

impl<'a> Iterator for PolygonPointsIterator<'a> {
    type Item = [f64; 2];

    fn next(&mut self) -> Option<Self::Item> {
        if self.index == self.num_segs + 1 {
            return None;
        }
        let point = match self.segment {
            SegmentRef::LineSegment(segment) => {
                let x = segment.start()[0]
                    + self.index as f64 / self.num_segs as f64
                        * (segment.stop()[0] - segment.start()[0]);
                let y = segment.start()[1]
                    + self.index as f64 / self.num_segs as f64
                        * (segment.stop()[1] - segment.start()[1]);
                [x, y]
            }
            SegmentRef::ArcSegment(segment) => {
                let angle = segment.start_angle()
                    + self.index as f64 / self.num_segs as f64 * segment.offset_angle();
                let mut center = segment.center();
                center.translate([
                    segment.radius() * angle.cos(),
                    segment.radius() * angle.sin(),
                ]);
                center
            }
        };

        self.index += 1;

        return Some(point);
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let length = self.num_segs + 1;
        (length, Some(length))
    }
}
