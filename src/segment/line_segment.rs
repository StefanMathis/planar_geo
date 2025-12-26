/*!
A module containing the [`LineSegment`] struct, one of the available variants
for the [`Segment`](crate::segment::Segment) enum.

The [`LineSegment`] is one of the "primitive" geometric types defined in this
crate (implementing the [`Primitive`] type). Wrapped in the
[`Segment`](crate::segment::Segment) enum, it serves as a fundamental building
block of composite geometric types such as the
[`SegmentChain`](crate::segment_chain::SegmentChain).
See the module documentation of [segment](crate::segment) for more.

Most users should interact with this module through the [`LineSegment`] type
itself; see its documentation for details on construction, invariants, and
usage.
*/

use std::f64::INFINITY;

use crate::{
    Rotation2,
    primitive::{Primitive, PrimitiveIntersections},
};

use approx::ulps_eq;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::Transformation;
use bounding_box::BoundingBox;

/**
A straight, directed connection between a start and an end / stop point.

*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/example_line_segment.svg \"Line segment\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/example_line_segment.svg)"
)]
/**

By definition, start and end point must not be identical, because otherwise the
line segment would degenerate to a point. The finite length defined by those two
points is the main difference between this type and the
[`Line`](crate::line::Line) struct. It is trivial to convert a [`LineSegment`]
to a [`Line`](crate::line::Line) via the corresponding [`From`]
/ [`Into`] implementations, however a direct conversion in the other direction
is not possible because the start and end point information is lost.

# Serialization and deserialization

When the `serde` feature is enabled, line segments can be serialized and
deserialized using the [`serde`] crate.
 */
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct LineSegment {
    start: [f64; 2],
    stop: [f64; 2],
}

impl LineSegment {
    /**
    Creates a new [`LineSegment`] instance if the `start` and `stop` are not
    equal (within the tolerances defined by `epsilon` and `max_ulps`).

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    assert!(LineSegment::new([0.0, 0.0], [1.0, 0.0], 0.0, 0).is_ok());

    // Points are identical
    assert!(LineSegment::new([0.0, 0.0], [0.0, 0.0], 0.0, 0).is_err());

    // Points are identical within the defined tolerance
    assert!(LineSegment::new([0.0, 0.0], [1.0, 0.0], 1.0, 0).is_err());
    ```
     */
    pub fn new(
        start: [f64; 2],
        stop: [f64; 2],
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        if ulps_eq!(start, stop, epsilon = epsilon, max_ulps = max_ulps) {
            return Err(crate::error::ErrorType::PointsIdentical { start, stop }.into());
        } else {
            return Ok(LineSegment { start, stop });
        }
    }

    /**
    Creates a [`LineSegment`] from a `start` point, an `angle` and its `length`.
    If the `length` is zero, the points are identical and the construction
    fails.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    assert!(LineSegment::from_start_angle_length([0.0, 0.0], 0.0, 1.0, 0.0, 0).is_ok());
    assert!(LineSegment::from_start_angle_length([0.0, 0.0], 0.0, -1.0, 0.0, 0).is_ok());

    // Points are identical
    assert!(LineSegment::from_start_angle_length([0.0, 0.0], 0.0, 0.0, 0.0, 0).is_err());

    // Points are identical within the defined tolerance
    assert!(LineSegment::from_start_angle_length([0.0, 0.0], 0.0, 1.0, 1.0, 0).is_err());
    ```
     */
    pub fn from_start_angle_length(
        start: [f64; 2],
        angle: f64,
        length: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        let mut stop = start.clone();
        stop.translate([length * angle.cos(), length * angle.sin()]);
        return Self::new(start, stop, epsilon, max_ulps);
    }

    /**
    Returns the smallest x-value of `self`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], 0.0, 0).expect("points not identical");
    assert_eq!(ls.xmin(), 0.0);
    ```
     */
    pub fn xmin(&self) -> f64 {
        if self.start[0] < self.stop[0] {
            return self.start[0];
        } else {
            return self.stop[0];
        }
    }

    /**
    Returns the largest x-value of `self`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], 0.0, 0).expect("points not identical");
    assert_eq!(ls.xmax(), 1.0);
    ```
     */
    pub fn xmax(&self) -> f64 {
        if self.start[0] > self.stop[0] {
            return self.start[0];
        } else {
            return self.stop[0];
        }
    }

    /**
    Returns the smallest y-value of `self`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], 0.0, 0).expect("points not identical");
    assert_eq!(ls.ymin(), -1.0);
    ```
     */
    pub fn ymin(&self) -> f64 {
        if self.start[1] < self.stop[1] {
            return self.start[1];
        } else {
            return self.stop[1];
        }
    }

    /**
    Returns the largest y-value of `self`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // Successfull construction
    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], 0.0, 0).expect("points not identical");
    assert_eq!(ls.ymax(), 0.0);
    ```
     */
    pub fn ymax(&self) -> f64 {
        if self.start[1] > self.stop[1] {
            return self.start[1];
        } else {
            return self.stop[1];
        }
    }

    /**
    Returns the slope of `self`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    // 45° slope
    let line = LineSegment::new([0.0, 0.0], [1.0, -1.0], 0.0, 0).unwrap();
    assert_eq!(line.slope(), -1.0);

    // Infinite slope (vertical line)
    let line = LineSegment::new([0.0, 0.0], [0.0, -1.0], 0.0, 0).unwrap();
    assert!(line.slope().is_infinite());
    ```
     */
    pub fn slope(&self) -> f64 {
        if self.stop[0] == self.start[0] {
            return INFINITY * (self.stop[1] - self.start[1]).signum();
        } else {
            return (self.stop[1] - self.start[1]) / (self.stop[0] - self.start[0]);
        }
    }

    /**
    Returns the angle of `self` from [`LineSegment::start`] to
    [`LineSegment::stop`].

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::PI;

    // 45°
    let line = LineSegment::new([0.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON,
                                DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(line.angle(), 0.25 * PI);

    // 180°
    let line = LineSegment::new([1.0, 1.0], [-1.0, 1.0], DEFAULT_EPSILON,
                                DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(line.angle(), PI);

    // 225°
    let line = LineSegment::new([-2.0, -8.0], [-4.0, -10.0], DEFAULT_EPSILON,
                                DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(line.angle(), -0.75 * PI);

    // 315°
    let line = LineSegment::new([5.0, 0.0], [6.0, -1.0], DEFAULT_EPSILON,
                                DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(line.angle(), -0.25 * PI);
    ```
     */
    pub fn angle(&self) -> f64 {
        return (self.stop[1] - self.start[1]).atan2(self.stop[0] - self.start[0]);
    }

    /**
    Returns the euclidian distance of the `point` to `self`. The euclidian
    distance is the length of the shortest line which can be drawn between the
    `point` and any point on `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let line = LineSegment::new([0.0, 0.0], [1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    assert_eq!(line.euclidian_distance_to_point([-1.0, 0.0]), 1.0);
    approx::assert_abs_diff_eq!(line.euclidian_distance_to_point([2.0, 0.0]), 2.0f64.sqrt());

    // Point is on the line
    assert_eq!(line.euclidian_distance_to_point([0.5, -0.5]), 0.0);
    ```
     */
    pub fn euclidian_distance_to_point(&self, point: [f64; 2]) -> f64 {
        let start = self.start;
        let stop = self.stop;

        let dx = stop[0] - start[0];
        let dy = stop[1] - start[1];

        // r is the ratio of the projection of point onto self. If it is between 0.0 and 1.0, the shortest distance
        // from self to point is not from the end point
        let r =
            ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / (dx.powi(2) + dy.powi(2));
        if r <= 0.0 {
            return (start[0] - point[0]).hypot(start[1] - point[1]);
        }
        if r >= 1.0 {
            return (stop[0] - point[0]).hypot(stop[1] - point[1]);
        }
        let s = ((start[1] - point[1]) * dx - (start[0] - point[0]) * dy) / (dx * dx + dy * dy);
        s.abs() * dx.hypot(dy)
    }

    /**
    Returns the length of `self` calculated with the Pythagorean theorem.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    let ls = LineSegment::new([0.0, 0.0], [-2.0, 0.0], 0.0, 0).unwrap();
    assert_eq!(ls.length(), 2.0);
     ```
     */
    pub fn length(&self) -> f64 {
        return self.length_sq().sqrt();
    }

    /**
    Returns the squared length of `self` calculated with the Pythagorean theorem.

    This function is more efficient / faster than [`LineSegment::length`], since
    it does not require the expensive square root function. For example, if the
    longer between two [`LineSegment`]s should be found, using this function
    is the preferred way to do it:

    ```
    use planar_geo::segment::LineSegment;

    let ls1 = LineSegment::new([0.0, 0.0], [1.0, 0.0], 0.0, 0).unwrap();
    let ls2 = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).unwrap();

    assert_eq!(ls1.length_sq() > ls2.length_sq(), ls1.length() > ls2.length());
    ```
     */
    pub fn length_sq(&self) -> f64 {
        return (self.stop[0] - self.start[0]).powi(2) + (self.stop[1] - self.start[1]).powi(2);
    }

    /**
    Returns the start point of `self`.
     */
    pub fn start(&self) -> [f64; 2] {
        return self.start;
    }

    /**
    Returns the stop point of `self`.
     */
    pub fn stop(&self) -> [f64; 2] {
        return self.stop;
    }

    /**
    Returns the number of points in the segment. For a [`LineSegment`], this is
    always 2.
     */
    pub fn number_points(&self) -> usize {
        return 2;
    }

    /**
    Reverses the direction of `self` - i.e., exchange its start and stop point.

    ```
    use planar_geo::segment::LineSegment;

    let mut ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], 0.0, 0).unwrap();

    assert_eq!(ls.start(), [0.0, 0.0]);
    assert_eq!(ls.stop(), [1.0, 0.0]);

    ls.reverse();

    assert_eq!(ls.start(), [1.0, 0.0]);
    assert_eq!(ls.stop(), [0.0, 0.0]);
    ```
     */
    pub fn reverse(&mut self) {
        let tmp = self.start;
        self.start = self.stop;
        self.stop = tmp;
    }

    /**
    Returns a point on the segment defined by its normalized position on it.

    For example, `normalized = 0` returns the start point, `normalized = 1`
    returns the end point and `normalized = 0.5` returns the middle point of the
    segment. The input `normalized` is clamped to [0, 1].

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    let line = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).unwrap();

    // Middle point of the segment
    assert_eq!(line.segment_point(0.5), [1.0, 0.0]);
    assert_eq!(line.segment_point(-1.0), [0.0, 0.0]); // Start point
    assert_eq!(line.segment_point(1.5), [2.0, 0.0]); // Stop point
    ```
     */
    pub fn segment_point(&self, normalized: f64) -> [f64; 2] {
        let normalized = normalized.clamp(0.0, 1.0);
        return [
            self.start[0] + normalized * (self.stop[0] - self.start[0]),
            self.start[1] + normalized * (self.stop[1] - self.start[1]),
        ]
        .into();
    }

    /**
    Returns the centroid / center of mass of `self`. In case of the
    [`LineSegment`], this is equal to the middle of the segment.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    let line = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).unwrap();

    assert_eq!(line.centroid(), [1.0, 0.0]);
    assert_eq!(line.centroid(), line.segment_point(0.5)); // Middle of the segment
    ```
     */
    pub fn centroid(&self) -> [f64; 2] {
        let x = 0.5 * (self.start()[0] + self.stop()[0]);
        let y = 0.5 * (self.start()[1] + self.stop()[1]);
        return [x, y];
    }

    /**
    Returns the points of a polygon chain which approximates `self`. The number
    of points is defined by the [`SegmentPolygonizer`] (see its docstring). The
    points are regularily distributed over the segment, which means that two
    subsequent points always have the same euclidian distance from each other.

    # Examples

    ```
    use planar_geo::prelude::*;

    let line = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).unwrap();

    let mut iter = line.polygonize(SegmentPolygonizer::NumberSegments(4));

    assert_eq!(iter.next(), Some(line.segment_point(0.0)));
    assert_eq!(iter.next(), Some(line.segment_point(0.25)));
    assert_eq!(iter.next(), Some(line.segment_point(0.5)));
    assert_eq!(iter.next(), Some(line.segment_point(0.75)));
    assert_eq!(iter.next(), Some(line.segment_point(1.0)));
    assert!(iter.next().is_none());
    ```
     */
    pub fn polygonize<'a>(
        &'a self,
        polygonizer: super::SegmentPolygonizer,
    ) -> super::PolygonPointsIterator<'a> {
        let num_segs = match polygonizer {
            super::SegmentPolygonizer::NumberSegments(segs) => segs,
            super::SegmentPolygonizer::MaximumSegmentLength(max_len) => {
                if max_len <= 0.0 {
                    0
                } else {
                    (self.length() / max_len).ceil() as usize
                }
            }
            super::SegmentPolygonizer::MaximumAngle(_) => 0,
        };
        return super::PolygonPointsIterator {
            index: 0,
            num_segs,
            segment: super::SegmentRef::LineSegment(self),
        };
    }

    /**
    Returns an iterator over the start and stop point of `self`

    This is a shorthand for `self.polygonize( Polygonizer::NumberSegments(1))`.

    # Examples

    ```
    use planar_geo::segment::LineSegment;

    let line = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).unwrap();

    let mut iter = line.points();

    assert_eq!(iter.next(), Some(line.segment_point(0.0)));
    assert_eq!(iter.next(), Some(line.segment_point(1.0)));
    assert!(iter.next().is_none());
    ```
     */
    pub fn points<'a>(&'a self) -> super::PolygonPointsIterator<'a> {
        return self.polygonize(super::SegmentPolygonizer::NumberSegments(1));
    }

    /**
    Finds the endpoint of the segments P and Q which is closest to the other segment.  This is
    a reasonable surrogate for the true intersection points in ill-conditioned cases (e.g.
    where two segments are nearly coincident, or where the endpoint of one segment lies almost
    on the other segment).

    This replaces the older CentralEndpoint heuristic, which chose the wrong endpoint in some
    cases where the segments had very distinct slopes and one endpoint lay almost on the other
    segment.
     */
    fn nearest_endpoint(&self, line_segment: &LineSegment) -> [f64; 2] {
        let mut nearest_pt = self.start;
        let mut min_dist = line_segment.euclidian_distance_to_point(self.start);

        let dist = line_segment.euclidian_distance_to_point(self.stop);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = self.stop;
        }
        let dist = self.euclidian_distance_to_point(line_segment.start);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = line_segment.start;
        }
        let dist = self.euclidian_distance_to_point(line_segment.stop);
        if dist < min_dist {
            nearest_pt = line_segment.stop;
        }
        nearest_pt
    }

    fn raw_line_intersection(&self, line_segment: &LineSegment) -> Option<[f64; 2]> {
        let p_min_x = self.start[0].min(self.stop[0]);
        let p_min_y = self.start[1].min(self.stop[1]);
        let p_max_x = self.start[0].max(self.stop[0]);
        let p_max_y = self.start[1].max(self.stop[1]);

        let q_min_x = line_segment.start[0].min(line_segment.stop[0]);
        let q_min_y = line_segment.start[1].min(line_segment.stop[1]);
        let q_max_x = line_segment.start[0].max(line_segment.stop[0]);
        let q_max_y = line_segment.start[1].max(line_segment.stop[1]);

        let int_min_x = p_min_x.max(q_min_x);
        let int_max_x = p_max_x.min(q_max_x);
        let int_min_y = p_min_y.max(q_min_y);
        let int_max_y = p_max_y.min(q_max_y);

        let mid_x = 0.5 * (int_min_x + int_max_x);
        let mid_y = 0.5 * (int_min_y + int_max_y);

        // condition ordinate values by subtracting midpoint
        let p1x = self.start[0] - mid_x;
        let p1y = self.start[1] - mid_y;
        let p2x = self.stop[0] - mid_x;
        let p2y = self.stop[1] - mid_y;
        let q1x = line_segment.start[0] - mid_x;
        let q1y = line_segment.start[1] - mid_y;
        let q2x = line_segment.stop[0] - mid_x;
        let q2y = line_segment.stop[1] - mid_y;

        // unrolled computation using homogeneous coordinates eqn
        let px = p1y - p2y;
        let py = p2x - p1x;
        let pw = p1x * p2y - p2x * p1y;

        let qx = q1y - q2y;
        let qy = q2x - q1x;
        let qw = q1x * q2y - q2x * q1y;

        let xw = py * qw - qy * pw;
        let yw = qx * pw - px * qw;
        let w = px * qy - qx * py;

        let x_int = xw / w;
        let y_int = yw / w;

        // check for parallel lines
        if (x_int.is_nan() || x_int.is_infinite()) || (y_int.is_nan() || y_int.is_infinite()) {
            None
        } else {
            // de-condition intersection point
            Some([x_int + mid_x, y_int + mid_y])
        }
    }

    /**
    This method computes the actual value of the intersection point.
    To obtain the maximum precision from the intersection calculation,
    the coordinates are normalized by subtracting the minimum
    ordinate values (in absolute value).  This has the effect of
    removing common significant digits from the calculation to
    maintain more bits of precision.
     */
    fn proper_intersection(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> [f64; 2] {
        // Computes a segment intersection using homogeneous coordinates.
        // Round-off error can cause the raw computation to fail,
        // (usually due to the segments being approximately parallel).
        // If this happens, a reasonable approximation is computed instead.
        let mut int_pt = self
            .raw_line_intersection(line_segment)
            .unwrap_or_else(|| self.nearest_endpoint(line_segment));

        // NOTE: At this point, JTS does a `Envelope::contains(coord)` check, but confusingly,
        // Envelope::contains(coord) in JTS is actually an *intersection* check, not a true SFS
        // `contains`, because it includes the boundary of the rect.
        if !(BoundingBox::from(self).approx_contains_point(int_pt, epsilon, max_ulps)
            && BoundingBox::from(line_segment).approx_contains_point(int_pt, epsilon, max_ulps))
        {
            // compute a safer result
            // copy the coordinate, since it may be rounded later
            int_pt = self.nearest_endpoint(line_segment);
        }
        int_pt
    }
}

impl Transformation for LineSegment {
    fn translate(&mut self, shift: [f64; 2]) {
        self.start.translate(shift);
        self.stop.translate(shift);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let t = Rotation2::new(angle);

        let pt = [self.start[0] - center[0], self.start[1] - center[1]];
        self.start = t * pt;
        self.start.translate([center[0], center[1]]);

        let pt = [self.stop[0] - center[0], self.stop[1] - center[1]];
        self.stop = t * pt;
        self.stop.translate([center[0], center[1]]);
    }

    fn scale(&mut self, factor: f64) {
        self.start.scale(factor);
        self.stop.scale(factor);
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        // Treatment of special case: vertical line
        if start[0] == stop[0] {
            self.start = [-self.start[0] + 2.0 * start[0], self.start[1]];
            self.stop = [-self.stop[0] + 2.0 * stop[0], self.stop[1]];

        // Treatment of special case: Horizontal line
        } else if start[1] == stop[1] {
            self.start = [self.start[0], -self.start[1] + 2.0 * start[1]];
            self.stop = [self.stop[0], -self.stop[1] + 2.0 * start[1]];

        // All other cases
        } else {
            // Solve the line equation
            let m = (stop[1] - start[1]) / (stop[0] - start[0]);
            let c = (stop[0] * start[1] - start[0] * stop[1]) / (stop[0] - start[0]);

            let d = (self.start[0] + (self.start[1] - c) * m) / (1.0 + m.powi(2));
            self.start = [
                2.0 * d - self.start[0],
                2.0 * d * m - self.start[1] + 2.0 * c,
            ];
            let d = (self.stop[0] + (self.stop[1] - c) * m) / (1.0 + m.powi(2));
            self.stop = [2.0 * d - self.stop[0], 2.0 * d * m - self.stop[1] + 2.0 * c];
        }
    }
}

impl From<&LineSegment> for BoundingBox {
    fn from(value: &LineSegment) -> BoundingBox {
        let [xmin, xmax] = if value.start()[0] > value.stop()[0] {
            [value.stop()[0], value.start()[0]]
        } else {
            [value.start()[0], value.stop()[0]]
        };
        let [ymin, ymax] = if value.start()[1] > value.stop()[1] {
            [value.stop()[1], value.start()[1]]
        } else {
            [value.start()[1], value.stop()[1]]
        };
        return BoundingBox::new(xmin, xmax, ymin, ymax);
    }
}

impl From<&LineSegment> for crate::CentroidData {
    fn from(value: &LineSegment) -> Self {
        // Apply the formula for a triangle with one vertex being (0,0)
        let x = (value.start()[0] + value.stop()[0]) / 3.0;
        let y = (value.start()[1] + value.stop()[1]) / 3.0;
        let area = 0.5 * (value.start()[0] * value.stop()[1] - value.stop()[0] * value.start()[1]);
        return Self { area, x, y };
    }
}

impl std::fmt::Display for LineSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        return write!(f, "{:?} -> {:?}", self.start, self.stop);
    }
}

impl crate::primitive::private::Sealed for LineSegment {}

impl Primitive for LineSegment {
    fn contains_point(&self, p: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        // Quick first check: If point p is outside the segment bounding box, it can't be contained.
        if !BoundingBox::from(self).approx_contains_point(p, epsilon, max_ulps) {
            return false;
        }

        // Check if the points of the segment and p are collinear. If they aren't, p can't be contained.
        // If they are, p is contained, since it is also inside the segment bounding box.
        if epsilon == 0.0 && max_ulps == 0 {
            return geometry_predicates::orient2d(self.start.into(), self.stop.into(), p.into())
                == 0.0;
        } else {
            // This is the solution from https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
            let len_sq = self.length_sq();
            let v = self.start();
            let w = self.stop();
            let t: f64 = ((p[0] - v[0]) * (w[0] - v[0]) + (p[1] - v[1]) * (w[1] - v[1])) / len_sq;
            let t = t.clamp(0.0, 1.0); // This is equal to t = Math.max(0, Math.min(1, t)) in the link

            // Calculate the (squared) distance between the segment and the given point
            let pt = [v[0] + t * (w[0] - v[0]), v[1] + t * (w[1] - v[1])];
            let dist = ((pt[0] - p[0]).powi(2) + (pt[1] - p[1]).powi(2)).sqrt();
            return approx::ulps_eq!(dist, 0.0, epsilon = epsilon, max_ulps = max_ulps);
        }
    }

    fn intersections_line(
        &self,
        line: &crate::line::Line,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        return line.intersections_line_segment(self, epsilon, max_ulps);
    }

    fn intersections_line_segment(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        // This is geo's line_intersection algorithm
        // (https://github.com/georust/geo/blob/3b0d5738f54bd8964f7d1f573bd63dc114587dc4/geo/src/algorithm/line_intersection.rs)
        // translated into the planar_geo crate.

        #[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
        enum Orientation {
            Clockwise,
            CounterClockwise,
            Collinear,
        }

        impl From<f64> for Orientation {
            fn from(value: f64) -> Self {
                if value < 0.0 {
                    return Orientation::Clockwise;
                } else if value > 0.0 {
                    return Orientation::CounterClockwise;
                } else {
                    return Orientation::Collinear;
                }
            }
        }

        // =====================================================================

        // If the segments are identical, they do not intersect by definition
        if std::ptr::eq(self, line_segment) {
            return PrimitiveIntersections::Zero;
        }

        // If the bounding boxes of the line segments don't interact with each other, then the segments cannot intersect at all.
        let bb1 = BoundingBox::from(self);
        let bb2 = BoundingBox::from(line_segment);
        if !bb1.intersects(&bb2) && !bb1.touches(&bb2) {
            return PrimitiveIntersections::Zero;
        }

        let self_start: Orientation = geometry_predicates::orient2d(
            self.start.into(),
            self.stop.into(),
            line_segment.start.into(),
        )
        .into();
        let self_stop: Orientation = geometry_predicates::orient2d(
            self.start.into(),
            self.stop.into(),
            line_segment.stop.into(),
        )
        .into();
        if matches!(
            (self_start, self_stop),
            (Orientation::Clockwise, Orientation::Clockwise)
                | (Orientation::CounterClockwise, Orientation::CounterClockwise)
        ) {
            return PrimitiveIntersections::Zero;
        }

        let line_segment_start: Orientation = geometry_predicates::orient2d(
            line_segment.start.into(),
            line_segment.stop.into(),
            self.start.into(),
        )
        .into();
        let line_segment_stop: Orientation = geometry_predicates::orient2d(
            line_segment.start.into(),
            line_segment.stop.into(),
            self.stop.into(),
        )
        .into();
        if matches!(
            (line_segment_start, line_segment_stop),
            (Orientation::Clockwise, Orientation::Clockwise)
                | (Orientation::CounterClockwise, Orientation::CounterClockwise)
        ) {
            return PrimitiveIntersections::Zero;
        }

        if matches!(
            (self_start, self_stop, line_segment_start, line_segment_stop),
            (
                Orientation::Collinear,
                Orientation::Collinear,
                Orientation::Collinear,
                Orientation::Collinear
            )
        ) {
            let bb1 = BoundingBox::from(self);
            let bb2 = BoundingBox::from(line_segment);

            return match (
                bb1.approx_contains_point(line_segment.start, epsilon, max_ulps),
                bb1.approx_contains_point(line_segment.stop, epsilon, max_ulps),
                bb2.approx_contains_point(self.start, epsilon, max_ulps),
                bb2.approx_contains_point(self.stop, epsilon, max_ulps),
            ) {
                (true, true, _, _) => {
                    PrimitiveIntersections::Two([line_segment.start, line_segment.stop])
                }
                (_, _, true, true) => PrimitiveIntersections::Two([self.start, self.stop]),
                (true, false, true, false) if line_segment.start == self.start => {
                    PrimitiveIntersections::One(line_segment.start)
                }
                (true, _, true, _) => PrimitiveIntersections::Two([line_segment.start, self.start]),
                (true, false, false, true) if line_segment.start == self.stop => {
                    PrimitiveIntersections::One(line_segment.start)
                }
                (true, _, _, true) => PrimitiveIntersections::Two([line_segment.start, self.stop]),
                (false, true, true, false) if line_segment.stop == self.start => {
                    PrimitiveIntersections::One(line_segment.stop)
                }
                (_, true, true, _) => PrimitiveIntersections::Two([line_segment.stop, self.start]),
                (false, true, false, true) if line_segment.stop == self.stop => {
                    PrimitiveIntersections::One(line_segment.stop)
                }
                (_, true, _, true) => PrimitiveIntersections::Two([line_segment.stop, self.stop]),
                _ => PrimitiveIntersections::Zero,
            };
        }

        // At this point we know that there is a single intersection point (since the lines are not
        // collinear).
        //
        // Check if the intersection is an endpoint. If it is, copy the endpoint as the
        // intersection point. Copying the point rather than computing it ensures the point has the
        // exact value, which is important for robustness. It is sufficient to simply check for an
        // endpoint which is on the other line, since at this point we know that the input lines
        // must intersect.
        if self_start == Orientation::Collinear
            || self_stop == Orientation::Collinear
            || line_segment_start == Orientation::Collinear
            || line_segment_stop == Orientation::Collinear
        {
            // Check for two equal endpoints.
            // This is done explicitly rather than by the orientation tests below in order to improve
            // robustness.
            let intersection =
                if self.start == line_segment.start || self.start == line_segment.stop {
                    self.start
                } else if self.stop == line_segment.start || self.stop == line_segment.stop {
                    self.stop
                    // Now check to see if any endpoint lies on the interior of the other segment.
                } else if self_start == Orientation::Collinear {
                    line_segment.start
                } else if self_stop == Orientation::Collinear {
                    line_segment.stop
                } else if line_segment_start == Orientation::Collinear {
                    self.start
                } else {
                    assert_eq!(line_segment_stop, Orientation::Collinear);
                    self.stop
                };
            PrimitiveIntersections::One(intersection)
        } else {
            let intersection = self.proper_intersection(line_segment, epsilon, max_ulps);
            return PrimitiveIntersections::One(intersection);
        }
    }

    fn intersections_arc_segment(
        &self,
        arc_segment: &crate::segment::ArcSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        // Algorithm from https://cp-algorithms.com/geometry/circle-line-intersection.html
        let center = arc_segment.center();

        let x1 = self.start[0] - center[0];
        let y1 = self.start[1] - center[1];
        let x2 = self.stop[0] - center[0];
        let y2 = self.stop[1] - center[1];

        // Calculate values A, B and C as used in the circle-line-intersection algorithm
        // "Koordinatenform", see https://de.wikipedia.org/wiki/Koordinatenform
        let a = y1 - y2;
        let b = x2 - x1;

        // Multiply C by -1 to account for the different definitions of the Wikipedia and the
        // cp-algorithms.com approach
        let c = -(x2 * y1 - x1 * y2);

        let mut intersections = PrimitiveIntersections::Zero;

        for pt in arc_segment
            .intersections_line_circle(a, b, c, epsilon, max_ulps)
            .into_iter()
            .map(From::from)
        {
            if self.contains_point(pt, epsilon, max_ulps) {
                intersections.push(pt);
            }
        }
        return intersections;
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        other.intersections_line_segment(self, epsilon, max_ulps)
    }
}
