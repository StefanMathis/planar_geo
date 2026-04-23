/*!
Defines the [`Contour`] type, a closed [`Composite`] used to describe enclosed
geometric boundaries.

Contours are primarily used as building blocks for
[`Shape`](crate::shape::Shape)s and other higher-level geometric abstractions.
Conceptually, a contour represents a closed outline composed of connected
segments.

Unlike [`Polysegment`], contours are immutable with respect to their segment
structure, ensuring that they always remain closed. This distinction allows
code working with contours to rely on closure as a guaranteed property rather
than a runtime condition.

Most functionality is provided directly by the [`Contour`] type; see its
documentation for details on invariants, construction, and usage.
*/

use std::collections::VecDeque;
use std::f64::consts::TAU;

use compare_variables::*;
use rayon::iter::ParallelIterator;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::composite::*;
use crate::polysegment::Polysegment;
use crate::primitive::Primitive;
use crate::segment::{Segment, SegmentRef, arc_segment::ArcSegment, line_segment::LineSegment};
use crate::{CentroidData, Transformation};
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS};

use bounding_box::{BoundingBox, ToBoundingBox};

/**
A closed [`Polysegment`] whose last segment is guaranteed to connect back to
the first.
*/
#[doc = ""]
#[cfg_attr(feature = "doc-images", doc = "![Contour example][example_contour]")]
#[cfg_attr(
    feature = "doc-images",
    embed_doc_image::embed_doc_image("example_contour", "docs/img/example_contour.svg")
)]
#[cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images not enabled**. Compile docs with
    `cargo doc --features 'doc-images'` and Rust version >= 1.54."
)]
/**

More precisely, the stop point of the last (back) segment is guaranteed to be
equal to the start point of the first (front) segment. This invariant is upheld
by the type itself and cannot be violated through the public API.

Internally, a contour is a newtype wrapper around a [`VecDeque<Segment>`], just
like [`Polysegment`]. However, unlike a polysegment, a contour does not
expose any methods that would allow adding or removing segments. This immutability
ensures that the contour remains closed at all times.

# Intersections

A contour can intersect itself (i.e. at least two of its segments
intersect each other). This can be tested using its [`Composite`] trait
implementation:

```
use planar_geo::prelude::*;

let polysegment = Polysegment::from_points(&[[0.0, 0.0], [2.0, 2.0], [0.0, 2.0], [2.0, 0.0]]);
let contour = Contour::new(polysegment);
assert_eq!(contour.intersections_contour(&contour, 0.0, 0).count(), 1);
```

One difference to a [`Polysegment`] is that the "touching" of the first and
the last segment is not counted as an intersection:

```
use planar_geo::prelude::*;

let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]);

// Last and first segment are "touching" at 0.0
assert_eq!(polysegment.intersections_polysegment(&polysegment, 0.0, 0).count(), 1);

let contour = Contour::new(polysegment);
assert_eq!(contour.intersections_contour(&contour, 0.0, 0).count(), 0);
```

# Constructing a contour

The standard constructor is [`new`](Contour::new) which takes a [`Polysegment`],
closes it ([`Polysegment::close`]) and makes the underlying deque contiguous
via [`VecDeque::make_contiguous`]. As a result, all segments of a contour are
always stored contiguously in memory and can be accessed as a single slice using
[`segments`](Contour::segments). A [`Polysegment`] can also be converted using
a [`From`] implementation (which uses [`new`](Contour::new)). Similarily, a
[`Contour`] can be converted back to a [`Polysegment`] via [`From`]. It is not
possible to add or remove segments to / from a [`Contour`]; it needs to be
converted into a [`Polysegment`], modified and then converted back.

For common geometric bodies, a variety of convenience constructors is available
- [`circle`](Contour::circle): Constructs a circle contour.
- [`rectangle`](Contour::rectangle): Constructs a rectangle contour.
- [`arrow_from_tail_length_angle`](Contour::arrow_from_tail_length_angle):
Constructs an arrow from the defined values.
- [`arrow_from_tail_head`](Contour::arrow_from_tail_head):
Constructs an arrow from the defined values.
- [`arrow_from_head_length_angle`](Contour::arrow_from_head_length_angle):
Constructs an arrow from the defined values.

# Access of individual segments

Extracting individual elements works in the same fashion as it does for a
[`Polysegment`] via either the [`get`](Contour::get) method or indexing. Of
course, the slice returned by [`segments`](Contour::segments) can also be used
and returns the same segment for an index.

# Serialization and deserialization

When the `serde` feature is enabled, a polysegment can be serialized and
deserialized using the [`serde`] crate. It uses the same serialized
representation as a [`VecDeque`].
 */
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Contour(Polysegment);

impl Contour {
    /**
    Creates a contour from a [`Polysegment`].

    This calls the [`From`] implementation, which performs the following
    operations on the [`Polysegment`]:
    1) [`Polysegment::close`]
    2) [`Polysegment::make_contiguous`]

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);

    // Two segments
    assert_eq!(polysegment.len(), 2);

    let contour = Contour::new(polysegment);

    // Now the polysegment has been closed and has three segments
    assert_eq!(contour.polysegment().len(), 3);
    ```
     */
    pub fn new(polysegment: Polysegment) -> Self {
        return polysegment.into();
    }

    /**
    Returns a reference to the underlying [`Polysegment`].
     */
    pub fn polysegment(&self) -> &Polysegment {
        return &self.0;
    }

    /**
    Returns a reference to the underlying [`VecDeque`].

    This allows using all methods of [`VecDeque`] which use a shared reference,
    e.g. its iterators, accessor methods etc.
     */
    pub fn vec_deque(&self) -> &VecDeque<Segment> {
        return self.0.vec_deque();
    }

    /**
    Returns a slice covering all segments of the underlying [`VecDeque`].

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let contour = Contour::new(polysegment);

    // Now the polysegment has been closed and has three segments
    assert_eq!(contour.as_slice().len(), contour.len());
    ```
     */
    pub fn as_slice(&self) -> &[Segment] {
        // Note: When making a contour out of a polysegment, the underlying
        // VecDeque is made contiguous, hence the first slice covers all
        // segments and the second one is empty
        return self.0.as_slices().0;
    }

    /**
    Provides a reference to the back element, or [`None`] if the contour is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let contour = Contour::new(Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]));
    let ls: Segment = LineSegment::new([0.0, 1.0], [0.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(contour.back(), Some(&ls));
    ```
     */
    pub fn back(&self) -> Option<&Segment> {
        return self.0.back();
    }

    /**
    Provides a reference to the front element, or [`None`] if the contour is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let contour = Contour::new(Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]));
    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(contour.front(), Some(&ls));
     */
    pub fn front(&self) -> Option<&Segment> {
        return self.0.front();
    }

    /**
    Returns whether the contour is empty or not.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Empty contour
    assert!(Contour::new(Polysegment::new()).is_empty());

    // Non-empty contour
    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let contour = Contour::new(polysegment);
    assert!(!contour.is_empty());
    ```
     */
    pub fn is_empty(&self) -> bool {
        return self.0.is_empty();
    }

    /**
    Provides a reference to the [`Segment`] at the given index.

    The segment at index 0 is the front of the contour.

    # Examples

    ```
     use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let contour = Contour::new(polysegment);
    assert!(contour.get(0).is_some());
    assert!(contour.get(3).is_none());
    ```
     */
    pub fn get(&self, index: usize) -> Option<&Segment> {
        return self.0.get(index);
    }

    /**
    Reverses all [`Segment`]s of the polysegment by switching their start and stop
    points (see [`Segment::reverse`]). To ensure the "connected" property holds
    true, the ordering of the segments itself is exchanged as well.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let mut contour = Contour::new(polysegment);
    assert_eq!(contour[1].start(), [1.0, 0.0]);
    assert_eq!(contour[1].stop(), [0.0, 1.0]);

    contour.reverse();
    assert_eq!(contour[1].start(), [0.0, 1.0]);
    assert_eq!(contour[1].stop(), [1.0, 0.0]);
    ```
     */
    pub fn reverse(&mut self) -> () {
        self.0.reverse();
    }

    /**
    Returns the number of [`Segment`]s in `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let contour = Contour::new(Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]));
    assert_eq!(contour.len(), 3);
    ```
     */
    pub fn len(&self) -> usize {
        return self.0.len();
    }

    /**
    Returns a front-to-back iterator over all [`Segment`]s of `self`.
     */
    pub fn segments(&self) -> std::collections::vec_deque::Iter<'_, Segment> {
        return self.0.segments();
    }

    /**
    Returns a parallel front-to-back iterator over all [`Segment`]s of `self`.
     */
    pub fn segments_par(&self) -> rayon::collections::vec_deque::Iter<'_, Segment> {
        return self.0.segments_par();
    }

    /**
    Returns an iterator over all points of the contour.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let contour = Contour::new(polysegment);

    let mut iter = contour.points();
    assert_eq!(iter.next(), Some([0.0, 0.0]));
    assert_eq!(iter.next(), Some([1.0, 0.0]));
    assert_eq!(iter.next(), Some([1.0, 1.0]));
    assert_eq!(iter.next(), None);
    ```
     */
    pub fn points(&self) -> PointIterator<'_> {
        return PointIterator::new(self.polysegment(), true, Polygonizer::default());
    }

    /**
    Returns the points of a polygon polysegment which approximates `self`. The
    individual segments are "polygonized" via [`Segment::polygonize`] and
    an [`SegmentPolygonizer`](crate::segment::SegmentPolygonizer) specified
    within [`Polygonizer`]. See the docstring of the latter fore more.
     */
    pub fn polygonize(&self, polygonizer: Polygonizer) -> PointIterator<'_> {
        return PointIterator::new(self.polysegment(), true, polygonizer);
    }

    /**
    TODO
     */
    pub fn overlaps_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        let segment: SegmentRef = segment.into();

        // If the segment is outside the bounding box of self, then the segment
        // and the contour surely don't overlap
        if !self
            .bounding_box()
            .approx_covers(&segment.bounding_box(), epsilon, max_ulps)
        {
            return false;
        }

        match segment {
            SegmentRef::LineSegment(line_segment) => {
                /*
                Sort the intersection points by their distance to
                line_segment.start(). Each pair of neighboring intersections
                describes a part of line_segment. If the middle points of these
                parts are all contained in self, then the entire line segment is
                also contained in self.
                 */
                let mut intersections: Vec<[f64; 2]> = self
                    .intersections_primitive_par(&segment, epsilon, max_ulps)
                    .map(|i| i.point)
                    .collect();
                intersections.push(line_segment.start());
                intersections.push(line_segment.stop());

                // Sort the intersections in ascending order by their distance
                // to line_segment.start().
                let s = line_segment.start();
                intersections.sort_unstable_by(|a, b| {
                    ((a[0] - s[0]).powi(2) + (a[1] - s[1]).powi(2))
                        .total_cmp(&((b[0] - s[0]).powi(2) + (b[1] - s[1]).powi(2)))
                });

                // Create individual segments out of the sorted intersections
                // and check if their middle point is contained in self
                for w in intersections.windows(2) {
                    if let Some(start) = w.get(0) {
                        if let Some(stop) = w.get(1) {
                            let mx = 0.5 * (start[0] + stop[0]);
                            let my = 0.5 * (start[1] + stop[1]);
                            if self.contains_point([mx, my], epsilon, max_ulps) {
                                return true;
                            }
                        }
                    }
                }
                return false;
            }
            SegmentRef::ArcSegment(arc_segment) => {
                /*
                Calculate the angles of all intersection points relative to the
                center of arc_segment, then sort them. This separates
                arc_segment in multiple arc segments. Calculate their middle
                point and check if that point is contained within self. If that is
                true for all partial segments, then the entire arc_segment is
                also contained in self.
                 */
                let c = arc_segment.center();
                let r = arc_segment.radius();

                let mut offset_angles: Vec<f64> = self
                    .intersections_primitive_par(&segment, epsilon, max_ulps)
                    .map(|i| ((i.point[1] - c[1]).atan2(i.point[0] - c[0])).rem_euclid(TAU))
                    .collect();

                offset_angles.sort_unstable_by(f64::total_cmp);

                for w in offset_angles.windows(2) {
                    if let Some(start) = w.get(0) {
                        if let Some(stop) = w.get(1) {
                            let mid_angle = 0.5 * (stop + start);
                            let mx = c[0] + r * mid_angle.cos();
                            let my = c[1] + r * mid_angle.sin();
                            if self.contains_point([mx, my], epsilon, max_ulps) {
                                return true;
                            }
                        }
                    }
                }
                return false;
            }
        }
    }

    /**
    TODO
     */
    pub fn overlaps_contour<'a>(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        let b_self = BoundingBox::from(self);
        let b_other = BoundingBox::from(other);

        // If the bounding boxes do not overlap, then self and other also can't overlap
        if !b_self.overlaps(&b_other) {
            return false;
        }

        return self
            .intersections_contour_par(other, epsilon, max_ulps)
            .any(|i| {
                if let Some(s_self) = self.segment(i.left) {
                    return other.overlaps_segment(s_self, epsilon, max_ulps);
                }
                if let Some(s_other) = other.segment(i.right) {
                    return self.overlaps_segment(s_other, epsilon, max_ulps);
                }
                return false;
            });
    }

    /**
    Cuts `self` into multiple polysegments by intersecting it with `other` and returns
    those them.

    The specified tolerances `epsilon` and `max_ulps` are used to determine
    intersections, see documentation of
    [`PrimitiveIntersections`](crate::primitive::PrimitiveIntersections).

    # Examples

    ```
    use planar_geo::prelude::*;

    let points = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let contour: Contour = Polysegment::from_points(points).into();
    let cut = Polysegment::from_points(&[[-1.0, 1.0], [3.0, 1.0]]);

    let separated_lines = contour.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    // This cut results in two separate polysegments
    assert_eq!(separated_lines.len(), 2);
    ```
     */
    pub fn intersection_cut(
        &self,
        other: &Polysegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Polysegment> {
        let mut lines = self
            .polysegment()
            .intersection_cut(other, epsilon, max_ulps);
        if lines.len() > 1 {
            match lines.pop() {
                Some(mut last_polyline) => {
                    let mut first_polyline =
                        std::mem::replace(lines.first_mut().unwrap(), Polysegment::new());

                    last_polyline.append(&mut first_polyline);
                    lines[0] = last_polyline;
                }
                None => (),
            }
        }

        return lines;
    }

    /**
    Calculates the area enclosed by `self`.

    This is an exact calculation which does not approximate arc segments as
    polylines. The algorithm works as follows:

    For each segment, calculate the individual area of a simple shape consisting
    of the segment and the origin.

    ## LineSegment
    Connect start and stop to the origin, then calculate the shape area as
    `0.5 * ((stop[0] - start[0]) * (origin[1] - start[1]) - (origin[0] - start[0]) * (stop[1] - start[1]))`.

    ## ArcSegment
    Separate the arc into the following three shapes:
    1) triangle start -> stop -> origin
    2) circular segment start -> stop -> arc center
    3) triangle start -> stop -> center.
    Calculate the areas of those segments, then sum them up: `A = A1 + A2 - A3`.

    All areas are then summed up, taking their sign into account.

    # Examples

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::PI;

    // Circle with origin at [10, -20] and a radius of 2
    let polygon = Contour::circle([10.0, -20.0], 2.0);
    assert_eq!(polygon.area(), PI*2.0_f64.powi(2));
    ```
    */
    pub fn area(&self) -> f64 {
        return self
            .segments()
            .map(|segment| match segment {
                Segment::LineSegment(line) => {
                    0.5 * (line.start()[0] * line.stop()[1] - line.stop()[0] * line.start()[1])
                }
                Segment::ArcSegment(arc) => {
                    let center = arc.center();

                    // A naive implementation would be:
                    // let area_1 =
                    //     0.5 * (arc.start()[0] * arc.stop()[1] - arc.stop()[0] * arc.start()[1]);
                    // let area_3 = 0.5
                    //     * ((arc.stop()[0] - arc.start()[0]) * (center[1] - arc.start()[1])
                    //         - (center[0] - arc.start()[0]) * (arc.stop()[1] - arc.start()[1]));
                    //
                    // area_1 and area_3 can be summarized further as:
                    let area_13 = 0.5
                        * (center[0] * (arc.stop()[1] - arc.start()[1])
                            - (arc.stop()[0] - arc.start()[0]) * center[1]);

                    let angle = arc.offset_angle();
                    let radius = arc.radius();
                    let area_2 = 0.5 * angle * radius.powi(2);

                    return area_13 + area_2;
                }
            })
            .sum::<f64>()
            .abs();
    }

    /**
    Returns the combined length of all segments of `self`.

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::{PI, TAU, SQRT_2};

    // Square with a side length of 1
    let points = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let polysegment = Polysegment::from_points(points);
    let contour = Contour::new(polysegment);

    assert_eq!(contour.length(), 4.0);
    ```
     */
    pub fn length(&self) -> f64 {
        return self.polysegment().length();
    }

    /**
    Creates a contour representing a circle.

    This is a shorthand for constructing an [`ArcSegment`], turning it into a
    [`Polysegment`], and converting the polysegment into a contour.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let contour = Contour::circle([10.0, -10.0], 2.0);
    assert_eq!(contour.area(), PI * 2.0f64.powi(2));
    ```
     */
    pub fn circle(center: [f64; 2], radius: f64) -> Self {
        match ArcSegment::circle(center, radius) {
            Ok(segment) => {
                let polysegment = Polysegment::from(Segment::from(segment));
                return Contour::new(polysegment);
            }
            Err(_) => return Contour::new(Polysegment::new()),
        }
    }

    /**
    Creates a contour representing a rectangle from two points defining the
    minimum and maximum x- and y-values

    This is a shorthand for constructing a [`Polysegment`] as a rectangular box
    of [`LineSegment`]s, and converting the polysegment into a contour.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polygon = Contour::rectangle([0.0, 1.0], [1.0, 0.0]);
    let verts: Vec<[f64; 2]> = polygon.points().collect();
    assert_eq!(verts.len(), 4);
    assert_eq!(verts[0], [0.0, 1.0]);
    assert_eq!(verts[1], [0.0, 0.0]);
    assert_eq!(verts[2], [1.0, 0.0]);
    assert_eq!(verts[3], [1.0, 1.0]);
    ```
     */
    pub fn rectangle(pt1: [f64; 2], pt2: [f64; 2]) -> Self {
        let xmin = pt1[0].min(pt2[0]);
        let xmax = pt1[0].max(pt2[0]);
        let ymin = pt1[1].min(pt2[1]);
        let ymax = pt1[1].max(pt2[1]);
        let lower_left = [xmin, ymin];
        let upper_left = [xmin, ymax];
        let lower_right = [xmax, ymin];
        let upper_right = [xmax, ymax];
        return Polysegment::from_points(&[upper_left, lower_left, lower_right, upper_right])
            .into();
    }

    /**
    Creates a contour forming an arrow from the `tail`, the arrow `length`, the
    `stem_width` and the arrow `angle` (in radians). The `head_size` is
    specified with the [`ArrowHeadSize`] enum (either as head height or side
    length).

    # Examples

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::PI;

    let tail = [0.0, 0.0];
    let length = 2.0;
    let angle = 0.0;
    let stem_width = 0.0;
    let arrow = Contour::arrow_from_tail_length_angle(tail, length, angle, stem_width, ArrowHeadSize::Height(1.0)).unwrap();
    let verts: Vec<[f64; 2]> = arrow.points().collect();

    assert_eq!(verts.len(), 6);
    approx::assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[1], [1.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[2], [1.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[3], [2.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[4], [1.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[5], [1.0, 0.0], epsilon = 1e-5);
    ```
     */
    pub fn arrow_from_tail_length_angle(
        tail: [f64; 2],
        length: f64,
        angle: f64,
        stem_width: f64,
        head_size: ArrowHeadSize,
    ) -> crate::error::Result<Self> {
        let height = head_size.height();

        // Sanity checks
        compare_variables!(0.0 <= length)?;
        compare_variables!(0.0 <= height)?;
        compare_variables!(0.0 <= stem_width)?;

        let sin = angle.sin();
        let cos = angle.cos();

        // Calculate the points defining the path
        let triangle_side_length = head_size.side_length();

        let x_flat_end_center = tail[0] + (length - height) * cos;
        let y_flat_end_center = tail[1] + (length - height) * sin;

        let x_tail_l = tail[0] - 0.5 * stem_width * sin;
        let y_tail_l = tail[1] + 0.5 * stem_width * cos;
        let x_tail_to_head_l = x_tail_l + (length - height) * cos;
        let y_tail_to_head_l = y_tail_l + (length - height) * sin;
        let x_head_outer_l = x_flat_end_center - 0.5 * triangle_side_length * sin;
        let y_head_outer_l = y_flat_end_center + 0.5 * triangle_side_length * cos;

        let x_pointed_end = tail[0] + length * cos;
        let y_pointed_end = tail[1] + length * sin;

        let x_tail_r = tail[0] + 0.5 * stem_width * sin;
        let y_tail_r = tail[1] - 0.5 * stem_width * cos;
        let x_tail_to_head_r = x_tail_r + (length - height) * cos;
        let y_tail_to_head_r = y_tail_r + (length - height) * sin;
        let x_head_outer_r = x_flat_end_center + 0.5 * triangle_side_length * sin;
        let y_head_outer_r = y_flat_end_center - 0.5 * triangle_side_length * cos;

        let mut vertices: Vec<[f64; 2]> = Vec::with_capacity(7);
        vertices.push([x_tail_l, y_tail_l]);
        vertices.push([x_tail_to_head_l, y_tail_to_head_l]);
        vertices.push([x_head_outer_l, y_head_outer_l]);
        vertices.push([x_pointed_end, y_pointed_end]);
        vertices.push([x_head_outer_r, y_head_outer_r]);
        vertices.push([x_tail_to_head_r, y_tail_to_head_r]);
        vertices.push([x_tail_r, y_tail_r]);

        return Ok(Polysegment::from_points(&vertices).into());
    }

    /**
    Creates an arrow from `head` and `tail`, otherwise identical to
    [`Contour::arrow_from_tail_length_angle`].

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::PI;

    let tail = [0.0, 0.0];
    let head = [2.0, 0.0];
    let stem_width = 0.0;
    let arrow = Contour::arrow_from_tail_head(tail, head, stem_width, ArrowHeadSize::Height(1.0)).unwrap();
    let verts: Vec<[f64; 2]> = arrow.points().collect();

    approx::assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[1], [1.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[2], [1.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[3], [2.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[4], [1.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[5], [1.0, 0.0], epsilon = 1e-5);
    ```
     */
    pub fn arrow_from_tail_head(
        tail: [f64; 2],
        head: [f64; 2],
        stem_width: f64,
        head_size: ArrowHeadSize,
    ) -> crate::error::Result<Self> {
        let length = ((head[0] - tail[0]).powi(2) + (head[1] - tail[1]).powi(2)).sqrt();
        let angle = (head[1] - tail[1]).atan2(head[0] - tail[0]);
        return Contour::arrow_from_tail_length_angle(tail, length, angle, stem_width, head_size);
    }

    /**
    Creates an arrow from `head` and `tail`, otherwise identical to
    [`Contour::arrow_from_tail_length_angle`].

    # Examples

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::PI;

    let tail = [0.0, 0.0];
    let length = 2.0;
    let angle = 0.0;
    let stem_width = 0.0;
    let arrow = Contour::arrow_from_head_length_angle(tail, length, angle, stem_width, ArrowHeadSize::Height(1.0)).unwrap();
    let verts: Vec<[f64; 2]> = arrow.points().collect();

    assert_eq!(verts.len(), 6);
    approx::assert_abs_diff_eq!(verts[0], [-2.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[1], [-1.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[2], [-1.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[3], [0.0, 0.0], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[4], [-1.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    approx::assert_abs_diff_eq!(verts[5], [-1.0, 0.0], epsilon = 1e-5);
    ```
     */
    pub fn arrow_from_head_length_angle(
        head: [f64; 2],
        length: f64,
        angle: f64,
        stem_width: f64,
        head_size: ArrowHeadSize,
    ) -> crate::error::Result<Self> {
        let x_tail = head[0] - length * angle.cos();
        let y_tail = head[1] - length * angle.sin();
        let tail = [x_tail, y_tail];
        return Contour::arrow_from_tail_length_angle(tail, length, angle, stem_width, head_size);
    }

    /// The algorithms are almost identical, hence the function is reused and
    /// a compile-time boolean is used to differentiate
    fn covers_or_contains_point<const COVERS: bool>(
        &self,
        point: [f64; 2],
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        // First coarse, but fast test: Check if the point is inside the bounding box
        let bb = BoundingBox::from(self);
        if COVERS {
            if !bb.approx_covers_point(point, epsilon, max_ulps) {
                return false;
            }
        } else {
            if !bb.contains_point(point) {
                return false;
            }
        }

        // Check if the point is located on the edge of the polysegment
        for segment in self.segments() {
            if segment.covers_point(point, epsilon, max_ulps) {
                if COVERS {
                    return true;
                } else {
                    return false;
                }
            }
        }

        let mut inside = false;
        for segment in self.segments() {
            match segment {
                Segment::LineSegment(l) => {
                    let [x1, y1] = l.start();
                    let [x2, y2] = l.stop();

                    // Skip horizontal line segments
                    if approx::ulps_eq!(y1, y2, epsilon = epsilon, max_ulps = max_ulps) {
                        continue;
                    }

                    if (y1 > point[1]) != (y2 > point[1]) {
                        let x_intersect = x1 + (point[1] - y1) * (x2 - x1) / (y2 - y1);
                        if x_intersect > point[0] {
                            inside = !inside;
                        }
                    }
                }
                Segment::ArcSegment(a) => {
                    for split in a.split_y_monotonic(epsilon, max_ulps) {
                        if let Some(partial_arc) = split {
                            let [cx, cy] = partial_arc.center();
                            let r = partial_arc.radius();

                            let dy = point[1] - cy;
                            let disc = r * r - dy * dy;

                            if disc < 0.0 {
                                continue;
                            }
                            let sqrt = disc.sqrt();

                            // Only one of these can be valid due to y-monotonicity
                            let candidates = [cx - sqrt, cx + sqrt];

                            let mut hit = None;

                            for x in candidates {
                                let theta = (point[1] - cy).atan2(x - cx);
                                if partial_arc.covers_angle(theta) {
                                    hit = Some(x);
                                    break;
                                }
                            }

                            if let Some(x_intersect) = hit {
                                if x_intersect > point[0] {
                                    inside = !inside;
                                }
                            }
                        } else {
                            // Once a None is hit, no other arc will follow
                            break;
                        }
                    }
                }
            }
        }
        return inside;
    }
}

impl crate::composite::private::Sealed for Contour {}

impl Composite for Contour {
    fn segment(&self, key: SegmentKey) -> Option<&crate::segment::Segment> {
        return self.polysegment().segment(key);
    }

    fn num_segments(&self) -> usize {
        return self.polysegment().num_segments();
    }

    fn centroid(&self) -> [f64; 2] {
        return CentroidData::from(self).into();
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self.0.iter();
    }

    fn par_iter<'a>(
        &'a self,
    ) -> impl ParallelIterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self.0.par_iter();
    }

    fn intersections_polysegment<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.polysegment()
            .intersections_polysegment(polysegment, epsilon, max_ulps)
    }

    fn intersections_polysegment_par<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.polysegment()
            .intersections_polysegment_par(polysegment, epsilon, max_ulps)
    }

    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a {
        let is_identical = std::ptr::eq(self, contour);
        self.polysegment()
            .intersections_polysegment(contour.polysegment(), epsilon, max_ulps)
            .filter(move |i| {
                if is_identical
                    && i.left.segment_idx == 0
                    && i.right.segment_idx + 1 == self.num_segments()
                {
                    if let Some(start) = self.polysegment().front() {
                        let pt = start.start();
                        return !approx::ulps_eq!(
                            pt,
                            i.point,
                            epsilon = epsilon,
                            max_ulps = max_ulps
                        );
                    }
                };
                return true;
            })
    }

    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        let is_identical = std::ptr::eq(self, contour);
        self.polysegment()
            .intersections_polysegment_par(contour.polysegment(), epsilon, max_ulps)
            .filter(move |i| {
                if is_identical
                    && i.left.segment_idx == 0
                    && i.right.segment_idx + 1 == self.num_segments()
                {
                    if let Some(start) = self.polysegment().front() {
                        let pt = start.start();
                        return !approx::ulps_eq!(
                            pt,
                            i.point,
                            epsilon = epsilon,
                            max_ulps = max_ulps
                        );
                    }
                };
                return true;
            })
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a {
        shape
            .intersections_contour(self, epsilon, max_ulps)
            .map(Intersection::switch)
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        shape
            .intersections_contour_par(self, epsilon, max_ulps)
            .map(Intersection::switch)
    }

    fn intersections_composite<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a
    where
        Self: Sized,
    {
        return other
            .intersections_contour(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn intersections_composite_par<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        Self: Sized,
    {
        return other
            .intersections_contour_par(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        return self.covers_or_contains_point::<true>(point, epsilon, max_ulps);
    }

    fn covers_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        let segment: SegmentRef = segment.into();
        // If the segment is outside the bounding box of self, then it is surely
        // not covered.
        if !self
            .bounding_box()
            .approx_covers(&segment.bounding_box(), epsilon, max_ulps)
        {
            return false;
        }

        // Are the start or stop point outside self?
        if !self.covers_point(segment.start(), epsilon, max_ulps)
            || !self.covers_point(segment.stop(), epsilon, max_ulps)
        {
            return false;
        }

        match segment {
            SegmentRef::LineSegment(line_segment) => {
                /*
                Sort the intersection points by their distance to
                line_segment.start(). Each pair of neighboring intersections
                describes a part of line_segment. If the middle points of these
                parts are all covered by self, then the entire line segment is
                also covered by self.
                 */
                let mut intersections: Vec<[f64; 2]> = self
                    .intersections_primitive_par(&segment, epsilon, max_ulps)
                    .map(|i| i.point)
                    .collect();

                // Sort the intersections in ascending order by their distance
                // to line_segment.start().
                let s = line_segment.start();
                intersections.sort_unstable_by(|a, b| {
                    ((a[0] - s[0]).powi(2) + (a[1] - s[1]).powi(2))
                        .total_cmp(&((b[0] - s[0]).powi(2) + (b[1] - s[1]).powi(2)))
                });

                // Create individual segments out of the sorted intersections
                // and check if their middle point is outside self.
                for w in intersections.windows(2) {
                    if let Some(start) = w.get(0) {
                        if let Some(stop) = w.get(1) {
                            let mx = 0.5 * (start[0] + stop[0]);
                            let my = 0.5 * (start[1] + stop[1]);
                            if !self.covers_point([mx, my], epsilon, max_ulps) {
                                return false;
                            }
                        }
                    }
                }
                return true;
            }
            SegmentRef::ArcSegment(arc_segment) => {
                /*
                Calculate the angles of all intersection points relative to the
                center of arc_segment, then sort them. This separates
                arc_segment in multiple arc segments. Calculate their middle
                point and check if that point is covered by self. If that is
                true for all partial segments, then the entire arc_segment is
                also covered by self.
                 */
                let c = arc_segment.center();
                let r = arc_segment.radius();

                let mut offset_angles: Vec<f64> = self
                    .intersections_primitive_par(&segment, epsilon, max_ulps)
                    .map(|i| ((i.point[1] - c[1]).atan2(i.point[0] - c[0])).rem_euclid(TAU))
                    .collect();

                offset_angles.sort_unstable_by(f64::total_cmp);

                for w in offset_angles.windows(2) {
                    if let Some(start) = w.get(0) {
                        if let Some(stop) = w.get(1) {
                            let mid_angle = 0.5 * (stop + start);
                            let mx = c[0] + r * mid_angle.cos();
                            let my = c[1] + r * mid_angle.sin();
                            if !self.covers_point([mx, my], epsilon, max_ulps) {
                                return false;
                            }
                        }
                    }
                }
                return true;
            }
        }
    }

    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        return self.covers_or_contains_point::<false>(point, epsilon, max_ulps);
    }

    fn contains_segment<'a, T: Into<crate::prelude::SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        let segment: SegmentRef = segment.into();

        // If the segment is outside the bounding box of self, then it is surely
        // not covered.
        if !self.bounding_box().contains(&segment.bounding_box()) {
            return false;
        }

        // Are the start or stop point outside self?
        if !self.contains_point(segment.start(), epsilon, max_ulps)
            || !self.contains_point(segment.stop(), epsilon, max_ulps)
        {
            return false;
        }

        return self
            .intersections_primitive(&segment, epsilon, max_ulps)
            .map(|i| i.point)
            .count()
            == 0;
    }

    fn covers_composite<'a, T: Composite + Sync>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        return other.covers_contour(self, epsilon, max_ulps);
    }

    fn contains_composite<'a, T: Composite + Sync>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        return other.contains_contour(self, epsilon, max_ulps);
    }
}

impl std::ops::Index<usize> for Contour {
    type Output = Segment;

    fn index(&self, index: usize) -> &Self::Output {
        &self.polysegment()[index]
    }
}

impl Transformation for Contour {
    fn translate(&mut self, shift: [f64; 2]) -> () {
        return self.0.translate(shift);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) -> () {
        return self.0.rotate(center, angle);
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        return self.0.line_reflection(start, stop);
    }

    fn point_reflection(&mut self, point: [f64; 2]) -> () {
        return self.0.point_reflection(point);
    }

    fn scale(&mut self, factor: f64) -> () {
        return self.0.scale(factor);
    }
}

impl ToBoundingBox for Contour {
    fn bounding_box(&self) -> BoundingBox {
        return self.0.bounding_box();
    }
}

impl From<&Contour> for CentroidData {
    fn from(value: &Contour) -> Self {
        return (&value.0).into();
    }
}

impl From<Polysegment> for Contour {
    fn from(mut polysegment: Polysegment) -> Self {
        if let Some(first) = polysegment.vec_deque().front() {
            if let Some(last) = polysegment.vec_deque().back() {
                // If start and stop are identical, the line_segment_from_points
                // returns None => the polysegment is already closed!
                if let Ok(ls) = LineSegment::new(
                    last.stop(),
                    first.start(),
                    DEFAULT_EPSILON,
                    DEFAULT_MAX_ULPS,
                ) {
                    polysegment.push_back(ls.into());
                }
            }
        }
        // This is necessary to make Contour::segments work.
        polysegment.make_contiguous();
        return Contour(polysegment);
    }
}

impl From<Contour> for Polysegment {
    fn from(polygon: Contour) -> Self {
        return polygon.0;
    }
}

impl From<BoundingBox> for Contour {
    fn from(bounding_box: BoundingBox) -> Self {
        return Polysegment::from(bounding_box).into();
    }
}

impl From<&BoundingBox> for Contour {
    fn from(bounding_box: &BoundingBox) -> Self {
        return Polysegment::from(bounding_box).into();
    }
}

impl From<Segment> for Contour {
    fn from(value: Segment) -> Self {
        return Polysegment::from(value).into();
    }
}

impl From<LineSegment> for Contour {
    fn from(value: LineSegment) -> Self {
        return Polysegment::from(value).into();
    }
}

impl From<ArcSegment> for Contour {
    fn from(value: ArcSegment) -> Self {
        return Polysegment::from(value).into();
    }
}

/**
Enum to define the size of an arrow head either by its height or its side length
(assuming a symmetric triangle with equal side lengths). This enum is mainly
used in the "arrow" constructors of [`Contour`].
 */
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArrowHeadSize {
    /// Height of the triangle forming the arrow length.
    Height(f64),
    /// Side length of the symmetric triangle forming the arrow head.
    SideLength(f64),
}

impl ArrowHeadSize {
    /**
    Returns the arrow head height covered in `self`.

    In case of the [`ArrowHeadSize::Height`] variant, this is simply the
    covered value. For the [`ArrowHeadSize::SideLength`] variant, it is
    calculated as `sqrt(3) / 2` (half the square root of 3).

    # Examples

    ```
    use planar_geo::contour::ArrowHeadSize;

    let arrow_head_size = ArrowHeadSize::Height(5.0);
    assert_eq!(arrow_head_size.height(), 5.0);

    let arrow_head_size = ArrowHeadSize::SideLength(2.0);
    assert_eq!(arrow_head_size.height(), 2.0 * 0.5 * 3.0f64.sqrt());
    ```
     */
    pub fn height(&self) -> f64 {
        match self {
            ArrowHeadSize::Height(h) => return *h,
            ArrowHeadSize::SideLength(len) => return 0.5 * 3.0f64.sqrt() * *len,
        }
    }

    /**
    Returns the arrow head side length covered in `self`.

    In case of the [`ArrowHeadSize::SideLength`] variant, this is simply the
    covered value. For the [`ArrowHeadSize::Height`] variant, it is
    calculated as `2 / sqrt(3)` (2 divided by square root of 3).

    # Examples

    ```
    use planar_geo::contour::ArrowHeadSize;

    let arrow_head_size = ArrowHeadSize::Height(5.0);
    assert_eq!(arrow_head_size.side_length(), 5.0 * 2.0 / 3.0f64.sqrt());

    let arrow_head_size = ArrowHeadSize::SideLength(2.0);
    assert_eq!(arrow_head_size.side_length(), 2.0);
    ```
     */
    pub fn side_length(&self) -> f64 {
        match self {
            ArrowHeadSize::Height(h) => return 2.0 / 3.0f64.sqrt() * *h,
            ArrowHeadSize::SideLength(len) => return *len,
        }
    }
}
