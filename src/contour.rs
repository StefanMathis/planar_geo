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

use compare_variables::*;
use num::Integer;
use rayon::iter::ParallelIterator;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::composite::*;
use crate::polysegment::Polysegment;
use crate::primitive::Primitive;
use crate::segment::{Segment, arc_segment::ArcSegment, line_segment::LineSegment};
use crate::{CentroidData, Transformation};
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS};

use bounding_box::BoundingBox;

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
    Returns a slice containing all segments of the underlying [`VecDeque`].

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let contour = Contour::new(polysegment);

    // Now the polysegment has been closed and has three segments
    assert_eq!(contour.segments().len(), contour.len());
    ```
     */
    pub fn segments(&self) -> &[Segment] {
        // Note: When making a contour out of a polysegment, the underlying
        // VecDeque is made contiguous, hence the first slice contains all
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
    pub fn iter(&self) -> std::collections::vec_deque::Iter<'_, Segment> {
        return self.0.iter();
    }

    /**
    Returns a parallel front-to-back iterator over all [`Segment`]s of `self`.
     */
    pub fn par_iter(&self) -> rayon::collections::vec_deque::Iter<'_, Segment> {
        return self.0.par_iter();
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
    Returns whether `self` contains an `other` contour.

    A contour contains another contour if the bounding box of the latter is
    within that of the former and if there is no intersection between the two
    The specified tolerances `epsilon` and `max_ulps` are used within the
    intersection algorithm.

    # Examples
    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ]);
    let outer = Contour::new(polysegment);

    // This contour is inside "outer"
    let polysegment = Polysegment::from_points(&[
        [0.1, 0.1],
        [0.9, 0.1],
        [0.9, 0.9],
        [0.1, 0.9],
    ]);
    let inner = Contour::new(polysegment);
    assert!(outer.contains(&inner, 0.0, 0));

    // This contour intersects with "outer"
    let polysegment = Polysegment::from_points(&[
        [0.1, 0.1],
        [1.1, 0.1],
        [1.1, 0.9],
        [0.1, 0.9],
    ]);
    let inner = Contour::new(polysegment);
    assert!(!outer.contains(&inner, 0.0, 0));
    ```
    */
    pub fn contains(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        // Check if the bounding box of self contains the bounding box of other.
        // If that's not the case, self cannot contain other by definition
        if !BoundingBox::from(self).contains(&BoundingBox::from(other)) {
            return false;
        }

        /*
        Check if all end points of segments of other are within self. If that is
        not the the case, the contour is not contained.
        */
        if other
            .points()
            .any(|point| !self.contains_point(point, epsilon, max_ulps))
        {
            return false;
        }

        /*
        Check for intersection between self and other. This check seems
        superfluous, as the previous point containing check already detects
        some intersections, but is actually necessary because a segment might be
        partially outside self while its end points are inside it. As this is
        a costly check to make, it is done at the end.
         */
        return self
            .polysegment()
            .intersections_polysegment_par(other.polysegment(), epsilon, max_ulps)
            .count()
            == 0;
    }

    // /**
    // TODO
    //  */
    // pub fn interiors_intersect(&self, other: &Self, epsilon: f64, max_ulps: u32)
    // -> bool {     let b_self = BoundingBox::from(self);
    //     let b_other = BoundingBox::from(other);

    //     // If the bounding box do not contain each other or intersect, then the
    //     // interiors cannot intersect as well
    //     if !b_self.contains(&b_other) && !b_other.contains(&b_self) &&
    // !b_self.intersects(&b_other) {         return false;
    //     }

    // }

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
            .iter()
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

    fn intersections_primitive<'a, T: Primitive>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.polysegment()
            .intersections_primitive(primitive, epsilon, max_ulps)
    }

    fn intersections_primitive_par<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.polysegment()
            .intersections_primitive_par(primitive, epsilon, max_ulps)
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

    /**
    This function check if the point is on the contour OR inside it.

    # Examples

    ```
    use planar_geo::prelude::*;

    let rect = Contour::rectangle([0.0, 1.0], [1.0, 0.0]);

    // Point on the polysegment forming the contour
    assert!(rect.contains_point([0.5, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    assert!(rect.polysegment().contains_point([0.5, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // Point inside the contour
    assert!(rect.contains_point([0.5, 0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    assert!(!rect.polysegment().contains_point([0.5, 0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // Point outside contour
    assert!(!rect.contains_point([2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    assert!(!rect.polysegment().contains_point([2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    ```
     */
    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        // First coarse, but fast test: Check if the point is inside the polygon
        let bb = BoundingBox::from(self);
        if !bb.approx_contains_point(point, epsilon, max_ulps) {
            return false;
        }

        // Check if the point is located on the edge of the polysegment
        for segment in self.0.vec_deque().iter() {
            if segment.contains_point(point, epsilon, max_ulps) {
                return true;
            }
        }

        // Use the ray casting algorithm
        if let Ok(ray) = LineSegment::new(
            [bb.xmin() - 1.0, point[1]],
            point,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        ) {
            let mut counter = 0;
            for segment in self.0.vec_deque().iter() {
                counter += segment
                    .intersections_primitive(&ray, epsilon, max_ulps)
                    .len();
            }
            return counter.is_odd();
        } else {
            return false;
        }
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

impl From<&Contour> for BoundingBox {
    fn from(value: &Contour) -> BoundingBox {
        return BoundingBox::from(&value.0);
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
    Returns the arrow head height contained in `self`.

    In case of the [`ArrowHeadSize::Height`] variant, this is simply the
    contained value. For the [`ArrowHeadSize::SideLength`] variant, it is
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
    Returns the arrow head side length contained in `self`.

    In case of the [`ArrowHeadSize::SideLength`] variant, this is simply the
    contained value. For the [`ArrowHeadSize::Height`] variant, it is
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
