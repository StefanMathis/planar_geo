/*!
Defines the [`Polysegment`] type, a foundational [`Composite`] used throughout
this crate.

Segment polysegments serve as the basis for higher-level geometric types such as
[`Contour`](crate::contour::Contour) and [`Shape`](crate::shape::Shape). They
represent connected sequences of segments and provide the core functionality
required by composite geometries.

Most users should interact with this module through the [`Polysegment`] type
itself; see its documentation for details on construction, invariants, and
usage.
*/

use std::collections::VecDeque;

use crate::primitive::Primitive;
use crate::segment::{ArcSegment, LineSegment, Segment, SegmentRef};
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS};
use crate::{Transformation, composite::*};
use approx::ulps_eq;
use bounding_box::{BoundingBox, ToBoundingBox};
use rayon::prelude::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/**
A sequence of [`Segment`]s where each segment is "connected" to its successor
(end / stop point of a segment is approximately equal to the start point of the
successor).

*/
#[doc = ""]
#[cfg_attr(
    feature = "doc-images",
    doc = "![Polysegment example][example_polysegment]"
)]
#[cfg_attr(
    feature = "doc-images",
    embed_doc_image::embed_doc_image("example_polysegment", "docs/img/example_polysegment.svg")
)]
#[cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images not enabled**. Compile docs with `cargo doc --features 'doc-images'` and Rust version >= 1.54."
)]
/**

The approximate equality of start and end point is defined by
[`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`]:

```ignore
approx::ulps_eq!(polysegment[i].stop(), polysegment[i+1].start(), epsilon = DEFAULT_EPSILON, max_ulps = DEFAULT_MAX_ULPS)
```

A polysegment can intersect itself (i.e. at least two of its segments
intersect each other). This can be tested using its [`Composite`] trait
implementation:

```
use planar_geo::prelude::*;

let polysegment = Polysegment::from_points(&[[0.0, 0.0], [2.0, 2.0], [0.0, 2.0], [2.0, 0.0]]);
assert_eq!(polysegment.intersections_polysegment(&polysegment, 0.0, 0).count(), 1);
```

# Constructing a polysegment

A polysegment is essentially a [newtype](https://doc.rust-lang.org/rust-by-example/generics/new_types.html)
wrapper around a [`VecDeque<Segment>`] and therefore exposes many of the same
methods. For example, a polysegment can be built via
[`push_front`](Polysegment::push_front) and
[`push_back`](Polysegment::push_back) just like a [`VecDeque`]. The
[`extend_front`](Polysegment::extend_front) and
[`extend_back`](Polysegment::extend_back) methods can be used to implicitly
add [`LineSegment`]s.

There are multiple constructors available:
- [`new`](Polysegment::new): A new, empty polysegment with no segment.
- [`with_capacity`](Polysegment::with_capacity): A new, empty polysegment
with a defined space for segments preallocated with no segment.
- [`from_points`](Polysegment::from_points): A polysegment consisting of
[`LineSegment`]s which connect the given points.
- [`from_iter`](Polysegment::from_iter): A polysegment consisting of the
segments given by an iterator. In case two subsequent segments are not
connected, a "filler" [`LineSegment`] is introduced between them.

# Modifying a polysegment

To uphold the "connection" property, it is not possible to manipulate arbitrary
segments, which is why methods like [`VecDeque::get_mut`] are missing.
It is however possible to manipulate the whole polysegment via the
[`Transformation`] implementation. Additionally, individual segments can be
removed from the ends of the polysegment with [`pop_front`](Polysegment::pop_front)
and [`pop_back`](Polysegment::pop_back).

# Access of individual segments

Accessing individual segments is possible via indexing (which panics when
out-of-bounds) and the [`get`](Polysegment::get) method (which returns `None`
when out-of-bounds). As in the underlying deque, the segments are not
necessarily contiguous in memory and can generally only be accessed via two
slices with the [`as_slices`](Polysegment::as_slices) method. [`Polysegment`]
does however expose the [`make_contiguous`](Polysegment::make_contiguous) to
store the segments contiguous in memory.

# Serialization and deserialization

When the `serde` feature is enabled, a polysegment can be serialized and
deserialized using the [`serde`] crate. It uses the same serialized
representation as a [`VecDeque`].
 */
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polysegment(VecDeque<Segment>);

impl Polysegment {
    /**
    Creates an empty [`Polysegment`].

    # Examples

    ```
    use planar_geo::prelude::Polysegment;

    let polysegment = Polysegment::new();
    ```
     */
    pub fn new() -> Self {
        return Self(VecDeque::new());
    }

    /**
    Creates an empty [`Polysegment`] with space for at least `capacity`
    [`Segment`].

    # Examples

    ```
    use planar_geo::prelude::Polysegment;

    let polysegment = Polysegment::with_capacity(10);
    ```
     */
    pub fn with_capacity(capacity: usize) -> Self {
        return Self(VecDeque::with_capacity(capacity));
    }

    /**
    Creates an [`Polysegment`] of [`LineSegment`]s from the given points.
    If two consecutive points are equal (within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`]), they are treated as a single
    vertex.

    # Examples

    ```
    use planar_geo::prelude::Polysegment;

    // Three points -> Two line segments
    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(polysegment.len(), 2);

    // Two consecutive points are equal
    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]]);
    assert_eq!(polysegment.len(), 1);
    ```
     */
    pub fn from_points(points: &[[f64; 2]]) -> Self {
        let mut segments: VecDeque<Segment> = VecDeque::with_capacity(points.len());
        for window in points.windows(2) {
            let start = window[0];
            let stop = window[1];
            match segments.back() {
                Some(last_segment) => {
                    if last_segment.stop() == stop {
                        continue;
                    }
                }
                None => (),
            }

            if let Ok(ls) = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                segments.push_back(ls.into());
            }
        }
        return Self(segments);
    }

    /**
    "Closes" the [`Polysegment`] by connecting the start point of the first /
    "front" segment with the stop point of the last / "back" segment with a
    [`LineSegment`] (if the two points aren't already equal).

    # Examples

    ```
    use planar_geo::prelude::Polysegment;

    // Three points -> Two line segments
    let mut polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(polysegment.len(), 2);

    // Close the polysegment
    polysegment.close();
    assert_eq!(polysegment.len(), 3);
    ```
     */
    pub fn close(&mut self) -> () {
        if let Some(start) = self.0.front() {
            if let Some(stop) = self.0.back() {
                if let Ok(ls) = LineSegment::new(
                    stop.stop(),
                    start.start(),
                    DEFAULT_EPSILON,
                    DEFAULT_MAX_ULPS,
                ) {
                    self.0.push_back(ls.into())
                }
            }
        }
    }

    /**
    Returns a reference to the underlying [`VecDeque`].

    This allows using all methods of [`VecDeque`] which use a shared reference,
    e.g. its iterators, accessor methods etc.
     */
    pub fn vec_deque(&self) -> &VecDeque<Segment> {
        return &self.0;
    }

    /**
    Calls the [`VecDeque::as_slices`] method of the underlying deque. See its
    docstring for more.
     */
    pub fn as_slices(&self) -> (&[Segment], &[Segment]) {
        return self.0.as_slices();
    }

    /**
    Calls the [`VecDeque::make_contiguous`] method of the underlying deque. See
    its docstring for more.
     */
    pub fn make_contiguous(&mut self) -> &[Segment] {
        return self.0.make_contiguous();
    }

    /**
    Appends a [`Segment`] to the back of the [`Polysegment`].

    If the polysegment already has a "back" segment ([`Polysegment::back`] returns
    [`Some`]), the start point of `segment` is compared to the stop point of
    the back segment. If they aren't equal within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`], a filler line segment is
    inserted between the two.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();
    assert_eq!(polysegment.len(), 0);

    polysegment.push_back(LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 1);

    // The start point of this segment matches the stop point of the already
    // inserted segment
    polysegment.push_back(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 2);

    // Start and stop point don't match -> A filler segment is introduced
    polysegment.push_back(LineSegment::new([2.0, 1.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 4);
    ```
     */
    pub fn push_back(&mut self, segment: Segment) {
        let opt_stop = self.back().map(Segment::stop);
        if let Some(stop) = opt_stop {
            if let Ok(ls) =
                LineSegment::new(stop, segment.start(), DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            {
                self.0.push_back(ls.into());
            }
        }
        self.0.push_back(segment);
    }

    /**
    Prepends a [`Segment`] to the front of the [`Polysegment`].

    If the polysegment already has a "front" segment ([`Polysegment::front`] returns
    [`Some`]), the stop point of `segment` is compared to the start point of
    the front segment. If they aren't equal within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`], a filler line segment is
    inserted between the two.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();
    assert_eq!(polysegment.len(), 0);

    polysegment.push_front(LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 1);

    // The stop point of this segment matches the start point of the already
    // inserted segment
    polysegment.push_front(LineSegment::new([0.0, 1.0], [0.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 2);

    // Start and stop point don't match -> A filler segment is introduced
    polysegment.push_front(LineSegment::new([2.0, 1.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 4);
    ```
     */
    pub fn push_front(&mut self, segment: Segment) {
        let opt_start = self.front().map(Segment::start);
        if let Some(start) = opt_start {
            if let Ok(ls) =
                LineSegment::new(segment.stop(), start, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            {
                self.0.push_front(ls.into());
            }
        }
        self.0.push_front(segment);
    }

    /**
    Removes the last / back element from the polysegment and returns it, or [`None`]
    if it is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();

    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    polysegment.push_back(ls.clone());
    assert_eq!(polysegment.len(), 1);

    assert_eq!(polysegment.pop_back(), Some(ls));
    assert_eq!(polysegment.pop_back(), None);
     */
    pub fn pop_back(&mut self) -> Option<Segment> {
        return self.0.pop_back();
    }

    /**
    Removes the first / front element from the polysegment and returns it, or [`None`]
    if it is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();

    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    polysegment.push_front(ls.clone());
    assert_eq!(polysegment.len(), 1);

    assert_eq!(polysegment.pop_front(), Some(ls));
    assert_eq!(polysegment.pop_front(), None);
     */
    pub fn pop_front(&mut self) -> Option<Segment> {
        return self.0.pop_front();
    }

    /**
    Provides a reference to the back element, or [`None`] if the polysegment is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let ls: Segment = LineSegment::new([1.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(polysegment.back(), Some(&ls));
    ```
     */
    pub fn back(&self) -> Option<&Segment> {
        return self.0.back();
    }

    /**
    Provides a reference to the front element, or [`None`] if the polysegment is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(polysegment.front(), Some(&ls));
     */
    pub fn front(&self) -> Option<&Segment> {
        return self.0.front();
    }

    /**
    Adds a [`LineSegment`] to the front of `self` which stops at `point` and
    starts at the current stop point of `self` - i.e., the `stop` point of the
    [`Segment`] returned from [`Polysegment::back`]. If `self` is empty or
    `point` is equal to `stop`, this is a no-op.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();
    assert_eq!(polysegment.len(), 0);

    // No-op, since the polysegment is empty
    polysegment.extend_back([0.0, 0.0]);
    assert_eq!(polysegment.len(), 0);

    // Now add a line
    polysegment.push_back(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 1);

    // Now adding the point works
    polysegment.extend_back([0.0, 0.0]);
    assert_eq!(polysegment.len(), 2);
    ```
     */
    pub fn extend_back(&mut self, point: [f64; 2]) {
        let endpoint = self.back().map(|s| s.stop());
        if let Some(endpoint) = endpoint {
            if let Ok(ls) = LineSegment::new(endpoint, point, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                self.push_back(ls.into());
            }
        }
    }

    /**
    Adds a [`LineSegment`] to the front of `self` which starts at `point` and
    stops at the current start point of `self` - i.e., the `start` point of the
    [`Segment`] returned from [`Polysegment::front`]. If `self` is empty or
    `point` is equal to `start`, this is a no-op.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();
    assert_eq!(polysegment.len(), 0);

    // No-op, since the polysegment is empty
    polysegment.extend_front([0.0, 0.0]);
    assert_eq!(polysegment.len(), 0);

    // Now add a line
    polysegment.push_front(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(polysegment.len(), 1);

    // Now adding the point works
    polysegment.extend_front([0.0, 0.0]);
    assert_eq!(polysegment.len(), 2);
    ```
     */
    pub fn extend_front(&mut self, point: [f64; 2]) {
        let endpoint = self.front().map(|s| s.start());
        if let Some(endpoint) = endpoint {
            if let Ok(ls) = LineSegment::new(point, endpoint, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                self.push_front(ls.into());
            }
        }
    }

    /**
    Returns whether the polysegment / the underlying [`VecDeque`] is empty or not.

    # Examples

    ```
     use planar_geo::prelude::*;

    let mut polysegment = Polysegment::new();
    assert!(polysegment.is_empty());

    polysegment.push_back(LineSegment::new([0.0, 0.0], [1.0, 0.0], 0.0, 0).unwrap().into());
    assert!(!polysegment.is_empty());
    ```
     */
    pub fn is_empty(&self) -> bool {
        return self.0.is_empty();
    }

    /**
    Provides a reference to the [`Segment`] at the given index.

    The segment at index 0 is the front of the polysegment.

    # Examples

    ```
     use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    assert!(polysegment.get(0).is_some());
    assert!(polysegment.get(3).is_none());
    ```
     */
    pub fn get(&self, index: usize) -> Option<&Segment> {
        return self.0.get(index);
    }

    /**
    Returns the number of [`Segment`]s in `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    assert_eq!(polysegment.len(), 2);
    ```
     */
    pub fn len(&self) -> usize {
        return self.0.len();
    }

    /**
    Returns a front-to-back iterator over all [`Segment`]s of `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let mut iter = polysegment.segments();
    assert!(iter.next().is_some());
    assert!(iter.next().is_some());
    assert!(iter.next().is_none());
    ```
     */
    pub fn segments(&self) -> std::collections::vec_deque::Iter<'_, Segment> {
        return self.0.iter();
    }

    /**
    Returns a parallel front-to-back iterator over all [`Segment`]s of `self`.

    # Examples

    ```
    use rayon::prelude::*;
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let vec: Vec<_> = polysegment.segments_par().collect();
    assert_eq!(vec.len(), 2);
    ```
     */
    pub fn segments_par(&self) -> rayon::collections::vec_deque::Iter<'_, Segment> {
        return self.0.par_iter();
    }

    /**
    Moves all [`Segment`]s of `other` into `self`, leaving `other` empty.

    If the first point of `other` is not equal to the last point of `self`, a
    filler line segment is introduced first

    # Panics

    Panics if the new number of elements in `self` overflows an [`usize`].

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment1 = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let mut polysegment2 = Polysegment::from_points(&[[1.0, 1.0], [0.0, 1.0]]);
    let mut polysegment3 = Polysegment::from_points(&[[0.0, 2.0], [0.0, 3.0]]);

    assert_eq!(polysegment1.len(), 2);
    assert_eq!(polysegment2.len(), 1);
    assert_eq!(polysegment3.len(), 1);

    polysegment1.append(&mut polysegment2);
    assert_eq!(polysegment1.len(), 3);

    polysegment1.append(&mut polysegment3);
    assert_eq!(polysegment1.len(), 5);
    ```
     */
    pub fn append(&mut self, other: &mut Polysegment) -> () {
        if let Some(s) = other.0.front() {
            let start = s.start();
            if let Some(s) = self.0.back() {
                let stop = s.stop();
                if let Ok(ls) = LineSegment::new(stop, start, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                    self.0.push_back(ls.into())
                }
            }
        }

        // Append the segments
        self.0.append(&mut other.0);
    }

    /**
    Returns an iterator over all points of the polysegment.

    # Examples

    ```
    use planar_geo::prelude::*;

    let polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);

    let mut iter = polysegment.points();
    assert_eq!(iter.next(), Some([0.0, 0.0]));
    assert_eq!(iter.next(), Some([1.0, 0.0]));
    assert_eq!(iter.next(), Some([1.0, 1.0]));
    assert_eq!(iter.next(), None);
    ```
     */
    pub fn points(&self) -> PointIterator<'_> {
        return PointIterator::new(self, false, Polygonizer::default());
    }

    /**
    Returns the points of a polygon polysegment which approximates `self`. The
    individual segments are "polygonized" via [`Segment::polygonize`] and
    an [`SegmentPolygonizer`](crate::segment::SegmentPolygonizer) specified
    within [`Polygonizer`]. See the docstring of the latter fore more.
     */
    pub fn polygonize(&self, polygonizer: Polygonizer) -> PointIterator<'_> {
        return PointIterator::new(self, false, polygonizer);
    }

    /**
    Returns the combined length of all segments of `self`.

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::{PI, TAU, SQRT_2};

    // Square with a side length of 1
    let points = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut polysegment = Polysegment::from_points(points);
    polysegment.close();
    assert_eq!(polysegment.length(), 4.0);

    // Circle with a radius of 1
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, TAU, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();
    let polysegment = Polysegment::from(segment);
    assert_eq!(polysegment.length(), TAU);

    // Half-circle segment
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();
    let mut polysegment = Polysegment::from(segment);
    polysegment.close();
    assert_eq!(polysegment.length(), PI + 2.0);

    // Triangle
    let points = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut polysegment = Polysegment::from_points(points);
    polysegment.close();
    assert_eq!(polysegment.length(), 2.0 + SQRT_2);
    ```
     */
    pub fn length(&self) -> f64 {
        return self.segments_par().map(|segment| segment.length()).sum();
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
    let line = Polysegment::from_points(points);
    let cut = Polysegment::from_points(&[[-1.0, 1.0], [3.0, 1.0]]);

    let separated_lines = line.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

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
        // Inner helper function
        fn create_cut_segment(
            segment: &Segment,
            start: [f64; 2],
            stop: [f64; 2],
            center: [f64; 2],
            epsilon: f64,
            max_ulps: u32,
        ) -> Option<Segment> {
            // Check to avoid degenerated segments, when start equals stop
            if !ulps_eq!(start, stop, epsilon = epsilon, max_ulps = max_ulps) {
                match segment {
                    Segment::LineSegment(_) => {
                        if let Ok(ls) = LineSegment::new(start, stop, epsilon, max_ulps) {
                            return Some(ls.into());
                        } else {
                            return None;
                        }
                    }
                    Segment::ArcSegment(_) => {
                        let normalized_stop = [stop[0] - center[0], stop[1] - center[1]];
                        let normalized_start = [start[0] - center[0], start[1] - center[1]];
                        let stop_angle = normalized_stop[1].atan2(normalized_stop[0]);
                        let start_angle = normalized_start[1].atan2(normalized_start[0]);
                        let offset_angle = stop_angle - start_angle;

                        if let Ok(arc) = ArcSegment::from_start_center_angle(
                            start,
                            center,
                            offset_angle,
                            epsilon,
                            max_ulps,
                        ) {
                            return Some(arc.into());
                        } else {
                            return None;
                        }
                    }
                };
            } else {
                return None;
            }
        }

        let mut lines: Vec<Polysegment> = Vec::new();
        let mut segments_of_current_line: Vec<Segment> = Vec::new();

        // Temporary variables. The values themselves are dummy values to satisfy the
        // borrow checker.
        let mut center: [f64; 2] = [0.0, 0.0];

        // Storage for intersections
        let mut intersections: Vec<[f64; 2]> = Vec::new();

        for segment_self in self.0.iter() {
            // Preparations for the current segment
            match segment_self {
                Segment::LineSegment(_) => (),
                Segment::ArcSegment(arc) => {
                    center = arc.center();
                }
            }

            // Reset the intersections vector
            intersections.clear();

            for segment_other in other.0.iter() {
                let intersections_segment_other =
                    segment_self.intersections_primitive(segment_other, epsilon, max_ulps);
                intersections.extend(
                    intersections_segment_other
                        .into_iter()
                        .map::<[f64; 2], _>(From::from),
                );
            }

            if intersections.is_empty() {
                // No intersections could be found -> add segment_self to the list of segments
                // forming the current line
                segments_of_current_line.push(segment_self.clone())
            } else {
                // Sort the intersections depending on their distance from the start of
                // segment_self (smalles first)
                intersections.sort_unstable_by(|a, b| {
                    let start = segment_self.start();

                    // Calculate the distance of a from segment_self.start and compare it to the
                    // distance of b from segment_self. Since we are only
                    // interested in the order, it is not necessary to normalize the distance.
                    // Do the same for b.
                    let dist_a = (start[0] - a[0]).powi(2) + (start[1] - a[1]).powi(2);
                    let dist_b = (start[0] - b[0]).powi(2) + (start[1] - b[1]).powi(2);

                    // When non-comparable values are used (e.g. NaN), the user
                    // has provided nonsensical values for the segment points
                    // and the result of the intersection cut will be
                    // meaningless anyway -> We can return any ordering
                    dist_a
                        .partial_cmp(&dist_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                for (idx, intersection) in intersections.iter().enumerate() {
                    if idx == 0 {
                        let start = segment_self.start();
                        let stop = *intersection;

                        if let Some(segment) = create_cut_segment(
                            &segment_self,
                            start,
                            stop,
                            center,
                            epsilon,
                            max_ulps,
                        ) {
                            if segments_of_current_line.is_empty() {
                                lines.push(Polysegment::from(segment));
                            } else {
                                let mut polysegment =
                                    Polysegment::with_capacity(segments_of_current_line.len());
                                for s in segments_of_current_line.into_iter() {
                                    polysegment.push_back(s);
                                }
                                polysegment.push_back(segment);
                                lines.push(polysegment);
                            }
                        } else {
                            if !segments_of_current_line.is_empty() {
                                let mut polysegment =
                                    Polysegment::with_capacity(segments_of_current_line.len());
                                for s in segments_of_current_line.into_iter() {
                                    polysegment.push_back(s);
                                }
                                lines.push(polysegment);
                            }
                        }

                        // Reinstatiate the vector for the next segment_self
                        segments_of_current_line = Vec::new();
                    } else {
                        // SAFETY: An out-of-bound index is catched by the if idx == 0 above
                        let start = unsafe { *intersections.get_unchecked(idx - 1) };
                        let stop = *intersection;

                        // Found some intersections => finish the current line and create a new one
                        if let Some(segment) = create_cut_segment(
                            &segment_self,
                            start,
                            stop,
                            center,
                            epsilon,
                            max_ulps,
                        ) {
                            lines.push(Polysegment::from(segment));
                        }
                    };
                }

                // Populate the "tail end" (the last cut)
                if let Some(start) = intersections.last() {
                    if let Some(segment) = create_cut_segment(
                        &segment_self,
                        start.clone(),
                        segment_self.stop(),
                        center,
                        epsilon,
                        max_ulps,
                    ) {
                        segments_of_current_line.push(segment);
                    }
                }
            };
        }

        // Store the last segments as a line into the lines vector
        let mut polysegment = Polysegment::with_capacity(segments_of_current_line.len());
        for s in segments_of_current_line.into_iter() {
            polysegment.push_back(s);
        }
        lines.push(polysegment);

        return lines;
    }

    /**
    Reverses all [`Segment`]s of the polysegment by switching their start and stop
    points (see [`Segment::reverse`]). To ensure the "connected" property holds
    true, the ordering of the segments itself is exchanged as well.

    This method calls [`make_contiguous`](Polysegment::make_contiguous) so all
    segments are in the first slice before reversing it.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(polysegment.front().unwrap().start(), [0.0, 0.0]);
    assert_eq!(polysegment.back().unwrap().stop(), [0.0, 1.0]);

    polysegment.reverse();
    assert_eq!(polysegment.front().unwrap().start(), [0.0, 1.0]);
    assert_eq!(polysegment.back().unwrap().stop(), [0.0, 0.0]);
    ```
     */
    pub fn reverse(&mut self) -> () {
        self.make_contiguous();

        // After making the deque contiguous, the second slice is empty
        let (slice, _) = self.0.as_mut_slices();
        for segment in slice.iter_mut() {
            segment.reverse();
        }
        slice.reverse();
    }

    /**
    Creates a rotated pattern from `self`. For each one of the specified
    `repetitions`, the individual segments of `self` are cloned, rotated around
    `center` and then pushed to the back of `self`. The rotation angle is
    `angle` times the index of the current repetition plus one.

    # Examples

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::FRAC_PI_2;

    let mut polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0]]);

    // 0 repetitions => The original polysegment is left unchanged
    polysegment.rotational_pattern([1.0, 1.0], FRAC_PI_2, 0);
    let points: Vec<[f64; 2]> = polysegment.points().collect();
    assert_eq!(points.len(), 2);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);

    // Two repetitions: The points are rotated by 90° and 180° respectively.
    // The polysegment has 6 points after the extension
    polysegment.rotational_pattern([1.0, 1.0], FRAC_PI_2, 2);
    let points: Vec<[f64; 2]> = polysegment.points().collect();
    assert_eq!(points.len(), 6);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[2], [2.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[3], [2.0, 1.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[4], [2.0, 2.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[5], [1.0, 2.0], epsilon = 1e-10);
    ```
    */
    pub fn rotational_pattern(&mut self, center: [f64; 2], angle: f64, repetitions: usize) -> () {
        let num_segments = self.num_segments();
        for rep in 0..repetitions {
            for seg_idx in 0..num_segments {
                // Index is always in bounds, because its maximum value is
                // the number of segments at the start of the outer loop
                // and no segments are removed from self.
                let mut segment = self[seg_idx].clone();
                segment.rotate(center, angle * (rep as f64 + 1.0));
                self.push_back(segment);
            }
        }
    }

    /**
    Creates a translated pattern from `self`. For each one of the specified
    `repetitions`, the individual segments of `self` are cloned, shifted and
    then pushed to the back of `self`. The shift vector is `shift` times the
    index of the current repetition plus one.

    # Examples

    ```
    use planar_geo::prelude::*;
    use std::f64::consts::FRAC_PI_2;

    let mut polysegment = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0]]);

    // 0 repetitions => The original polysegment is left unchanged
    polysegment.translational_pattern([1.0, 1.0], 0);
    let points: Vec<[f64; 2]> = polysegment.points().collect();
    assert_eq!(points.len(), 2);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);

    // Two repetitions: A "stair" polysegment is created
    polysegment.translational_pattern([1.0, 1.0], 2);
    let points: Vec<[f64; 2]> = polysegment.points().collect();
    assert_eq!(points.len(), 6);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[2], [1.0, 1.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[3], [2.0, 1.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[4], [2.0, 2.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[5], [3.0, 2.0], epsilon = 1e-10);
    ```
    */
    pub fn translational_pattern(&mut self, shift: [f64; 2], repetitions: usize) -> () {
        let num_segments = self.num_segments();
        for rep in 0..repetitions {
            let f = rep as f64 + 1.0; // f for "factor"
            for seg_idx in 0..num_segments {
                // Index is always in bounds, because its maximum value is
                // the number of segments at the start of the outer loop
                // and no segments are removed from self.
                let mut segment = self[seg_idx].clone();
                segment.translate([shift[0] * f, shift[1] * f]);
                self.push_back(segment);
            }
        }
    }
}

impl FromIterator<Segment> for Polysegment {
    fn from_iter<T: IntoIterator<Item = Segment>>(iter: T) -> Self {
        let mut polysegment = Polysegment::new();
        for segment in iter {
            polysegment.push_back(segment);
        }
        return polysegment;
    }
}

impl crate::composite::private::Sealed for Polysegment {}

impl Composite for Polysegment {
    fn segment(&self, key: SegmentKey) -> Option<&crate::segment::Segment> {
        return self.get(key.segment_idx);
    }

    fn num_segments(&self) -> usize {
        return self.len();
    }

    fn centroid(&self) -> [f64; 2] {
        return crate::CentroidData::from(self).into();
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self
            .segments()
            .enumerate()
            .map(|(i, s)| (SegmentKey::from_segment_idx(i), s));
    }

    fn par_iter<'a>(
        &'a self,
    ) -> impl ParallelIterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self
            .segments_par()
            .enumerate()
            .map(|(i, s)| (SegmentKey::from_segment_idx(i), s));
    }

    fn intersections_polysegment<'a>(
        &'a self,
        other: &'a Self,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a {
        let same_polysegment_len = if std::ptr::eq(self, other) {
            Some(self.num_segments())
        } else {
            None
        };

        // Naive implementation: Loop over all segments of segments_self. Check each
        // segment of self for an intersection with all segments of
        // segments_other.
        return other
            .0
            .iter()
            .enumerate()
            .flat_map(move |(right_idx, right_seg)| {
                intersections_between_polysegment_and_segment_priv(
                    self.0.iter(),
                    right_idx,
                    right_seg,
                    same_polysegment_len,
                    epsilon,
                    max_ulps,
                )
            });
    }

    fn intersections_polysegment_par<'a>(
        &'a self,
        other: &'a Self,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        let same_polysegment_len = if std::ptr::eq(self, other) {
            Some(self.num_segments())
        } else {
            None
        };

        // Naive implementation: Loop over all segments of segments_self. Check each
        // segment of self for an intersection with all segments of
        // segments_other. The outer iteration is done in parallel,
        // therefore the intersections are out of order.
        return other
            .0
            .par_iter()
            .enumerate()
            .flat_map(move |(right_idx, right_seg)| {
                intersections_between_polysegment_and_segment_priv_par(
                    self.0.par_iter(),
                    right_idx,
                    right_seg,
                    same_polysegment_len,
                    epsilon,
                    max_ulps,
                )
            });
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> {
        shape
            .intersections_polysegment(self, epsilon, max_ulps)
            .map(Intersection::switch)
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> {
        shape
            .intersections_polysegment_par(self, epsilon, max_ulps)
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
            .intersections_polysegment(self, epsilon, max_ulps)
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
            .intersections_polysegment_par(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        return self
            .intersections_primitive(&point, epsilon, max_ulps)
            .next()
            .is_some();
    }

    fn covers_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        let segment: SegmentRef = segment.into();
        return covers_segment(self.0.iter(), segment, epsilon, max_ulps);
    }

    fn contains_point(&self, _: [f64; 2], _: f64, _: u32) -> bool {
        // A polysegment has no surface area and therefore cannot contain a point.
        return false;
    }

    fn contains_segment<'a, T: Into<SegmentRef<'a>>>(&self, _: T, _: f64, _: u32) -> bool {
        return false;
    }

    fn covers_composite<'a, T: Composite + Sync>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        return other.covers_polysegment(self, epsilon, max_ulps);
    }

    fn contains_composite<'a, T: Composite + Sync>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> bool {
        return other.contains_polysegment(self, epsilon, max_ulps);
    }
}

impl std::ops::Index<usize> for Polysegment {
    type Output = Segment;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl From<Segment> for Polysegment {
    fn from(value: Segment) -> Self {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(value);
        return polysegment;
    }
}

impl From<LineSegment> for Polysegment {
    fn from(value: LineSegment) -> Self {
        return Polysegment::from(Segment::LineSegment(value));
    }
}

impl From<ArcSegment> for Polysegment {
    fn from(value: ArcSegment) -> Self {
        return Polysegment::from(Segment::ArcSegment(value));
    }
}

impl crate::Transformation for Polysegment {
    fn rotate(&mut self, center: [f64; 2], angle: f64) -> () {
        self.0.par_iter_mut().for_each(|segment| {
            segment.rotate(center, angle);
        })
    }

    fn translate(&mut self, shift: [f64; 2]) -> () {
        self.0.par_iter_mut().for_each(|segment| {
            segment.translate(shift);
        })
    }

    fn scale(&mut self, factor: f64) -> () {
        self.0.par_iter_mut().for_each(|segment| {
            segment.scale(factor);
        })
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        self.0.par_iter_mut().for_each(|segment| {
            segment.line_reflection(start, stop);
        })
    }
}

impl ToBoundingBox for Polysegment {
    fn bounding_box(&self) -> BoundingBox {
        // Use the bounding box of the first segment as a starting point
        return self
            .segments_par()
            .skip(1)
            .map(|segment| BoundingBox::from(segment))
            .reduce(
                || {
                    if let Some(s) = self.get(0) {
                        BoundingBox::from(s)
                    } else {
                        BoundingBox::new(0.0, 0.0, 0.0, 0.0)
                    }
                },
                |prev, curr| prev.union(&curr),
            );
    }
}

impl IntoIterator for Polysegment {
    type Item = Segment;
    type IntoIter = std::collections::vec_deque::IntoIter<Segment>;

    fn into_iter(self) -> Self::IntoIter {
        return self.0.into_iter();
    }
}

impl From<Polysegment> for VecDeque<Segment> {
    fn from(line: Polysegment) -> Self {
        return line.0;
    }
}

impl From<BoundingBox> for Polysegment {
    fn from(bounding_box: BoundingBox) -> Self {
        return Self::from(&bounding_box);
    }
}

impl From<&BoundingBox> for Polysegment {
    fn from(bounding_box: &BoundingBox) -> Self {
        return Polysegment::from_points(&[
            [bounding_box.xmin(), bounding_box.ymin()],
            [bounding_box.xmax(), bounding_box.ymin()],
            [bounding_box.xmax(), bounding_box.ymax()],
            [bounding_box.xmin(), bounding_box.ymax()],
            [bounding_box.xmin(), bounding_box.ymin()],
        ]);
    }
}

impl From<&Polysegment> for crate::CentroidData {
    fn from(value: &Polysegment) -> Self {
        return value
            .segments_par()
            .map(|segment| match segment {
                Segment::LineSegment(line) => Self::from(line),
                Segment::ArcSegment(arc) => Self::from(arc),
            })
            .reduce(
                || Self {
                    area: 0.0,
                    x: 0.0,
                    y: 0.0,
                },
                |prev, curr| prev.union(&curr),
            );
    }
}

/**
Calculates the area of a polygon defined by the points given by an iterator.
If the polygon orientation is mathematically positive (points are given
counterclockwise), then the result is also positive. Otherwise, it is negative.

This implementation uses the "Shoelace formula", as e.g. described here:
<https://en.wikipedia.org/wiki/Shoelace_formula>.

```
use planar_geo::polysegment::area_signed;

let pt1 = [0.0, 0.0];
let pt2 = [1.0, 0.0];
let pt3 = [1.0, 1.0];
assert_eq!(
    area_signed([pt1, pt2, pt3].into_iter()),
    0.5
);
assert_eq!(
    area_signed([pt2, pt1, pt3].into_iter()),
    -0.5
);
```
 */
pub fn area_signed<'a, I>(mut points: I) -> f64
where
    I: Iterator<Item = [f64; 2]>,
{
    // Keep the first vertex
    let first_vertex = match points.next() {
        Some(vertex) => vertex,
        None => return 0.0,
    };
    let mut previous_vertex = first_vertex.clone();

    let mut area = 0.0;

    // The polysegmented element covers the end value x_n*y_0 - x_0*y_n, where n
    // equals the last iterator element
    for current_vertex in points.chain(std::iter::once(first_vertex)) {
        area += (previous_vertex[1] + current_vertex[1]) * (previous_vertex[0] - current_vertex[0]);
        previous_vertex = current_vertex.clone();
    }

    return 0.5 * area; // Correction factor of 0.5 comes from the area formel of a triangle.
}

pub(crate) fn covers_segment<'a, 'b, I>(
    iterator: I,
    segment: SegmentRef<'b>,
    epsilon: f64,
    max_ulps: u32,
) -> bool
where
    I: Iterator<Item = &'a Segment>,
{
    match segment {
        SegmentRef::LineSegment(line_segment) => {
            // Multiple subsequent line segments are combined into a single
            // one if they form a single line segment (i.e. their angles are
            // identical)
            let mut combined: Option<LineSegment> = None;
            for seg in iterator {
                if let Segment::LineSegment(sl) = seg {
                    combined = match combined {
                        Some(c) => {
                            if ulps_eq!(
                                c.angle(),
                                sl.angle(),
                                epsilon = epsilon,
                                max_ulps = max_ulps
                            ) {
                                // Add the new line segment, if it has the same angle as
                                // "combined"
                                if let Ok(new_combined) =
                                    LineSegment::new(c.start(), sl.stop(), epsilon, max_ulps)
                                {
                                    if new_combined.covers(line_segment, epsilon, max_ulps) {
                                        return true;
                                    }
                                    Some(new_combined)
                                } else {
                                    None
                                }
                            } else {
                                // Replace the old combined segment with the new line segment
                                if sl.covers(line_segment, epsilon, max_ulps) {
                                    return true;
                                }
                                Some(sl.clone())
                            }
                        }
                        None => {
                            if sl.covers(line_segment, epsilon, max_ulps) {
                                return true;
                            }
                            Some(sl.clone())
                        }
                    };
                } else {
                    // Arc segment breaks the chain
                    combined = None;
                }
            }
            return false;
        }
        SegmentRef::ArcSegment(arc_segment) => {
            // Multiple subsequent arc segments are combined into a single
            // one if they form a single arc segment
            // (i.e. their radius and center are identical).
            let mut combined: Option<ArcSegment> = None;
            for seg in iterator {
                if let Segment::ArcSegment(sa) = seg {
                    combined = match combined {
                        Some(c) => {
                            if ulps_eq!(
                                c.center(),
                                sa.center(),
                                epsilon = epsilon,
                                max_ulps = max_ulps
                            ) && ulps_eq!(
                                c.radius(),
                                sa.radius(),
                                epsilon = epsilon,
                                max_ulps = max_ulps
                            ) {
                                // Add the new line segment, if it has the same angle as
                                // "combined"
                                if let Ok(new_combined) =
                                    ArcSegment::from_center_radius_start_offset_angle(
                                        c.center(),
                                        c.radius(),
                                        c.start_angle(),
                                        c.offset_angle() + sa.offset_angle(),
                                        epsilon,
                                        max_ulps,
                                    )
                                {
                                    if new_combined.covers(arc_segment, epsilon, max_ulps) {
                                        return true;
                                    }
                                    Some(new_combined)
                                } else {
                                    None
                                }
                            } else {
                                // Replace the old combined segment with the new line segment
                                if sa.covers(arc_segment, epsilon, max_ulps) {
                                    return true;
                                }
                                Some(sa.clone())
                            }
                        }
                        None => {
                            if sa.covers(arc_segment, epsilon, max_ulps) {
                                return true;
                            }
                            Some(sa.clone())
                        }
                    };
                } else {
                    // Line segment breaks the chain
                    combined = None;
                }
            }
            return false;
        }
    }
}

fn intersections_between_polysegment_and_segment_priv<'a, L>(
    left: L,
    right_idx: usize,
    right_seg: &'a Segment,
    same_polysegment_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> impl Iterator<Item = Intersection> + 'a
where
    L: Iterator<Item = &'a Segment> + 'a,
{
    left.enumerate()
        .filter_map(move |(left_idx, left_seg)| {
            segment_intersections(
                left_seg,
                left_idx,
                right_seg,
                right_idx,
                same_polysegment_len,
                epsilon,
                max_ulps,
            )
        })
        .flatten()
}

fn intersections_between_polysegment_and_segment_priv_par<'a, L>(
    left: L,
    right_idx: usize,
    right_seg: &'a Segment,
    same_polysegment_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> impl ParallelIterator<Item = Intersection> + 'a
where
    L: IndexedParallelIterator<Item = &'a Segment> + 'a,
{
    left.enumerate()
        .filter_map(move |(left_idx, left_seg)| {
            segment_intersections(
                left_seg,
                left_idx,
                right_seg,
                right_idx,
                same_polysegment_len,
                epsilon,
                max_ulps,
            )
            .map(|iter| iter.par_bridge())
        })
        .flatten()
}

fn segment_intersections<'a>(
    left_seg: &'a Segment,
    left_idx: usize,
    right_seg: &'a Segment,
    right_idx: usize,
    same_polysegment_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> Option<impl Iterator<Item = Intersection> + 'a> {
    /*
    If the two segments are identical, they have no intersection by definition.

    If slices_are_identical, an additional optimization is possible:
    It is sufficient to only check segments where left_idx < right_idx.
    This is possible because a naive loop implementation visits every segment pair twice:
    a) Left segment is in the outer loop, right segment is in the inner loop
    b) Left segment is in the outer loop, right segment is in the inner loop
     */
    if same_polysegment_len.is_some() {
        if left_idx == right_idx {
            return None;
        }
        if left_idx >= right_idx {
            return None;
        }
    } else {
        if left_seg == right_seg {
            return None;
        }
    }

    let intersection_iter = left_seg
        .intersections_primitive(right_seg, epsilon, max_ulps)
        .into_iter()
        .filter_map(move |point| {
            let point: [f64; 2] = point.into();

            /*
            If the two slices are identical, check whether left_seg is a
            successor / predecessor of right_seg. If yes, filter out the
            connection points.
             */

            if let Some(len) = same_polysegment_len {
                /*
                This check evaluates whether the left segment is a predecessor
                of the right segment. If that is the case, filter out
                "intersections" which are actually just the connection point
                between the two segments.
                 */
                if left_idx + 1 == right_idx || (left_idx == 0 && right_idx == len - 1) {
                    if approx::ulps_eq!(
                        point,
                        left_seg.stop(),
                        epsilon = epsilon,
                        max_ulps = max_ulps
                    ) {
                        return None;
                    }
                }
            }

            Some(Intersection {
                point,
                left: SegmentKey::from_segment_idx(left_idx),
                right: SegmentKey::from_segment_idx(right_idx),
            })
        });

    return Some(intersection_iter);
}
