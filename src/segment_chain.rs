/*!
Defines the [`SegmentChain`] type, a foundational [`Composite`] used throughout
this crate.

Segment chains serve as the basis for higher-level geometric types such as
[`Contour`](crate::contour::Contour) and [`Shape`](crate::shape::Shape). They
represent connected sequences of segments and provide the core functionality
required by composite geometries.

Most users should interact with this module through the [`SegmentChain`] type
itself; see its documentation for details on construction, invariants, and
usage.
*/

use std::collections::VecDeque;

use crate::Transformation;
use crate::composite::*;
use crate::primitive::Primitive;
use crate::segment::Segment;
use crate::segment::arc_segment::ArcSegment;
use crate::segment::line_segment::LineSegment;
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS};
use approx::ulps_eq;
use bounding_box::BoundingBox;
use rayon::prelude::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/**
A sequence of [`Segment`]s where each segment is "connected" to its successor
(end / stop point of a segment is approximately equal to the start point of the
successor).

*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/example_segment_chain.svg \"Segment chain\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/example_segment_chain.svg)"
)]
/**

The approximate equality of start and end point is defined by
[`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`]:

```ignore
approx::ulps_eq!(chain[i].stop(), chain[i+1].start(), epsilon = DEFAULT_EPSILON, max_ulps = DEFAULT_MAX_ULPS)
```

A segment chain can intersect itself (i.e. at least two of its segments
intersect each other). This can be tested using its [`Composite`] trait
implementation:

```
use planar_geo::prelude::*;

let chain = SegmentChain::from_points(&[[0.0, 0.0], [2.0, 2.0], [0.0, 2.0], [2.0, 0.0]]);
assert_eq!(chain.intersections_segment_chain(&chain, 0.0, 0).count(), 1);
```

# Constructing a segment chain

A segment chain is essentially a [newtype](https://doc.rust-lang.org/rust-by-example/generics/new_types.html)
wrapper around a [`VecDeque<Segment>`] and therefore exposes many of the same
methods. For example, a segment chain can be built via
[`push_front`](SegmentChain::push_front) and
[`push_back`](SegmentChain::push_back) just like a [`VecDeque`]. The
[`extend_front`](SegmentChain::extend_front) and
[`extend_back`](SegmentChain::extend_back) methods can be used to implicitly
add [`LineSegment`]s.

There are multiple constructors available:
- [`new`](SegmentChain::new): A new, empty segment chain with no segment.
- [`with_capacity`](SegmentChain::with_capacity): A new, empty segment chain
with a defined space for segments preallocated with no segment.
- [`from_points`](SegmentChain::from_points): A segment chain consisting of
[`LineSegment`]s which connect the given points.
- [`from_iter`](SegmentChain::from_iter): A segment chain consisting of the
segments given by an iterator. In case two subsequent segments are not
connected, a "filler" [`LineSegment`] is introduced between them.

# Modifying a segment chain

To uphold the "connection" property, it is not possible to manipulate arbitrary
segments, which is why methods like [`VecDeque::get_mut`] are missing.
It is however possible to manipulate the whole segment chain via the
[`Transformation`] implementation. Additionally, individual segments can be
removed from the ends of the chain with [`pop_front`](SegmentChain::pop_front)
and [`pop_back`](SegmentChain::pop_back).

# Access of individual segments

Accessing individual segments is possible via indexing (which panics when
out-of-bounds) and the [`get`](SegmentChain::get) method (which returns `None`
when out-of-bounds). As in the underlying deque, the segments are not
necessarily contiguous in memory and can generally only be accessed via two
slices with the [`as_slices`](SegmentChain::as_slices) method. [`SegmentChain`]
does however expose the [`make_contiguous`](SegmentChain::make_contiguous) to
store the segments contiguous in memory.

# Serialization and deserialization

When the `serde` feature is enabled, a segment chain can be serialized and
deserialized using the [`serde`] crate. It uses the same serialized
representation as a [`VecDeque`].
 */
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SegmentChain(VecDeque<Segment>);

impl SegmentChain {
    /**
    Creates an empty [`SegmentChain`].

    # Examples

    ```
    use planar_geo::prelude::SegmentChain;

    let chain = SegmentChain::new();
    ```
     */
    pub fn new() -> Self {
        return Self(VecDeque::new());
    }

    /**
    Creates an empty [`SegmentChain`] with space for at least `capacity`
    [`Segment`].

    # Examples

    ```
    use planar_geo::prelude::SegmentChain;

    let chain = SegmentChain::with_capacity(10);
    ```
     */
    pub fn with_capacity(capacity: usize) -> Self {
        return Self(VecDeque::with_capacity(capacity));
    }

    /**
    Creates an [`SegmentChain`] of [`LineSegment`]s from the given points.
    If two consecutive points are equal (within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`]), they are treated as a single
    vertex.

    # Examples

    ```
    use planar_geo::prelude::SegmentChain;

    // Three points -> Two line segments
    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(chain.len(), 2);

    // Two consecutive points are equal
    let chain = SegmentChain::from_points(&[[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]]);
    assert_eq!(chain.len(), 1);
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
    "Closes" the [`SegmentChain`] by connecting the start point of the first /
    "front" segment with the stop point of the last / "back" segment with a
    [`LineSegment`] (if the two points aren't already equal).

    # Examples

    ```
    use planar_geo::prelude::SegmentChain;

    // Three points -> Two line segments
    let mut chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(chain.len(), 2);

    // Close the chain
    chain.close();
    assert_eq!(chain.len(), 3);
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
    Appends a [`Segment`] to the back of the [`SegmentChain`].

    If the chain already has a "back" segment ([`SegmentChain::back`] returns
    [`Some`]), the start point of `segment` is compared to the stop point of
    the back segment. If they aren't equal within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`], a filler line segment is
    inserted between the two.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();
    assert_eq!(chain.len(), 0);

    chain.push_back(LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 1);

    // The start point of this segment matches the stop point of the already
    // inserted segment
    chain.push_back(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 2);

    // Start and stop point don't match -> A filler segment is introduced
    chain.push_back(LineSegment::new([2.0, 1.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 4);
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
    Prepends a [`Segment`] to the front of the [`SegmentChain`].

    If the chain already has a "front" segment ([`SegmentChain::front`] returns
    [`Some`]), the stop point of `segment` is compared to the start point of
    the front segment. If they aren't equal within the tolerances defined by
    [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`], a filler line segment is
    inserted between the two.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();
    assert_eq!(chain.len(), 0);

    chain.push_front(LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 1);

    // The stop point of this segment matches the start point of the already
    // inserted segment
    chain.push_front(LineSegment::new([0.0, 1.0], [0.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 2);

    // Start and stop point don't match -> A filler segment is introduced
    chain.push_front(LineSegment::new([2.0, 1.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 4);
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
    Removes the last / back element from the chain and returns it, or [`None`]
    if it is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();

    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    chain.push_back(ls.clone());
    assert_eq!(chain.len(), 1);

    assert_eq!(chain.pop_back(), Some(ls));
    assert_eq!(chain.pop_back(), None);
     */
    pub fn pop_back(&mut self) -> Option<Segment> {
        return self.0.pop_back();
    }

    /**
    Removes the first / front element from the chain and returns it, or [`None`]
    if it is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();

    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    chain.push_front(ls.clone());
    assert_eq!(chain.len(), 1);

    assert_eq!(chain.pop_front(), Some(ls));
    assert_eq!(chain.pop_front(), None);
     */
    pub fn pop_front(&mut self) -> Option<Segment> {
        return self.0.pop_front();
    }

    /**
    Provides a reference to the back element, or [`None`] if the chain is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let ls: Segment = LineSegment::new([1.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(chain.back(), Some(&ls));
    ```
     */
    pub fn back(&self) -> Option<&Segment> {
        return self.0.back();
    }

    /**
    Provides a reference to the front element, or [`None`] if the chain is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let ls: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(chain.front(), Some(&ls));
     */
    pub fn front(&self) -> Option<&Segment> {
        return self.0.front();
    }

    /**
    Adds a [`LineSegment`] to the front of `self` which stops at `point` and
    starts at the current stop point of `self` - i.e., the `stop` point of the
    [`Segment`] returned from [`SegmentChain::back`]. If `self` is empty or
    `point` is equal to `stop`, this is a no-op.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();
    assert_eq!(chain.len(), 0);

    // No-op, since the chain is empty
    chain.extend_back([0.0, 0.0]);
    assert_eq!(chain.len(), 0);

    // Now add a line
    chain.push_back(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 1);

    // Now adding the point works
    chain.extend_back([0.0, 0.0]);
    assert_eq!(chain.len(), 2);
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
    [`Segment`] returned from [`SegmentChain::front`]. If `self` is empty or
    `point` is equal to `start`, this is a no-op.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();
    assert_eq!(chain.len(), 0);

    // No-op, since the chain is empty
    chain.extend_front([0.0, 0.0]);
    assert_eq!(chain.len(), 0);

    // Now add a line
    chain.push_front(LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into());
    assert_eq!(chain.len(), 1);

    // Now adding the point works
    chain.extend_front([0.0, 0.0]);
    assert_eq!(chain.len(), 2);
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
    Returns whether the chain / the underlying [`VecDeque`] is empty or not.

    # Examples

    ```
     use planar_geo::prelude::*;

    let mut chain = SegmentChain::new();
    assert!(chain.is_empty());

    chain.push_back(LineSegment::new([0.0, 0.0], [1.0, 0.0], 0.0, 0).unwrap().into());
    assert!(!chain.is_empty());
    ```
     */
    pub fn is_empty(&self) -> bool {
        return self.0.is_empty();
    }

    /**
    Provides a reference to the [`Segment`] at the given index.

    The segment at index 0 is the front of the chain.

    # Examples

    ```
     use planar_geo::prelude::*;

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    assert!(chain.get(0).is_some());
    assert!(chain.get(3).is_none());
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

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    assert_eq!(chain.len(), 2);
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

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let mut iter = chain.iter();
    assert!(iter.next().is_some());
    assert!(iter.next().is_some());
    assert!(iter.next().is_none());
    ```
     */
    pub fn iter(&self) -> std::collections::vec_deque::Iter<'_, Segment> {
        return self.0.iter();
    }

    /**
    Returns a parallel front-to-back iterator over all [`Segment`]s of `self`.

    # Examples

    ```
    use rayon::prelude::*;
    use planar_geo::prelude::*;

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let vec: Vec<_> = chain.par_iter().collect();
    assert_eq!(vec.len(), 2);
    ```
     */
    pub fn par_iter(&self) -> rayon::collections::vec_deque::Iter<'_, Segment> {
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

    let mut chain1 = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let mut chain2 = SegmentChain::from_points(&[[1.0, 1.0], [0.0, 1.0]]);
    let mut chain3 = SegmentChain::from_points(&[[0.0, 2.0], [0.0, 3.0]]);

    assert_eq!(chain1.len(), 2);
    assert_eq!(chain2.len(), 1);
    assert_eq!(chain3.len(), 1);

    chain1.append(&mut chain2);
    assert_eq!(chain1.len(), 3);

    chain1.append(&mut chain3);
    assert_eq!(chain1.len(), 5);
    ```
     */
    pub fn append(&mut self, other: &mut SegmentChain) -> () {
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
    Returns an iterator over all points of the chain.

    # Examples

    ```
    use planar_geo::prelude::*;

    let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);

    let mut iter = chain.points();
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
    Returns the points of a polygon chain which approximates `self`. The
    individual segments are "polygonized" via [`Segment::polygonize`] and
    an [`SegmentPolygonizer`] specified within [`Polygonizer`]. See the
    docstring of the latter fore more.
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
    let mut segment_chain = SegmentChain::from_points(points);
    segment_chain.close();
    assert_eq!(segment_chain.length(), 4.0);

    // Circle with a radius of 1
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, TAU, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();
    let segment_chain = SegmentChain::from(segment);
    assert_eq!(segment_chain.length(), TAU);

    // Half-circle segment
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();
    let mut segment_chain = SegmentChain::from(segment);
    segment_chain.close();
    assert_eq!(segment_chain.length(), PI + 2.0);

    // Triangle
    let points = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut segment_chain = SegmentChain::from_points(points);
    segment_chain.close();
    assert_eq!(segment_chain.length(), 2.0 + SQRT_2);
    ```
     */
    pub fn length(&self) -> f64 {
        return self.par_iter().map(|segment| segment.length()).sum();
    }

    /**
    Cuts `self` into multiple chains by intersecting it with `other` and returns
    those them.

    The specified tolerances `epsilon` and `max_ulps` are used to determine
    intersections, see documentation of
    [`PrimitiveIntersections`](crate::primitive::PrimitiveIntersections).

    # Examples

    ```
    use planar_geo::prelude::*;

    let points = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let line = SegmentChain::from_points(points);
    let cut = SegmentChain::from_points(&[[-1.0, 1.0], [3.0, 1.0]]);

    let separated_lines = line.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    // This cut results in two separate chains
    assert_eq!(separated_lines.len(), 2);
    ```
     */
    pub fn intersection_cut(
        &self,
        other: &SegmentChain,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<SegmentChain> {
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

        let mut lines: Vec<SegmentChain> = Vec::new();
        let mut segments_of_current_line: Vec<Segment> = Vec::new();

        // Temporary variables. The values themselves are dummy values to satisfy the borrow checker.
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
                // No intersections could be found -> add segment_self to the list of segments forming the current line
                segments_of_current_line.push(segment_self.clone())
            } else {
                // Sort the intersections depending on their distance from the start of segment_self (smalles first)
                intersections.sort_unstable_by(|a, b| {
                    let start = segment_self.start();

                    // Calculate the distance of a from segment_self.start and compare it to the distance of b from segment_self.
                    // Since we are only interested in the order, it is not necessary to normalize the distance.
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
                                lines.push(SegmentChain::from(segment));
                            } else {
                                let mut chain =
                                    SegmentChain::with_capacity(segments_of_current_line.len());
                                for s in segments_of_current_line.into_iter() {
                                    chain.push_back(s);
                                }
                                chain.push_back(segment);
                                lines.push(chain);
                            }
                        } else {
                            if !segments_of_current_line.is_empty() {
                                let mut chain =
                                    SegmentChain::with_capacity(segments_of_current_line.len());
                                for s in segments_of_current_line.into_iter() {
                                    chain.push_back(s);
                                }
                                lines.push(chain);
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
                            lines.push(SegmentChain::from(segment));
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
        let mut chain = SegmentChain::with_capacity(segments_of_current_line.len());
        for s in segments_of_current_line.into_iter() {
            chain.push_back(s);
        }
        lines.push(chain);

        return lines;
    }

    /**
    Reverses all [`Segment`]s of the chain by switching their start and stop
    points (see [`Segment::reverse`]). To ensure the "connected" property holds
    true, the ordering of the segments itself is exchanged as well.

    This method calls [`make_contiguous`](SegmentChain::make_contiguous) so all
    segments are in the first slice before reversing it.

    # Examples

    ```
    use planar_geo::prelude::*;

    let mut chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    assert_eq!(chain.front().unwrap().start(), [0.0, 0.0]);
    assert_eq!(chain.back().unwrap().stop(), [0.0, 1.0]);

    chain.reverse();
    assert_eq!(chain.front().unwrap().start(), [0.0, 1.0]);
    assert_eq!(chain.back().unwrap().stop(), [0.0, 0.0]);
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

    let mut chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0]]);

    // 0 repetitions => The original chain is left unchanged
    chain.rotational_pattern([1.0, 1.0], FRAC_PI_2, 0);
    let points: Vec<[f64; 2]> = chain.points().collect();
    assert_eq!(points.len(), 2);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);

    // Two repetitions: The points are rotated by 90° and 180° respectively.
    // The chain has 6 points after the extension
    chain.rotational_pattern([1.0, 1.0], FRAC_PI_2, 2);
    let points: Vec<[f64; 2]> = chain.points().collect();
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

    let mut chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0]]);

    // 0 repetitions => The original chain is left unchanged
    chain.translational_pattern([1.0, 1.0], 0);
    let points: Vec<[f64; 2]> = chain.points().collect();
    assert_eq!(points.len(), 2);
    approx::assert_abs_diff_eq!(points[0], [0.0, 0.0], epsilon = 1e-10);
    approx::assert_abs_diff_eq!(points[1], [1.0, 0.0], epsilon = 1e-10);

    // Two repetitions: A "stair" segment chain is created
    chain.translational_pattern([1.0, 1.0], 2);
    let points: Vec<[f64; 2]> = chain.points().collect();
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

impl FromIterator<Segment> for SegmentChain {
    fn from_iter<T: IntoIterator<Item = Segment>>(iter: T) -> Self {
        let mut chain = SegmentChain::new();
        for segment in iter {
            chain.push_back(segment);
        }
        return chain;
    }
}

impl Composite for SegmentChain {
    type SegmentKey = SegmentIdx;

    fn segment(&self, key: Self::SegmentKey) -> Option<&crate::segment::Segment> {
        return self.get(key.0);
    }

    fn num_segments(&self) -> usize {
        return self.len();
    }

    fn centroid(&self) -> [f64; 2] {
        return crate::CentroidData::from(self).into();
    }

    fn intersections_primitive<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, ()>> + 'a {
        self.iter()
            .enumerate()
            .map(move |(idx, s)| {
                s.intersections_primitive(primitive, epsilon, max_ulps)
                    .into_iter()
                    .map(move |point| Intersection {
                        point,
                        left: SegmentIdx(idx),
                        right: (),
                    })
            })
            .flatten()
    }

    fn intersections_primitive_par<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ()>> + 'a {
        self.par_iter()
            .enumerate()
            .map(move |(idx, s)| {
                s.intersections_primitive(primitive, epsilon, max_ulps)
                    .into_iter()
                    .par_bridge()
                    .map(move |point| Intersection {
                        point,
                        left: SegmentIdx(idx),
                        right: (),
                    })
            })
            .flatten()
    }

    fn intersections_segment_chain<'a>(
        &'a self,
        other: &'a Self,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        let same_segment_chain_len = if std::ptr::eq(self, other) {
            Some(self.num_segments())
        } else {
            None
        };

        // Naive implementation: Loop over all segments of segments_self. Check each segment of self
        // for an intersection with all segments of segments_other.
        return other
            .0
            .iter()
            .enumerate()
            .flat_map(move |(right_idx, right_seg)| {
                intersections_between_segment_chain_and_segment_priv(
                    self.0.iter(),
                    right_idx,
                    right_seg,
                    same_segment_chain_len,
                    epsilon,
                    max_ulps,
                )
            });
    }

    fn intersections_segment_chain_par<'a>(
        &'a self,
        other: &'a Self,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        let same_segment_chain_len = if std::ptr::eq(self, other) {
            Some(self.num_segments())
        } else {
            None
        };

        // Naive implementation: Loop over all segments of segments_self. Check each segment of self
        // for an intersection with all segments of segments_other. The outer iteration is done in parallel,
        // therefore the intersections are out of order.
        return other
            .0
            .par_iter()
            .enumerate()
            .flat_map(move |(right_idx, right_seg)| {
                intersections_between_segment_chain_and_segment_priv_par(
                    self.0.par_iter(),
                    right_idx,
                    right_seg,
                    same_segment_chain_len,
                    epsilon,
                    max_ulps,
                )
            });
    }

    fn intersections_contour<'a>(
        &'a self,
        contour: &'a crate::prelude::Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> {
        return self.intersections_segment_chain(contour.segment_chain(), epsilon, max_ulps);
    }

    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a crate::prelude::Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        return self.intersections_segment_chain_par(contour.segment_chain(), epsilon, max_ulps);
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> {
        shape
            .intersections_segment_chain(self, epsilon, max_ulps)
            .map(Intersection::switch)
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> {
        shape
            .intersections_segment_chain_par(self, epsilon, max_ulps)
            .map(Intersection::switch)
    }

    fn intersections_composite<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, T::SegmentKey>> + 'a
    where
        Self: Sized,
    {
        return other
            .intersections_segment_chain(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn intersections_composite_par<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, T::SegmentKey>> + 'a
    where
        Self: Sized,
        <T as Composite>::SegmentKey: Send,
    {
        return other
            .intersections_segment_chain_par(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        return self
            .intersections_primitive(&point, epsilon, max_ulps)
            .next()
            .is_some();
    }
}

impl std::ops::Index<usize> for SegmentChain {
    type Output = Segment;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl From<Segment> for SegmentChain {
    fn from(value: Segment) -> Self {
        let mut chain = SegmentChain::new();
        chain.push_back(value);
        return chain;
    }
}

impl From<LineSegment> for SegmentChain {
    fn from(value: LineSegment) -> Self {
        return SegmentChain::from(Segment::LineSegment(value));
    }
}

impl From<ArcSegment> for SegmentChain {
    fn from(value: ArcSegment) -> Self {
        return SegmentChain::from(Segment::ArcSegment(value));
    }
}

impl Transformation for SegmentChain {
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

impl From<&SegmentChain> for BoundingBox {
    fn from(value: &SegmentChain) -> BoundingBox {
        // Use the bounding box of the first segment as a starting point
        return value
            .par_iter()
            .skip(1)
            .map(|segment| BoundingBox::from(segment))
            .reduce(
                || {
                    if let Some(s) = value.get(0) {
                        BoundingBox::from(s)
                    } else {
                        BoundingBox::new(0.0, 0.0, 0.0, 0.0)
                    }
                },
                |prev, curr| prev.union(&curr),
            );
    }
}

impl IntoIterator for SegmentChain {
    type Item = Segment;
    type IntoIter = std::collections::vec_deque::IntoIter<Segment>;

    fn into_iter(self) -> Self::IntoIter {
        return self.0.into_iter();
    }
}

fn intersections_between_segment_chain_and_segment_priv<'a, L>(
    left: L,
    right_idx: usize,
    right_seg: &'a Segment,
    same_segment_chain_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> impl Iterator<Item = Intersection<SegmentIdx, SegmentIdx>> + 'a
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
                same_segment_chain_len,
                epsilon,
                max_ulps,
            )
        })
        .flatten()
}

fn intersections_between_segment_chain_and_segment_priv_par<'a, L>(
    left: L,
    right_idx: usize,
    right_seg: &'a Segment,
    same_segment_chain_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> impl ParallelIterator<Item = Intersection<SegmentIdx, SegmentIdx>> + 'a
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
                same_segment_chain_len,
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
    same_segment_chain_len: Option<usize>,
    epsilon: f64,
    max_ulps: u32,
) -> Option<impl Iterator<Item = Intersection<SegmentIdx, SegmentIdx>> + 'a> {
    /*
    If the two segments are identical, they have no intersection by definition.

    If slices_are_identical, an additional optimization is possible:
    It is sufficient to only check segments where left_idx < right_idx.
    This is possible because a naive loop implementation visits every segment pair twice:
    a) Left segment is in the outer loop, right segment is in the inner loop
    b) Left segment is in the outer loop, right segment is in the inner loop
     */
    if same_segment_chain_len.is_some() {
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

            if let Some(len) = same_segment_chain_len {
                /*
                This check evaluates whether the left segment is a predecessor
                of the right segment. If that is the case, filter out
                "intersections" which are actually just the connection point
                between the two segments.
                 */
                if left_idx + 1 == right_idx || (left_idx == 0 && right_idx == len - 1) {
                    if ulps_eq!(
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
                left: SegmentIdx(left_idx),
                right: SegmentIdx(right_idx),
            })
        });

    return Some(intersection_iter);
}

impl From<SegmentChain> for VecDeque<Segment> {
    fn from(line: SegmentChain) -> Self {
        return line.0;
    }
}

impl From<BoundingBox> for SegmentChain {
    fn from(bounding_box: BoundingBox) -> Self {
        return Self::from(&bounding_box);
    }
}

impl From<&BoundingBox> for SegmentChain {
    fn from(bounding_box: &BoundingBox) -> Self {
        return SegmentChain::from_points(&[
            [bounding_box.xmin(), bounding_box.ymin()],
            [bounding_box.xmax(), bounding_box.ymin()],
            [bounding_box.xmax(), bounding_box.ymax()],
            [bounding_box.xmin(), bounding_box.ymax()],
            [bounding_box.xmin(), bounding_box.ymin()],
        ]);
    }
}

impl From<&SegmentChain> for crate::CentroidData {
    fn from(value: &SegmentChain) -> Self {
        return value
            .par_iter()
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
https://en.wikipedia.org/wiki/Shoelace_formula

```
use planar_geo::segment_chain::area_signed;

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

    // The chained element covers the end value x_n*y_0 - x_0*y_n, where n equals the last iterator element
    for current_vertex in points.chain(std::iter::once(first_vertex)) {
        area += (previous_vertex[1] + current_vertex[1]) * (previous_vertex[0] - current_vertex[0]);
        previous_vertex = current_vertex.clone();
    }

    return 0.5 * area; // Correction factor of 0.5 comes from the area formel of a triangle.
}
