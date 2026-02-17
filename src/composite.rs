/*!
This module contains the definition of the [`Composite`] trait and the
[`Intersection`], the [`SegmentIdx`] and the [`ShapeIdx`] structs.

The geometric types in this crate are either "primitives" (which implement the
[`Primitive`] trait) or "composites" (which are defined from multiple primitives
and implement the [`Composite`](crate::composite::Composite) trait). See the
trait documentation for details.
 */

use rayon::prelude::ParallelIterator;

use crate::contour::Contour;
use crate::primitive::Primitive;
use crate::segment::{Segment, SegmentPolygonizer};
use crate::segment_chain::SegmentChain;
use crate::shape::Shape;

/**
Intersection point and indices of the "left" and "right" composite types.

This struct is used as the [`Iterator::Item`] of the iterators created by the
intersection methods of the [`Composite`] trait. It contains the intersection
point itself as well as the segment keys of the involved composite types.
These keys ([`SegmentIdx`] or [`ShapeIdx`]) can be used to retrieve the segments
from the composite types where the intersection occurred (see
[`Composite::segment`]).
 */
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Intersection<L, R> {
    /// The intersection point itself.
    pub point: [f64; 2],
    /// Key to retrieve the intersected segment of the "left" side composite
    /// (`self` for the intersection functions in [`Composite`]).
    pub left: L,
    /// Key to retrieve the intersected segment of the "right" side composite
    /// (the second argument of the intersection functions in [`Composite`]).
    pub right: R,
}

impl<L, R> Intersection<L, R> {
    /**
    Switches the `left` and ``right` fields of `Self`.
     */
    pub fn switch(self) -> Intersection<R, L> {
        return Intersection {
            point: self.point,
            left: self.right,
            right: self.left,
        };
    }
}

impl From<Intersection<(), ()>> for [f64; 2] {
    fn from(value: Intersection<(), ()>) -> Self {
        value.point
    }
}

impl From<[f64; 2]> for Intersection<(), ()> {
    fn from(value: [f64; 2]) -> Self {
        Intersection {
            point: value,
            left: (),
            right: (),
        }
    }
}

/**
A key to access a segment from a [`SegmentChain`] or a [`Contour`].

This is just a wrapper around an [`usize`] used to index into the underlying
vector of [`Segment`](crate::segment::Segment)s which stores the data of the
aforementioned [`Composite`] types. Its main function is to provide type safety
when used e.g. as part of an [`Intersection`] or when accessing segments via the
[`Composite::segment`] method.

# Examples

```
use planar_geo::prelude::*;

let raw_idx = 2;
let chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
assert_eq!(chain.segment(SegmentIdx(2)), chain.get(raw_idx));
```
 */
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SegmentIdx(pub usize);

impl From<usize> for SegmentIdx {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<SegmentIdx> for usize {
    fn from(value: SegmentIdx) -> Self {
        value.0
    }
}

/**
A key to access a segment from a [`Shape`].

A [`Shape`] is composed of multiple [`Contour`]s, which in turn are composed of
multiple [`Segment`](crate::segment::Segment)s. The field
[`ShapeIdx::contour_idx`] specifies the contour where the adressed segment is
located, while [`ShapeIdx::segment_idx`] retrieves the segment itself from the
corresponding contour.

# Examples

```
use planar_geo::prelude::*;

let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let contour = Contour::new(SegmentChain::from_points(vertices));
let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
let hole = Contour::new(SegmentChain::from_points(vertices));
let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

// Comparison of using the key vs. manual retrieval
let contour_idx = 1;
let segment_idx = 2;
let manually_retrieved = shape.contours().get(contour_idx).map(|c|c.segment_chain().get(segment_idx)).flatten();

let key = ShapeIdx {
    contour_idx,
    segment_idx: SegmentIdx(segment_idx)
};
assert_eq!(shape.segment(key), manually_retrieved);
```
 */
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShapeIdx {
    /// Index of the contour where the segment is located.
    pub contour_idx: usize,
    /// Index of the segment itself within the contour specified by
    /// [`ShapeIdx::contour_idx`].
    pub segment_idx: SegmentIdx,
}

impl ShapeIdx {
    /**
    Creates a new key out of any two values which can be converted into [`usize`].

    # Examples

    ```
    use planar_geo::composite::{ShapeIdx, SegmentIdx};

    // Create with from types used inside the struct directly
    let idx = ShapeIdx::new(1usize, SegmentIdx(2));
    assert_eq!(idx.contour_idx, 1);
    assert_eq!(idx.segment_idx, SegmentIdx(2));

    // From two u32
    let idx = ShapeIdx::new(0u16, 3u16);
    assert_eq!(idx.contour_idx, 0);
    assert_eq!(idx.segment_idx, SegmentIdx(3));
    ```
     */
    pub fn new<C: Into<usize>, S: Into<usize>>(contour_idx: C, segment_idx: S) -> Self {
        return Self {
            contour_idx: contour_idx.into(),
            segment_idx: SegmentIdx(segment_idx.into()),
        };
    }
}

impl<L: PartialEq, R: PartialEq> approx::AbsDiffEq for Intersection<L, R> {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        return self.point.abs_diff_eq(&other.point, epsilon)
            && self.left == other.left
            && self.right == other.right;
    }
}

impl<L: PartialEq, R: PartialEq> approx::RelativeEq for Intersection<L, R> {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        return self.point.relative_eq(&other.point, epsilon, max_relative)
            && self.left == other.left
            && self.right == other.right;
    }
}

impl<L: PartialEq, R: PartialEq> approx::UlpsEq for Intersection<L, R> {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        return self.point.ulps_eq(&other.point, epsilon, max_ulps)
            && self.left == other.left
            && self.right == other.right;
    }
}

/**
A trait for "composite" types: [`SegmentChain`]s, [`Contour`]s and [`Shape`]s.

This [`Composite`] trait provides a common innterface for intersection
calculation between a composite and a primitive or other composites.

Dufferent to primitives, there can be an arbitrary number of intersections
between two intersections, which is why the corresponding methods return
an iterator. Since intersection calculation between composites can be easily
parallelized, there is a serial and a parallel variant (the latter returning
a [`ParallelIterator`] and having a `_par` suffix).`

These iterators always return an [`Intersection`] object. The "left" side of the
object refers to the first argument `self`, whereas the "right" side
refers to the type of the other composite. Hence, the type of the left side is
an associated type [`Composite::SegmentKey`] which depends on the implementor,
while the right side is determined by the other composite. See the docstring of
[`Intersection`] for more.

In contrast to primitives, composites can self-intersect. The self-intersection
points can be calculated by using `self` also as the second argument `other`:

```
use planar_geo::prelude::*;

let sc = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]);

let mut iter = sc.intersections_segment_chain(&sc, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

// Intersection in the "cross" middle
assert_eq!(iter.next().unwrap().point, [0.5, 0.5]);

// The chain start and end point "touch", this is also counted as an intersection
assert_eq!(iter.next().unwrap().point, [0.0, 0.0]);
assert!(iter.next().is_none());

// By contrast, a contour is closed by default, hence start and end point are
// connected and therefore are not an intersection
let c = Contour::new(sc.clone());
let mut iter = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

// Intersection in the "cross" middle
assert_eq!(iter.next().unwrap().point, [0.5, 0.5]);

// No other intersection
assert!(iter.next().is_none());
```

The time complexity of calculating all intersections between two composites a
and b is O(n_a*n_b), where n_a/b is the value returned by
[`Composite::num_segments`] for a and b respectively. This is due to the fact
that the composite intersection functions compare each segment of a with each
segment of b (via [`Primitive::intersections_primitive`]).

As with the intersection methods of [`Primitive`], each intersection function
takes `epsilon` and `max_ulps` as additional arguments to specify a certain
tolerance for intersection detection.
 */
pub trait Composite {
    /**
    The segment key type needed to retrieve a segment from `Self` via
    [`Composite::segment`]. This type is the "left" part of the [`Intersection`]
    returned by the various intersection methods of the [`Composite`] trait and
    can be used to directly retrieve the segment where the intersection occurred.
     */
    type SegmentKey;

    /**
    Returns the segment associated with the given `key`, if it exists.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));

    // Drop the reference before reusing contour later
    {
        let segment = hole.segment(SegmentIdx(2));
        assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9], 0.0, 0).unwrap())));
    }

    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    // Self::SegmentKey of Shape is ShapeIdx. This accessor retrieves the third
    // segment of the second contour / first hole of the shape
    let segment = shape.segment(ShapeIdx {contour_idx: 1, segment_idx: SegmentIdx(2)});
    assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9], 0.0, 0).unwrap())));
    ```
     */
    fn segment(&self, key: Self::SegmentKey) -> Option<&crate::segment::Segment>;

    /**
    Returns the number of segments in the composite type.

    For a [`SegmentChain`] and a [`Contour`], this is equal to the length
    of the underlying [`Vec<Segment>`], whereas for a [`Shape`], this is equal
    to the number of all segments of all underlying [`Contour`]s.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));
    assert_eq!(contour.num_segments(), 4);

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));
    assert_eq!(hole.num_segments(), 4);

    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");
    assert_eq!(shape.num_segments(), 8);
    ```
     */
    fn num_segments(&self) -> usize;

    /**
    Calculates the centroid of `self`, i.e. its center of mass.

    The implementation of this function relies on the fact that the centroids
    of simple geometric bodies (such as the [`Segment`]s making up a
    [`SegmentChain`]) can be calculated from simple formulae. Those centroids
    can then be combined into that of a complex structure using the following
    formulae:
    `x = ∑ (xi * Ai) / ∑ Ai`
    `y = ∑ (yi * Ai) / ∑ Ai`
    where `xi` and `yi` are the centroid parameters and `Ai` is the surface area
    of the "segment area" `i` defined by connecting the segment end points to
    the origin.

    The centroids of the segment areas ar are calculated as follows:

    ## LineSegment
    Connect start and stop to the origin, then calculate the shape area as:
    `0.5 * ((stop[0] - start[0]) * (origin[1] - start[1]) - (origin[0] - start[0]) * (stop[1] - start[1]))`.
    The centroid is calculated as:
    `x = (start[0] + stop[0] + origin[0]) / 3` and `y = (start[1] + stop[1] + origin[1]) / 3`.

    ## ArcSegment
    Separate the arc into the following three shapes:
    1) triangle start -> stop -> origin
    2) circular segment start -> stop -> arc center
    3) triangle start -> stop -> center
    Calculate the areas and centroids of those segments, then summarize them as
    follows: `A = A1 + A2 - A3`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
    approx::assert_abs_diff_eq!(contour.centroid(), [0.5, 1.0], epsilon = 1e-3);

    let hole_1 = Contour::rectangle([0.6, 1.8], [0.9, 0.2]);
    let hole_2 = Contour::rectangle([0.1, 1.8], [0.4, 0.2]);
    let shape = Shape::new(vec![contour, hole_1, hole_2]).expect("holes do not intersect outer contour");

    approx::assert_abs_diff_eq!(shape.centroid(), [0.5, 1.0], epsilon = 1e-3);
    ```
     */
    fn centroid(&self) -> [f64; 2];

    /**
    Returns an iterator over all intersections of `self` with the `primitive`.

    The right side of the [`Intersection`]s created by the returned iterator
    is simply an empty tuple, because no index is needed to retrieve the
    primitive.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let chain = SegmentChain::from_points(vertices);

    let line = Line::from_point_angle([0.5, 0.5], 0.0);
    let mut intersections = chain.intersections_primitive(&line, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentIdx(1), right: ()})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentIdx(3), right: ()})
    );
    assert!(intersections.next().is_none());
    ```
     */
    fn intersections_primitive<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, ()>> + 'a;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `primitive`.

    This is the parallelized variant of [`Composite::intersections_primitive`].
    See its docstring for more information and examples.
     */
    fn intersections_primitive_par<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ()>> + 'a;

    /**
    Returns a iterator over all intersections of `self` with the `segment_chain`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let left = SegmentChain::from_points(vertices);

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let right = SegmentChain::from_points(vertices);

    let mut intersections = left.intersections_segment_chain(&right, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentIdx(1), right: SegmentIdx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentIdx(3), right: SegmentIdx(1)})
    );
    assert!(intersections.next().is_none());
    ```
    */
    fn intersections_segment_chain<'a>(
        &'a self,
        segment_chain: &'a SegmentChain,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `segment_chain`.

    This is the parallelized variant of [`Composite::intersections_segment_chain`].
    See its docstring for more information and examples.
     */
    fn intersections_segment_chain_par<'a>(
        &'a self,
        segment_chain: &'a SegmentChain,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a;

    /**
    Returns a iterator over all intersections of `self` with the `contour`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let chain = SegmentChain::from_points(vertices);

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));

    let mut intersections = chain.intersections_contour(&contour, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentIdx(1), right: SegmentIdx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentIdx(1), right: SegmentIdx(3)})
    );
    assert!(intersections.next().is_none());
    ```
    */
    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `contour`.

    This is the parallelized variant of [`Composite::intersections_contour`].
    See its docstring for more information and examples.
     */
    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a;

    /**
    Returns an iterator over all intersections of `self` with the `shape`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));
    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let chain = SegmentChain::from_points(vertices);

    let mut intersections = chain.intersections_shape(&shape, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {
            point: [1.0, 0.5],
            left: SegmentIdx(1),
            right: ShapeIdx {contour_idx: 0, segment_idx: SegmentIdx(1)}
        })
    );
    ```
     */
    fn intersections_shape<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> + 'a;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `contour`.

    This is the parallelized variant of
    [`Composite::intersections_shape`]. See its docstring for more
    information and examples.
     */
    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> + 'a;

    /**
    Returns the intersection between `self` and another type implementing
    [`Composite`].

    This is a generalized interface to all specialized intersection functions.
    For example, if `self` is a [`SegmentChain`], the implementation of this
    function boils down to:

    ```ignore
    impl Composite for SegmentChain {
        // Implementations of the other methods ...

        fn intersections_composite<'a, T: Composite>(
            &'a self,
            other: &'a T,
            epsilon: f64,
            max_ulps: u32,
        ) -> PrimitiveIntersections {
            other.intersections_segment_chain(self, epsilon, max_ulps)
        }
    }
    ```
     */
    fn intersections_composite<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, T::SegmentKey>> + 'a
    where
        Self: Sized;

    /**
    Returns a parallelized iterator over all intersections of `self` with `other`.

    This is the parallelized variant of [`Composite::intersections_composite`].
    See its docstring for more information and examples.
     */
    fn intersections_composite_par<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, T::SegmentKey>> + 'a
    where
        Self: Sized,
        <T as Composite>::SegmentKey: Send;

    /**
    Returns whether the given point is contained within the composite or not.

    The definition of "contained" depends on the composite type:
    - [`SegmentChain`]: A point is contained within it if it intersects any of
    its segments.
    - [`Contour`]: Like the [`SegmentChain`] definition, but additionally a
    point is also contained if it is within the enclosed surface.
    - [`Shape`]: Like the contour definition, but a point is not considered
    contained if it is inside one of the "hole" contours.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));
    let shape = Shape::new(vec![contour.clone(), hole]).expect("valid inputs");

    let pt = [0.5, 0.5];

    // pt is not on the segment chain defining the contour ...
    assert!(!contour.segment_chain().contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // ... but it is within the contour
    assert!(contour.contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // It is however not inside the shape, because it is inside the hole
    assert!(!shape.contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    ```
     */
    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool;
}

/**
This enum defines how a [`SegmentChain`] / [`Contour`] should be polygonized
(approximated by a polygon). It is the counterpart of [`SegmentPolygonizer`] and
built upon it.
 */
#[derive(Debug, Clone)]
pub enum Polygonizer {
    /**
    This variant defines a single [`SegmentPolygonizer`] for all
    [`ArcSegment`](crate::segment::ArcSegment)s and another single
    [`SegmentPolygonizer`] for all [`LineSegment`](crate::segment::LineSegment)s
    of a [`SegmentChain`] / [`Contour`].
     */
    PerType {
        /// [`SegmentPolygonizer`] for all arc segments.
        arc: SegmentPolygonizer,
        /// [`SegmentPolygonizer`] for all line segments.
        straight: SegmentPolygonizer,
    },
    /**
    This variant allows specifying individual [`SegmentPolygonizer`]s for each
    segment of a [`SegmentChain`] / [`Contour`] via a
    [`map`](Polygonizer::Individual::map). The map uses the segment index as its
    key. If no entry can be found within the map for a particular segment index,
    the fallback [`default`](Polygonizer::Individual::default) is used to
    polygonize the segment.
     */
    Individual {
        /// Fallback [`SegmentPolygonizer`] if there is no entry for a
        /// particular segment in the [`map`](Polygonizer::Individual::map).
        default: SegmentPolygonizer,
        /**
        Specifies the [`SegmentPolygonizer`] for each segment via an segment
        index -> polygonizer relationship.
         */
        map: std::collections::HashMap<usize, SegmentPolygonizer>,
    },
}

impl Polygonizer {
    /**
    Returns the [`SegmentPolygonizer`] for a `segment` with the given `index`.
     */
    pub fn segment_polygonizer(&self, segment: &Segment, index: usize) -> SegmentPolygonizer {
        match self {
            Polygonizer::PerType { arc, straight } => match segment {
                Segment::LineSegment(_) => return *straight,
                Segment::ArcSegment(_) => return *arc,
            },
            Polygonizer::Individual { default, map } => return *map.get(&index).unwrap_or(default),
        };
    }
}

impl Default for Polygonizer {
    fn default() -> Self {
        return Polygonizer::PerType {
            arc: SegmentPolygonizer::default(),
            straight: SegmentPolygonizer::default(),
        };
    }
}

/**
An iterator over all points contained within the
 */
#[derive(Clone, Debug)]
pub struct PointIterator<'a> {
    segment_chain: &'a SegmentChain,
    index: usize,
    skip_last_vertex: bool,
    polygonizer: Polygonizer,
    sub_iterator: Option<crate::segment::PolygonPointsIterator<'a>>,
}

impl<'a> PointIterator<'a> {
    pub(crate) fn new(
        segment_chain: &'a SegmentChain,
        skip_last_vertex: bool,
        polygonizer: Polygonizer,
    ) -> Self {
        let sub_iterator = segment_chain.front().map(|first_segment| {
            let sub_polygonizer = polygonizer.segment_polygonizer(first_segment, 0);

            first_segment.polygonize(sub_polygonizer)
        });
        return Self {
            segment_chain,
            index: 0,
            skip_last_vertex,
            polygonizer,
            sub_iterator,
        };
    }
}

impl<'a> Iterator for PointIterator<'a> {
    type Item = [f64; 2];
    fn next(&mut self) -> Option<[f64; 2]> {
        match self.sub_iterator.as_mut()?.next() {
            Some(pt) => return Some(pt),
            None => {
                // Replace the sub-iterator with that of the next segment
                self.index += 1;

                let segment = match self.segment_chain.get(self.index) {
                    Some(s) => s,
                    None => {
                        self.sub_iterator = None;
                        return None;
                    }
                };

                let sub_polygonizer = self.polygonizer.segment_polygonizer(segment, self.index);

                let mut sub_iterator = segment.polygonize(sub_polygonizer);
                sub_iterator.start_at_second_point();

                // If the final vertex of the segment_chain should be skipped, adjust the new sub iterator here
                if self.index + 1 == self.segment_chain.len() && self.skip_last_vertex {
                    sub_iterator.skip_last_vertex();
                }

                self.sub_iterator = Some(sub_iterator);
                return self.next();
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let length = self.segment_chain.len();
        (length, None)
    }
}
