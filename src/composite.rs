/*!
This module contains the definition of the [`Composite`] trait as well as the
[`Intersection`] and [`SegmentKey`] structs.

The geometric types in this crate are either "primitives" (which implement the
[`Primitive`] trait) or "composites" (which are defined from multiple primitives
and implement the [`Composite`] trait). See the trait documentation for details.
 */

use rayon::prelude::ParallelIterator;

use crate::contour::Contour;
use crate::polysegment::Polysegment;
use crate::primitive::Primitive;
use crate::segment::{Segment, SegmentPolygonizer};
use crate::shape::Shape;

/**
A key to access a [`Segment`] with the [`Composite::segment`] trait method.

The key consists of two indices: [`SegmentKey::contour_idx`] and
[`SegmentKey::segment_idx`]. The former index is used to select the
[`Contour`] from a [`Shape`]. In the implementations of other [`Composite`]
types ([`Polysegment`] and [`Contour`]), it is simply ignored.
[`SegmentKey::segment_idx`] is then used to access the [`Segment`] from the
[`Polysegment`] or [`Contour`].

If the [`contour_idx`](SegmentKey::contour_idx) is not needed, the convenience
constructor [`SegmentKey::from_segment_idx`] can be used (which sets the
[`contour_idx`](SegmentKey::contour_idx) to 0). The [`Default`] implementation
sets both indices to 0.

# Examples

```
use planar_geo::prelude::*;

let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let contour = Contour::new(Polysegment::from_points(vertices));

// Indices
let contour_idx = 1;
let segment_idx = 2;

// Comparison of using the key vs. manual retrieval
let manually_retrieved = contour.polysegment().get(segment_idx);

// Using the from_segment_idx constructor, the contour_idx is simply set to 0.
let key = SegmentKey::from_segment_idx(segment_idx);
assert_eq!(contour.segment(key), manually_retrieved);

let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
let hole = Contour::new(Polysegment::from_points(vertices));
let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

// Comparison of using the key vs. manual retrieval
let manually_retrieved = shape.contours().get(contour_idx).map(|c|c.polysegment().get(segment_idx)).flatten();

let key = SegmentKey {
    contour_idx,
    segment_idx,
};
assert_eq!(shape.segment(key), manually_retrieved);
```
 */
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct SegmentKey {
    /// Index of the [`Contour`] where the [`Segment`] is located. Ignored when
    ///the key is not used on a [`Shape`].
    pub contour_idx: usize,
    /// Index of the [`Segment`] itself within the [`Contour`] or
    /// [`Polysegment`].
    pub segment_idx: usize,
}

impl SegmentKey {
    /**
    Returns a new [`SegmentKey`] instance from its components.
     */
    pub fn new(contour_idx: usize, segment_idx: usize) -> Self {
        return Self {
            contour_idx,
            segment_idx,
        };
    }

    /**
    Returns a [`SegmentKey`] where [`SegmentKey::contour_idx`] is set to 0.
     */
    pub fn from_segment_idx(segment_idx: usize) -> Self {
        return Self {
            contour_idx: 0,
            segment_idx,
        };
    }
}

/**
Intersection between two [`Segment`]s of two geometric types consisting of the
[`Intersection::point`] itself and the keys to the involved segments.

This type is returned by all intersection methods of the [`Composite`] trait
and the geometry type enums in [crate::geometry]. It consists of the
[`Intersection::point`] and two [`SegmentKey`]s which can be used to retrieve
the segments which intersect each other. The [`Intersection::left`] key gets the
[`Segment`] of the first argument to the intersection method (usually `self`),
the [`Intersection::right`] key gets that of the second argument (often called
`&other`). If one of the arguments is not a [`Composite`], the corresponding key
has no use and is initialized to its default values (0 for both indices).

This struct implements [`approx::AbsDiffEq`], [`approx::RelativeEq`] and
[`approx::UlpsEq`] and can therefore be used in approximate comparisons:

```
use planar_geo::prelude::*;

let i1 = Intersection {
    point: [1.0, 0.0],
    left: SegmentKey::new(0, 1),
    right: SegmentKey::new(1, 2),
};
let i2 = Intersection {
    point: [1.1, 0.0],
    left: SegmentKey::new(0, 1),
    right: SegmentKey::new(1, 2),
};
assert_ne!(i1, i2);
approx::assert_abs_diff_eq!(i1, i2, epsilon = 0.5);

// Intersection point identical to i1, but left key different -> Unequal by default
let i3 = Intersection {
    point: [1.0, 0.0],
    left: SegmentKey::new(0, 0),
    right: SegmentKey::new(1, 2),
};
assert_ne!(i1, i3);
approx::assert_abs_diff_ne!(i1, i3, epsilon = 0.5);
```
 */
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Intersection {
    /// The intersection point itself.
    pub point: [f64; 2],
    /// Key to retrieve the intersected segment of the "left" side composite
    /// (`&self` for intersection methods).
    pub left: SegmentKey,
    /// Key to retrieve the intersected segment of the "right" side composite
    /// (second argument for intersection methods).
    pub right: SegmentKey,
}

impl Intersection {
    /**
    Exchanges [`Intersection::left`] and [`Intersection::right`].
     */
    pub fn switch(self) -> Intersection {
        Intersection {
            point: self.point,
            left: self.right,
            right: self.left,
        }
    }
}

impl From<Intersection> for [f64; 2] {
    fn from(value: Intersection) -> Self {
        value.point
    }
}

impl From<[f64; 2]> for Intersection {
    fn from(value: [f64; 2]) -> Self {
        Intersection {
            point: value,
            left: Default::default(),
            right: Default::default(),
        }
    }
}

impl approx::AbsDiffEq for Intersection {
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

impl approx::RelativeEq for Intersection {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        return self.point.relative_eq(&other.point, epsilon, max_relative)
            && self.left == other.left
            && self.right == other.right;
    }
}

impl approx::UlpsEq for Intersection {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        return self.point.ulps_eq(&other.point, epsilon, max_ulps)
            && self.left == other.left
            && self.right == other.right;
    }
}

impl Default for Intersection {
    fn default() -> Self {
        Self {
            point: Default::default(),
            left: Default::default(),
            right: Default::default(),
        }
    }
}

pub(crate) mod private {
    /// Sealed trait for [`Primitive`]
    pub trait Sealed {}
}

/**
A trait for "composite" types: [`Polysegment`]s, [`Contour`]s and [`Shape`]s.

This [`Composite`] trait provides a common innterface for intersection
calculation between a composite and a primitive or other composites.

Dufferent to primitives, there can be an arbitrary number of intersections
between two intersections, which is why the corresponding methods return
an iterator. Since intersection calculation between composites can be easily
parallelized, there is a serial and a parallel variant (the latter returning
a [`ParallelIterator`] and having a `_par` suffix).`

These iterators always return an [`Intersection`] object. The "left" side of the
object refers to the first argument `self`, whereas the "right" side
refers to the type of the other composite. See the docstring of [`Intersection`]
for more.

In contrast to primitives, composites can self-intersect. The self-intersection
points can be calculated by using `self` also as the second argument `other`:

```
use planar_geo::prelude::*;

let sc = Polysegment::from_points(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]);

let mut iter = sc.intersections_polysegment(&sc, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

// Intersection in the "cross" middle
assert_eq!(iter.next().unwrap().point, [0.5, 0.5]);

// The polysegment start and end point "touch", this is also counted as an intersection
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
pub trait Composite: private::Sealed {
    /**
    Returns the segment associated with the given `key`, if it exists.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));

    // Drop the reference before reusing contour later
    {
        let segment = hole.segment(SegmentKey::from_segment_idx(2));
        assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9], 0.0, 0).unwrap())));
    }

    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    // SegmentKey of Shape is SegmentKey. This accessor retrieves the third
    // segment of the second contour / first hole of the shape
    let segment = shape.segment(SegmentKey::new(1, 2));
    assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9], 0.0, 0).unwrap())));
    ```
     */
    fn segment(&self, key: SegmentKey) -> Option<&crate::segment::Segment>;

    /**
    Returns the number of segments in the composite type.

    For a [`Polysegment`] and a [`Contour`], this is equal to the length
    of the underlying [`Vec<Segment>`], whereas for a [`Shape`], this is equal
    to the number of all segments of all underlying [`Contour`]s.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    assert_eq!(contour.num_segments(), 4);

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
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
    [`Polysegment`]) can be calculated from simple formulae. Those centroids
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

    If eager collection of the returned [`Intersection`]s into a [`Vec`] is not
    an issue, consider using the more generic [`Composite::intersections`]
    instead (which can deal with both [`Primitive`]s and [`Composite`]s).

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let polysegment = Polysegment::from_points(vertices);

    let line = Line::from_point_angle([0.5, 0.5], 0.0);
    let mut intersections = polysegment.intersections_primitive(&line, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentKey::from_segment_idx(1), right: Default::default()})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentKey::from_segment_idx(3), right: Default::default()})
    );
    assert!(intersections.next().is_none());
    ```
     */
    fn intersections_primitive<'a, T: Primitive>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a;

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns a iterator over all intersections of `self` with the `polysegment`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let left = Polysegment::from_points(vertices);

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let right = Polysegment::from_points(vertices);

    let mut intersections = left.intersections_polysegment(&right, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentKey::from_segment_idx(3), right: SegmentKey::from_segment_idx(1)})
    );
    assert!(intersections.next().is_none());

    // Same result can be achieved with the generic method
    let mut intersections = left.intersections_composite(&right, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentKey::from_segment_idx(3), right: SegmentKey::from_segment_idx(1)})
    );
    assert!(intersections.next().is_none());
    ```
    */
    fn intersections_polysegment<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `polysegment`.

    This is the parallelized variant of [`Composite::intersections_polysegment`].
    See its docstring for more information and examples.
     */
    fn intersections_polysegment_par<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns a iterator over all intersections of `self` with the `contour`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));

    let mut intersections = polysegment.intersections_contour(&contour, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(3)})
    );
    assert!(intersections.next().is_none());

    // Same result can be achieved with the generic method
    let mut intersections = polysegment.intersections_composite(&contour, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [1.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(1)})
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {point: [0.0, 0.5], left: SegmentKey::from_segment_idx(1), right: SegmentKey::from_segment_idx(3)})
    );
    assert!(intersections.next().is_none());
    ```
    */
    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a;

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns an iterator over all intersections of `self` with the `shape`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let mut intersections = polysegment.intersections_shape(&shape, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {
            point: [1.0, 0.5],
            left: SegmentKey::from_segment_idx(1),
            right: SegmentKey::new(0, 1)
        })
    );

    // Same result can be achieved with the generic method
    let mut intersections = polysegment.intersections_composite(&shape, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {
            point: [1.0, 0.5],
            left: SegmentKey::from_segment_idx(1),
            right: SegmentKey::new(0, 1)
        })
    );
    ```
     */
    fn intersections_shape<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a;

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns the intersections between `self` and another type implementing
    [`Composite`].

    This is a generalized interface to all specialized intersection functions.
    For example, if `self` is a [`Polysegment`], the implementation of this
    function boils down to:

    ```ignore
    impl Composite for Polysegment {
        // Implementations of the other methods ...

        fn intersections_composite<'a, T: Composite>(
            &'a self,
            other: &'a T,
            epsilon: f64,
            max_ulps: u32,
        ) -> PrimitiveIntersections {
            other.intersections_polysegment(self, epsilon, max_ulps)
        }
    }
    ```

    It is recommended to use this function instead of the specialized methods
    such as [`Composite::intersections_polysegment`] for composite intersection
    to simplify the usage of this trait.

    If eager collection of the returned [`Intersection`]s into a [`Vec`] is not
    an issue, consider using the even more generic [`Composite::intersections`]
    instead (which can deal with both [`Primitive`]s and [`Composite`]s).
     */
    fn intersections_composite<'a, T: Composite>(
        &'a self,
        other: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection> + 'a
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
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        Self: Sized;

    /**
    Returns the intersections between a [`Composite`] and any other geometric
    type.

    This method is based on
    [`GeometryRef::intersections`](crate::geometry::GeometryRef::intersections).
    It can handle any geometric type ([`Primitive`] or [`Composite`]) defined in
    this crate and can therefore be seen as a combination of
    [`Composite::intersections_primitive`] and
    [`Composite::intersections_composite`]. Its downside is that it is not lazy
    (because the aforementioned methods return different iterator types which
    need to be collected into a [`Vec<Intersection>`]). If that is an issue,
    consider using the specialized methods instead (which are lazy).

    [`Composite::intersections_par`] is a parallelized variant of this method.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let intersections = polysegment.intersections(&shape, 0.0, 0);
    assert_eq!(intersections.len(), 4);
    ```
     */
    fn intersections<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        let geo_ref: crate::geometry::GeometryRef = other.into();
        return geo_ref.intersections_composite(self, epsilon, max_ulps);
    }

    /**
    Returns the intersections between a [`Composite`] and any other geometric
    type.

    This is the parallelized variant of [`Composite::intersections`].
    See its docstring for more information.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let intersections = polysegment.intersections_par(&shape, 0.0, 0);
    assert_eq!(intersections.len(), 4);
    ```
     */
    fn intersections_par<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        let geo_ref: crate::geometry::GeometryRef = other.into();
        return geo_ref.intersections_composite_par(self, epsilon, max_ulps);
    }

    /**
    Returns whether the given point is contained within the composite or not.

    The definition of "contained" depends on the composite type:
    - [`Polysegment`]: A point is contained within it if it intersects any of
    its segments.
    - [`Contour`]: Like the [`Polysegment`] definition, but additionally a
    point is also contained if it is within the enclosed surface.
    - [`Shape`]: Like the contour definition, but a point is not considered
    contained if it is inside one of the "hole" contours.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    let shape = Shape::new(vec![contour.clone(), hole]).expect("valid inputs");

    let pt = [0.5, 0.5];

    // pt is not on the polysegment defining the contour ...
    assert!(!contour.polysegment().contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // ... but it is within the contour
    assert!(contour.contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // It is however not inside the shape, because it is inside the hole
    assert!(!shape.contains_point(pt, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    ```
     */
    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool;
}

/**
This enum defines how a [`Polysegment`] / [`Contour`] should be polygonized
(approximated by a polygon). It is the counterpart of [`SegmentPolygonizer`] and
built upon it.
 */
#[derive(Debug, Clone)]
pub enum Polygonizer {
    /**
    This variant defines a single [`SegmentPolygonizer`] for all
    [`ArcSegment`](crate::segment::ArcSegment)s and another single
    [`SegmentPolygonizer`] for all [`LineSegment`](crate::segment::LineSegment)s
    of a [`Polysegment`] / [`Contour`].
     */
    PerType {
        /// [`SegmentPolygonizer`] for all arc segments.
        arc: SegmentPolygonizer,
        /// [`SegmentPolygonizer`] for all line segments.
        straight: SegmentPolygonizer,
    },
    /**
    This variant allows specifying individual [`SegmentPolygonizer`]s for each
    segment of a [`Polysegment`] / [`Contour`] via a
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
    polysegment: &'a Polysegment,
    index: usize,
    skip_last_vertex: bool,
    polygonizer: Polygonizer,
    sub_iterator: Option<crate::segment::PolygonPointsIterator<'a>>,
}

impl<'a> PointIterator<'a> {
    pub(crate) fn new(
        polysegment: &'a Polysegment,
        skip_last_vertex: bool,
        polygonizer: Polygonizer,
    ) -> Self {
        let sub_iterator = polysegment.front().map(|first_segment| {
            let sub_polygonizer = polygonizer.segment_polygonizer(first_segment, 0);

            first_segment.polygonize(sub_polygonizer)
        });
        return Self {
            polysegment,
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

                let segment = match self.polysegment.get(self.index) {
                    Some(s) => s,
                    None => {
                        self.sub_iterator = None;
                        return None;
                    }
                };

                let sub_polygonizer = self.polygonizer.segment_polygonizer(segment, self.index);

                let mut sub_iterator = segment.polygonize(sub_polygonizer);
                sub_iterator.start_at_second_point();

                // If the final vertex of the polysegment should be skipped, adjust the new
                // sub iterator here
                if self.index + 1 == self.polysegment.len() && self.skip_last_vertex {
                    sub_iterator.skip_last_vertex();
                }

                self.sub_iterator = Some(sub_iterator);
                return self.next();
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let length = self.polysegment.len();
        (length, None)
    }
}
