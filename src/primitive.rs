/*!
This module contains the [`Primitive`] trait and the [`PrimitiveIntersections`]
enum.

The geometric types in this crate are either "primitives" (which implement the
[`Primitive`] trait) or "composites" (which are made up from multiple primitives
and implement the [`Composite`](crate::composite::Composite) trait). Primitives
are simple types such as points (`[f64; 2]`), [`Line`]s,
[`LineSegment`]s or [`ArcSegment`]s. The [`Primitive`] trait provides a common
interface for intersection calculation between primitives which always return
a [`PrimitiveIntersections`]. The intersection algorithms between the composite
types are then based upon these.
 */

use crate::{
    geometry::GeometryRef,
    line::Line,
    segment::{ArcSegment, LineSegment, Segment, SegmentRef},
};
use approx::UlpsEq;
use bounding_box::BoundingBox;

/**
Result of an intersection calculation between two [`Primitive`]s.

This enum is usually created by one of the intersection functions of the
[`Primitive`] trait. It offers the following trait implementations:
- The [`std::ops::Index`] trait for treating this enum as an
array of length 0, 1 or 2 (depending on the variant).
- The [`IntoIterator`] trait to iterate over all intersection points `[f64; 2]`.
- The [`PartialEq`] trait to compare intersections. The order
of elements in the [`PrimitiveIntersections::Two`] variant does not matter:
```
use planar_geo::primitive::PrimitiveIntersections;

assert_eq!(
    PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 1.0]]),
    PrimitiveIntersections::Two([[1.0, 1.0], [0.0, 0.0]])
)
```
- The [`AbsDiffEq`](approx::AbsDiffEq), [`RelativeEq`](approx::RelativeEq) and
[`UlpsEq`] traits from the [approx] crate to allow for approximate equality
comparison. This should generally be preferred over comparing for exact equality
via the [`PartialEq`] trait when dealing with floating point values.
 */
#[derive(Clone, Copy, Debug)]
pub enum PrimitiveIntersections {
    /// The primitives don't intersect at all.
    Zero,
    /// The primitives have a single intersection point.
    One([f64; 2]),
    /// The primitives intersect in two points.
    Two([[f64; 2]; 2]),
}

impl PartialEq for PrimitiveIntersections {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::One(l), Self::One(r)) => l == r,
            (Self::Two(l), Self::Two(r)) => {
                (l[0] == r[0] && l[1] == r[1]) || (l[0] == r[1] && l[1] == r[0])
            }
            _ => core::mem::discriminant(self) == core::mem::discriminant(other),
        }
    }
}

impl IntoIterator for PrimitiveIntersections {
    type Item = [f64; 2];
    type IntoIter = PrimitiveIntersectionsIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        let intersections = match self {
            PrimitiveIntersections::Zero => [None, None],
            PrimitiveIntersections::One(i) => [Some(i), None],
            PrimitiveIntersections::Two(ai) => [Some(ai[0]), Some(ai[1])],
        };
        return PrimitiveIntersectionsIntoIter {
            intersections,
            counter: 0,
        };
    }
}

impl approx::AbsDiffEq for PrimitiveIntersections {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        match (self, other) {
            (Self::One(l), Self::One(r)) => l.abs_diff_eq(r, epsilon),
            (Self::Two(l), Self::Two(r)) => {
                (l[0].abs_diff_eq(&r[0], epsilon) && l[1].abs_diff_eq(&r[1], epsilon))
                    || (l[0].abs_diff_eq(&r[1], epsilon) && l[1].abs_diff_eq(&r[0], epsilon))
            }
            _ => core::mem::discriminant(self) == core::mem::discriminant(other),
        }
    }
}

impl approx::RelativeEq for PrimitiveIntersections {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        match (self, other) {
            (Self::One(l), Self::One(r)) => l.relative_eq(r, epsilon, max_relative),
            (Self::Two(l), Self::Two(r)) => {
                (l[0].relative_eq(&r[0], epsilon, max_relative)
                    && l[1].relative_eq(&r[1], epsilon, max_relative))
                    || (l[0].relative_eq(&r[1], epsilon, max_relative)
                        && l[1].relative_eq(&r[0], epsilon, max_relative))
            }
            _ => core::mem::discriminant(self) == core::mem::discriminant(other),
        }
    }
}

impl approx::UlpsEq for PrimitiveIntersections {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        match (self, other) {
            (Self::One(l), Self::One(r)) => l.ulps_eq(r, epsilon, max_ulps),
            (Self::Two(l), Self::Two(r)) => {
                (l[0].ulps_eq(&r[0], epsilon, max_ulps) && l[1].ulps_eq(&r[1], epsilon, max_ulps))
                    || (l[0].ulps_eq(&r[1], epsilon, max_ulps)
                        && l[1].ulps_eq(&r[0], epsilon, max_ulps))
            }
            _ => core::mem::discriminant(self) == core::mem::discriminant(other),
        }
    }
}

/**
An iterator over the intersection points in a [`PrimitiveIntersections`] enum.
It is created via the [`IntoIterator::into_iter`] implementation of
[`PrimitiveIntersections`].
 */
pub struct PrimitiveIntersectionsIntoIter {
    intersections: [Option<[f64; 2]>; 2],
    counter: usize,
}

impl Iterator for PrimitiveIntersectionsIntoIter {
    type Item = [f64; 2];

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.intersections.get_mut(self.counter)?.take()?;
        self.counter += 1;
        return Some(i);
    }
}

impl std::ops::Index<usize> for PrimitiveIntersections {
    type Output = [f64; 2];

    fn index(&self, index: usize) -> &Self::Output {
        self.get(index).expect("index is out of bounds")
    }
}

impl PrimitiveIntersections {
    /**
    Returns an iterator over references of all intersections of `self`.

    # Examples

    ```
    use planar_geo::primitive::PrimitiveIntersections;

    assert_eq!(PrimitiveIntersections::Zero.iter().count(), 0);
    assert_eq!(PrimitiveIntersections::One([0.0, 0.0]).iter().count(), 1);
    assert_eq!(PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]]).iter().count(), 2);
    ```
     */
    pub fn iter(&self) -> PrimitiveIntersectionsIter<'_> {
        return PrimitiveIntersectionsIter {
            intersections: self,
            counter: 0,
        };
    }

    /**
    Returns a reference to the intersection at the given index:
    - If `self` is [`PrimitiveIntersections::Zero`], any index returns [`None`].
    - If `self` is [`PrimitiveIntersections::One`], index 0 returns the contained
    intersection, all other indices return [`None`].
    - If `self` is [`PrimitiveIntersections::Two`], indices 0 and 1 return the respective
    intersection, all other indices return [`None`].

    # Examples

    ```
    use planar_geo::primitive::PrimitiveIntersections;

    let zero = PrimitiveIntersections::Zero;
    let one = PrimitiveIntersections::One([0.0, 0.0]);
    let two = PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]]);

    assert_eq!(zero.get(0), None);

    assert_eq!(one.get(0), Some(&[0.0, 0.0]));
    assert_eq!(one.get(1), None);

    assert_eq!(two.get(0), Some(&[0.0, 0.0]));
    assert_eq!(two.get(1), Some(&[1.0, 0.0]));
    assert_eq!(two.get(2), None);
    ```
     */
    pub fn get(&self, index: usize) -> Option<&[f64; 2]> {
        match self {
            PrimitiveIntersections::Zero => None,
            PrimitiveIntersections::One(i) => {
                if index == 0 {
                    return Some(i);
                } else {
                    return None;
                }
            }
            PrimitiveIntersections::Two(ai) => ai.get(index),
        }
    }

    /**
    Like [`PrimitiveIntersections::get`], but returns a mutable reference.
     */
    pub fn get_mut(&mut self, index: usize) -> Option<&mut [f64; 2]> {
        match self {
            PrimitiveIntersections::Zero => None,
            PrimitiveIntersections::One(i) => {
                if index == 0 {
                    return Some(i);
                } else {
                    return None;
                }
            }
            PrimitiveIntersections::Two(ai) => ai.get_mut(index),
        }
    }

    /**
    Returns the number of intersections stored in `self`.

    # Examples

    ```
    use planar_geo::primitive::PrimitiveIntersections;

    assert_eq!(PrimitiveIntersections::Zero.len(), 0);
    assert_eq!(PrimitiveIntersections::One([0.0, 0.0]).len(), 1);
    assert_eq!(PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]]).len(), 2);
    ```
     */
    pub fn len(&self) -> usize {
        match self {
            PrimitiveIntersections::Zero => 0,
            PrimitiveIntersections::One(_) => 1,
            PrimitiveIntersections::Two(_) => 2,
        }
    }

    /**
    Tries to add another intersection to [`PrimitiveIntersections`].

    If `self` is [`PrimitiveIntersections::Two`], no further intersection
    can be added and `false` is returned. Otherwise, `true` is returned
     */
    pub(crate) fn push(&mut self, intersection: [f64; 2]) -> bool {
        match self {
            PrimitiveIntersections::Zero => {
                *self = PrimitiveIntersections::One(intersection);
                return true;
            }
            PrimitiveIntersections::One(i) => {
                *self = PrimitiveIntersections::Two([*i, intersection]);
                return true;
            }
            PrimitiveIntersections::Two(_) => {
                return false;
            }
        }
    }
}

/**
An iterator over the intersection points in a [`PrimitiveIntersections`] enum.
It is created via the [`PrimitiveIntersections::iter`] method.
 */
pub struct PrimitiveIntersectionsIter<'a> {
    intersections: &'a PrimitiveIntersections,
    counter: usize,
}

impl<'a> Iterator for PrimitiveIntersectionsIter<'a> {
    type Item = &'a [f64; 2];

    fn next(&mut self) -> Option<Self::Item> {
        match self.intersections {
            PrimitiveIntersections::Zero => return None,
            PrimitiveIntersections::One(i) => {
                if self.counter == 0 {
                    self.counter += 1;
                    Some(i)
                } else {
                    return None;
                }
            }
            PrimitiveIntersections::Two(ai) => {
                let i = match self.counter {
                    0 => &ai[0],
                    1 => &ai[1],
                    _ => return None,
                };
                self.counter += 1;
                return Some(i);
            }
        }
    }
}

pub(crate) mod private {
    /// Sealed trait for [`Primitive`]
    pub trait Sealed {}
}

/**
A trait for "primitive" geometric types like points (`[f64; 2]`), [`Line`]s,
[`LineSegment`]s and [`ArcSegment`]s.

This trait provides methods for coverage checking and intersection calculation
between "primitive" geometric types representing one- or zero-dimensional
objects. It is not meant to be implemented for other types, hence it is
[sealed](https://rust-lang.github.io/api-guidelines/future-proofing.html).

# Coverage

A primitive "covers" another one if all points of the latter are part of the
former. The generic [`Primitive::covers`] method provides a unified interface to
check this. Specialized `covers_` methods for specific types are available as
well and are particularily useful when using [`Primitive`] to define a trait
object (which cannot use generic functions).

By definition, a primitive always covers itself.

```
use planar_geo::prelude::*;

let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
assert!(ls.covers(&[0.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
assert!(!ls.covers(&[0.5, 0.1], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
assert!(ls.covers(&ls, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
```

# Intersection

Two primitives "intersect" if they share at least one point. If a primitive
covers another one, only the start and end points of the covered primitive
(if available) are reported as intersections. The generic
[`Primitive::intersections_primitive`] method provides a unified interface to
find all intersection points between two primitives. However, specialized
intersection methods forspecific types are available as well
and are particularily useful when using [`Primitive`] to define a trait object
(which cannot use generic functions).

All intersection functions first check if the bounding boxes of the two
primitives overlap (short-circuiting the evaluation if they don't). Hence, it is
not necessary to check this before calling an intersection method.

By definition, a primitive does not self-intersect, but it does intersect with
equal primitives:

```
use planar_geo::prelude::*;

let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

// Self-intersection
assert_eq!(
    ls.intersections_primitive(&ls, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
    PrimitiveIntersections::Zero
);

// Intersections with equal primitive
let ls_cloned = ls.clone();
assert_eq!(
    ls.intersections_primitive(&ls_cloned, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
    PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]])
);
```
 */
pub trait Primitive: private::Sealed + Sync {
    /**
    Returns `true` if `self` covers the given point.

    Since floating point values are prone to rounding errors and precision
    issues (i.e. the number 0.2 cannot be represented exactly as a floating
    point number at all), the implementations for the different primitive types
    check whether the point is covered within a tolerance band defined by
    `epsilon` and `max_ulps`. See the [crate-level documentation](crate).

    # Examples

    ```
    use planar_geo::prelude::*;

    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let p1 = [0.5, -0.5];
    let p2 = [0.5, 0.5];
    let p3 = [1.5, 0.5];

    // Default tolerances for "intuitive" behaviour
    assert!(ls.covers_point(p1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    assert!(!ls.covers_point(p2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    assert!(!ls.covers_point(p3, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));

    // By increasing the absolute tolerance `epsilon`, the behaviour of the
    // function is changed
    assert!(ls.covers_point(p1, 10.0, DEFAULT_MAX_ULPS));
    assert!(ls.covers_point(p2, 10.0, DEFAULT_MAX_ULPS));
    assert!(ls.covers_point(p3, 10.0, DEFAULT_MAX_ULPS));
    ```
     */
    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool;

    /**
    Returns if `self` covers the given [`ArcSegment`].

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let s1: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m).unwrap().into();
    let quarter = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], 0.5*PI, e, m).unwrap();
    assert!(s1.covers_arc_segment(&quarter, e, m));

    let s2: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m).unwrap().into();
    assert!(!s2.covers_arc_segment(&quarter, e, m));
    ```
     */
    fn covers_arc_segment(&self, arc_segment: &ArcSegment, epsilon: f64, max_ulps: u32) -> bool;

    /**
    Returns if `self` covers the given [`LineSegment`].

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let s1: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m).unwrap().into();
    let ls_start_to_middle = LineSegment::new([0.0, 0.0], [0.5, 0.5], e, m).unwrap();
    assert!(s1.covers_line_segment(&ls_start_to_middle, e, m));

    let s2: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m).unwrap().into();
    assert!(!s2.covers_line_segment(&ls_start_to_middle, e, m));
    ```
     */
    fn covers_line_segment(&self, line_segment: &LineSegment, epsilon: f64, max_ulps: u32) -> bool;

    /**
    Returns if `self` covers the given [`Line`].

    Since a [`Line`] has infinite length, it can only be covered by another
    [`Line`] which is identical. For any other [`Primitive`], this function
    returns false.

    # Examples

    ```
    use planar_geo::prelude::*;

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let l1 = Line::from_point_angle([0.0, 0.0], 0.0);
    let l2 = Line::from_point_angle([1.0, 0.0], 0.0);
    let l3 = Line::from_point_angle([0.0, 0.0], 1.0);

    // l1 and l2 are identical
    assert!(l1.covers_line(&l2, e, m));
    assert!(l2.covers_line(&l1, e, m));

    // l3 is different
    assert!(!l1.covers_line(&l3, e, m));
    ```
    */
    fn covers_line(&self, line: &Line, epsilon: f64, max_ulps: u32) -> bool;

    /**
    Returns if `self` covers `other` (which can be any geometric type).

    Internally, this function converts `self` and `other` to [`GeometryRef`]s
    and then matches them to select the specific `covers_` function for the
    pairing (e.g. [`Primitive::covers_line`] if `other` is a [`Line`]). To
    avoid matching for maximum performance, consider using the specific method
    directly.
     */
    fn covers<'a, T>(&self, other: T, epsilon: f64, max_ulps: u32) -> bool
    where
        Self: Sized,
        for<'b> &'b Self: Into<GeometryRef<'b>>,
        T: Into<GeometryRef<'a>>,
    {
        let geo_self: GeometryRef = self.into();
        let geo_other: GeometryRef = other.into();
        match geo_self {
            GeometryRef::Point(this) => {
                if let GeometryRef::Point(pt_other) = geo_other {
                    return this.covers_point(*pt_other, epsilon, max_ulps);
                }
                return false;
            }
            GeometryRef::BoundingBox(this) => {
                let bb_other = BoundingBox::from(&geo_other);
                return this.approx_covers(&bb_other, epsilon, max_ulps);
            }
            GeometryRef::ArcSegment(this) => match geo_other {
                GeometryRef::Point(o) => {
                    return this.covers_point(*o, epsilon, max_ulps);
                }
                GeometryRef::ArcSegment(o) => {
                    return this.covers_arc_segment(o, epsilon, max_ulps);
                }
                GeometryRef::Segment(segment) => match segment {
                    Segment::LineSegment(_) => return false,
                    Segment::ArcSegment(o) => {
                        return this.covers_arc_segment(o, epsilon, max_ulps);
                    }
                },
                _ => return false,
            },
            GeometryRef::LineSegment(this) => match geo_other {
                GeometryRef::Point(o) => {
                    return this.covers_point(*o, epsilon, max_ulps);
                }
                GeometryRef::LineSegment(o) => {
                    return this.covers_line_segment(o, epsilon, max_ulps);
                }
                GeometryRef::Segment(segment) => match segment {
                    Segment::LineSegment(o) => {
                        return this.covers_line_segment(o, epsilon, max_ulps);
                    }
                    Segment::ArcSegment(_) => return false,
                },
                _ => return false,
            },
            GeometryRef::Line(this) => match geo_other {
                GeometryRef::Point(o) => return this.covers_point(*o, epsilon, max_ulps),
                GeometryRef::Line(o) => {
                    return this.covers_line(o, epsilon, max_ulps);
                }
                _ => return false,
            },
            GeometryRef::Segment(segment) => match segment {
                Segment::LineSegment(this) => match geo_other {
                    GeometryRef::Point(o) => {
                        return this.covers_point(*o, epsilon, max_ulps);
                    }
                    GeometryRef::LineSegment(ls_other) => {
                        return this.covers_line_segment(ls_other, epsilon, max_ulps);
                    }
                    GeometryRef::Segment(segment) => match segment {
                        Segment::LineSegment(ls_other) => {
                            return this.covers_line_segment(ls_other, epsilon, max_ulps);
                        }
                        Segment::ArcSegment(_) => return false,
                    },
                    _ => return false,
                },
                Segment::ArcSegment(this) => match geo_other {
                    GeometryRef::Point(o) => {
                        return this.covers_point(*o, epsilon, max_ulps);
                    }
                    GeometryRef::ArcSegment(as_other) => {
                        return this.covers_arc_segment(as_other, epsilon, max_ulps);
                    }
                    GeometryRef::Segment(segment) => match segment {
                        Segment::LineSegment(_) => return false,
                        Segment::ArcSegment(as_other) => {
                            return this.covers_arc_segment(as_other, epsilon, max_ulps);
                        }
                    },
                    _ => return false,
                },
            },
            _ => false, // Won't be reached anyway, since self is a primitive
        }
    }

    /**
    Returns the intersections between `self` and a point `[f64; 2]`

    This function wraps the given point in [`PrimitiveIntersections::One`] if
    [`Primitive::covers_point`] returned `true` and returns
    [`PrimitiveIntersections::Zero`] otherwise.

    This method is mainly used to implement
    [`Primitive::intersections_primitive`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let ls = LineSegment::new([0.0, 0.0], [1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let p1 = [0.5, -0.5];
    let p2 = [0.5, 0.5];
    let p3 = [1.5, 0.5];

    // Default tolerances for "intuitive" behaviour
    assert_eq!(ls.intersections_point(p1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS), PrimitiveIntersections::One(p1));
    assert_eq!(ls.intersections_point(p2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS), PrimitiveIntersections::Zero);
    assert_eq!(ls.intersections_point(p3, DEFAULT_EPSILON, DEFAULT_MAX_ULPS), PrimitiveIntersections::Zero);

    // By increasing the absolute tolerance `epsilon`, the behaviour of the
    // function is changed
    assert_eq!(ls.intersections_point(p1, 10.0, DEFAULT_MAX_ULPS), PrimitiveIntersections::One(p1));
    assert_eq!(ls.intersections_point(p2, 10.0, DEFAULT_MAX_ULPS), PrimitiveIntersections::One(p2));
    assert_eq!(ls.intersections_point(p3, 10.0, DEFAULT_MAX_ULPS), PrimitiveIntersections::One(p3));
    ```
     */
    fn intersections_point(
        &self,
        point: [f64; 2],
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        if self.covers_point(point, epsilon, max_ulps) {
            return PrimitiveIntersections::One(point);
        }
        return PrimitiveIntersections::Zero;
    }

    /**
    Returns the intersections between `self` and a [`Line`].

    This is a straightforward calculation, but there are two degenerate cases
    which need to be treated explicitly:
    - Intersection between two [`Line`]s which are identical
    - Intersection between a [`Line`] and [`LineSegment`] where the latter
    is contained in the former.
    In both cases, the returned result is [`PrimitiveIntersections::Zero`] (see
    examples).

    This method is mainly used to implement
    [`Primitive::intersections_primitive`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let line_1 = Line::from_point_angle([0.0, 0.0], 1.0);
    let line_2 = Line::from_point_angle([0.0, 0.0], 0.0);
    let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    // Non-degenerate cases
    assert_eq!(
        line_2.intersections_line(&line_1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::One([0.0, 0.0])
    );
    assert_eq!(
        ls.intersections_line(&line_1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::One([0.0, 0.0])
    );

    // Degenerate cases
    assert_eq!(
        line_1.intersections_line(&line_1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );
    assert_eq!(
        ls.intersections_line(&line_2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );
    ```
     */
    fn intersections_line(
        &self,
        line: &Line,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections;

    /**
    Returns the intersections between `self` and a [`LineSegment`].

    This is a straightforward calculation, but there are two degenerate cases
    which need to be treated explicitly:
    - Intersection between a [`Line`] and [`LineSegment`] where the latter
    is contained in the former.
    - Intersection between two [`LineSegment`]s which overlap.
    In the former case, the result is defined to be
    [`PrimitiveIntersections::Zero`]. In the latter case, the "common endpoints"
    are returned (see examples).

    This method is mainly used to implement
    [`Primitive::intersections_primitive`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let ls1 = LineSegment::new([0.0, 0.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let ls2 = LineSegment::new([0.5, 0.5], [0.5, -0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let ls3 = LineSegment::new([1.0, 0.0], [1.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let ls4 = LineSegment::new([1.0, 0.0], [3.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    // Non-degenerate case
    assert_eq!(
        ls1.intersections_line_segment(&ls2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::One([0.5, 0.0])
    );

    // Degenerate cases
    assert_eq!(
        ls1.intersections_line_segment(&ls3, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[1.0, 0.0], [1.5, 0.0]])
    );
    assert_eq!(
        ls1.intersections_line_segment(&ls4, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[1.0, 0.0], [2.0, 0.0]])
    );
    ```
     */
    fn intersections_line_segment(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections;

    /**
    Returns the intersections between `self` and an [`ArcSegment`].

    Similar to the intersection of two [`LineSegment`]s which overlap (see
    discussion in the docstring of [`Primitive::intersections_line_segment`]),
    this function has to deal with the case of two overlapping [`ArcSegment`]s.
    As with the [`LineSegment`] case, the "common endpoints" are returned (see
    examples).

    This method is mainly used to implement
    [`Primitive::intersections_primitive`], consider using that method instead.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let arc1 = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let line = Line::from_point_angle([0.0, 0.0], 0.5 * PI);
    let arc2 = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.5*PI, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    // Regular intersection
    approx::assert_abs_diff_eq!(
        line.intersections_arc_segment(&arc1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::One([0.0, 2.0]), epsilon = DEFAULT_EPSILON
    );

    // Degenerate case
    approx::assert_abs_diff_eq!(
        arc1.intersections_arc_segment(&arc2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[0.0, 2.0], [-2.0, 0.0]]), epsilon = DEFAULT_EPSILON
    );
    ```
    */
    fn intersections_arc_segment(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections;

    /**
    Returns the intersection between `self` and a [`Segment`].

    Since [`Segment`] is an enum of either [`LineSegment`] or [`ArcSegment`],
    this function behaves in the same way as
     [`Primitive::intersections_line_segment`] or
    [`Primitive::intersections_arc_segment`] respectively (especially for
    degenerate cases, see the function docstrings).

    This method is mainly used to implement
    [`Primitive::intersections_primitive`], consider using that method instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let ls: Segment = LineSegment::new([0.0, 0.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();
    let arc: Segment = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, -1.0, 1.0, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap().into();

    assert_eq!(
        ls.intersections_segment(&arc, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::One([1.0, 0.0])
    );
    ```
     */
    fn intersections_segment<'a, T>(
        &self,
        segment: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections
    where
        T: Into<SegmentRef<'a>>,
    {
        let segment: SegmentRef = segment.into();
        match segment {
            SegmentRef::LineSegment(line_segment) => {
                self.intersections_line_segment(line_segment, epsilon, max_ulps)
            }
            SegmentRef::ArcSegment(arc_segment) => {
                self.intersections_arc_segment(arc_segment, epsilon, max_ulps)
            }
        }
    }

    /**
    Returns the intersection between `self` and another type implementing
    [`Primitive`].

    This is a generalized interface to all the special intersection functions.
    For example, if `other` is a [`LineSegment`], the implementation of this
    function boils down to:

    ```ignore
    impl Primitive for LineSegment {
        // Implementations of the other methods ...

        fn intersections_primitive<T: Primitive>(
            &self,
            other: &T,
            epsilon: f64,
            max_ulps: u32,
        ) -> PrimitiveIntersections {
            other.intersections_line_segment(self, epsilon, max_ulps)
        }
    }
    ```

    It is recommended to use this function instead of the specialized methods
    such as [`Primitive::intersections_line_segment`] for primitive intersection
    to simplify the usage of this trait.
     */
    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections
    where
        Self: Sized;

    /**
    Returns the intersections between a [`Primitive`] and any other geometric
    type.

    This method is based on
    [`GeometryRef::intersections`](GeometryRef::intersections).
    It can handle any geometric type ([`Primitive`] or
    [`Composite`](crate::composite::Composite)) defined in this crate, but it
    needs to allocate a vector for the results. For intersections between two
    primitives, [`Primitive::intersections_primitive`] is therefore to be
     preferred. If `other` is a [`Composite`](crate::composite::Composite), the
    [`Composite::intersections_primitive`](crate::composite::Composite::intersections_primitive)
    method offers a non-allocating (lazy) alternative.

    The main advantage of this method is its genericness. If allocating a vector
    for the results is not an issue, consider using this method instead of the
    specialized variants mentioned above to simplify the interface to this
    trait.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let arc1 = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let arc2 = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.5*PI, PI, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    let intersections = arc1.intersections(&arc2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    assert_eq!(intersections.len(), 2);
    approx::assert_abs_diff_eq!(intersections[0].point, [-2.0, 0.0], epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(intersections[1].point, [0.0, 2.0], epsilon = DEFAULT_EPSILON);

    // Comparison to specialized function:
    approx::assert_abs_diff_eq!(
        arc1.intersections_primitive(&arc2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[0.0, 2.0], [-2.0, 0.0]]), epsilon = DEFAULT_EPSILON
    );
    ```
     */
    fn intersections<'a, T>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<crate::composite::Intersection>
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>,
    {
        let geo_ref: GeometryRef = other.into();
        return geo_ref.intersections_primitive(self, epsilon, max_ulps);
    }
}

impl private::Sealed for [f64; 2] {}

impl Primitive for [f64; 2] {
    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        return self.ulps_eq(&point, epsilon, max_ulps);
    }

    fn covers_arc_segment(&self, _arc_segment: &ArcSegment, _epsilon: f64, _max_ulps: u32) -> bool {
        return false;
    }

    fn covers_line_segment(
        &self,
        _line_segment: &LineSegment,
        _epsilon: f64,
        _max_ulps: u32,
    ) -> bool {
        return false;
    }

    fn covers_line(&self, _line: &Line, _epsilon: f64, _max_ulps: u32) -> bool {
        return false;
    }

    fn intersections_line(
        &self,
        line: &Line,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        line.intersections_point(*self, epsilon, max_ulps)
    }

    fn intersections_line_segment(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        line_segment.intersections_point(*self, epsilon, max_ulps)
    }

    fn intersections_arc_segment(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        arc_segment.intersections_point(*self, epsilon, max_ulps)
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        other.intersections_point(*self, epsilon, max_ulps)
    }
}
