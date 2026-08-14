/*!
This module contains the [`Primitive`] trait and the [`PrimitiveIntersections`]
enum.

The geometric types in this crate are either "primitives" (which implement the
[`Primitive`] trait) or "composites" (which are made up from multiple primitives
and implement the [`Composite`](crate::composite::Composite) trait). Primitives
are simple types such as points (`[f64; 2]`), [`Line`]s,
[`LineSegment`]s or [`ArcSegment`]s. The [`Primitive`] trait provides a common
interface for intersection and coverage calculation.
 */

use crate::{
    DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE, ToleranceContext, WithTolerance,
    geometry::GeometryRef,
    line::Line,
    segment::{ArcSegment, LineSegment, Segment, SegmentRef},
};
use approx::RelativeEq;

/**
Result of an intersection calculation between two [`Primitive`]s.

By definition, any two [`Primitive`]s can either have zero, one or two
intersections, see [`Primitive::intersections_primitive`]. This enum covers all
these cases in its variants.

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
- The [`AbsDiffEq`](approx::AbsDiffEq), [`RelativeEq`] and
[`UlpsEq`](approx::UlpsEq) traits from the [approx] crate to allow for
approximate equality comparison. This should generally be preferred over
comparing for exact equality via the [`PartialEq`] trait when dealing with
floating point values.
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

    fn ulps_eq(&self, other: &Self, epsilon: f64, ulps_eq: u32) -> bool {
        match (self, other) {
            (Self::One(l), Self::One(r)) => l.ulps_eq(r, epsilon, ulps_eq),
            (Self::Two(l), Self::Two(r)) => {
                (l[0].ulps_eq(&r[0], epsilon, ulps_eq) && l[1].ulps_eq(&r[1], epsilon, ulps_eq))
                    || (l[0].ulps_eq(&r[1], epsilon, ulps_eq)
                        && l[1].ulps_eq(&r[0], epsilon, ulps_eq))
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
    - If `self` is [`PrimitiveIntersections::One`], index 0 returns the
    contained intersection, all other indices return [`None`].
    - If `self` is [`PrimitiveIntersections::Two`], indices 0 and 1 return the
    respective intersection, all other indices return [`None`].

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
            PrimitiveIntersections::Two(i) => i.get(index),
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

/**
A trait for "primitive" geometric types like points (`[f64; 2]`), [`Line`]s,
[`LineSegment`]s and [`ArcSegment`]s.

This trait provides methods for coverage checking and intersection calculation
between "primitive" geometric types representing one- or zero-dimensional
objects. It is not meant to be implemented for other types, hence it is
[sealed](https://rust-lang.github.io/api-guidelines/future-proofing.html).

# Coverage

A [`Primitive`] [`covers`](Primitive::covers) another one if all points of the
latter are part of the former. By definition, a primitive always covers itself.

```
use planar_geo::prelude::*;

let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();;
assert!(ls.covers(&[0.5, 0.0]));
assert!(!ls.covers(&[0.5, 0.1]));
assert!(ls.covers(&ls));
```

# Intersection

Two geometric entities "intersect" if they share at least one point. The
[`Primitive`] provides two methods for finding intersection points:
[`Primitive::intersections_primitive`] and [`Primitive::intersections`].

As the name suggests, [`Primitive::intersections_primitive`] can only be used if
both entities implement the [`Primitive`] trait. By definition, two primitives
can either have zero, one or two intersections by definition, therefore
`intersections_primitive` returns the stack-allocated [`PrimitiveIntersections`]
enum. If applicable, this method is to be preferred over
[`Primitive::intersections`].

[`Primitive::intersections`] is the generic fallback which can handle a
[`Composite`](crate::composite::Composite) as second argument (besides
[`Primitive`]s). It returns all found intersections in a [`Vec`]. If a lazy
iterator is needed instead, using one of the specialized `intersections_`
methods from the [`Composite`](crate::composite::Composite) trait.

Intersection calculations perform inexpensive preliminary checks to reject
geometries that cannot intersect. Callers therefore generally do not need to
perform their own geometric prechecks before calculating intersections.
 */
pub trait Primitive: crate::private::Sealed + Sync {
    /**
    Returns whether every point of `other` lies within or on the boundary of
    `self`. By definition, a primitive always covers itself.

    By default, [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_RELATIVE`] are used for
    floating-point comparisons. For custom tolerances, use
    [`WithTolerance::with_tolerance`].

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let ls1 = LineSegment::new([0.0, 0.0], [3.0, 0.0]).expect("valid inputs");
    let line = Line::from(&ls1);
    let ls2 = LineSegment::new([1.0, 0.0], [2.0, 0.0]).expect("valid inputs");
    let ls3 = LineSegment::new([1.0, 1.0], [1.0, -1.0]).expect("valid inputs");
    let arc_segment = ArcSegment::from_center_radius_start_sweep_angle([1.0, 0.0], 1.0, 0.0, PI).expect("valid inputs");

    assert!(ls1.covers(&ls2));
    assert!(!ls2.covers(&ls1));
    assert!(!ls1.covers(&ls3));

    // An infinitely long line can cover a finite line segment
    assert!(line.covers(&ls1));
    assert!(line.covers(&ls2));

    // A finite line segment cannot cover an infinite line
    assert!(!ls1.covers(&line));
    assert!(!ls2.covers(&line));

    // Arc and line segment cannot cover each other
    assert!(!arc_segment.covers(&ls1));
    assert!(!ls1.covers(&arc_segment));

    // ls1 does not cover the point [0.0, 0.1] with default tolerances, but
    // with sufficiently large custom tolerances, it does
    assert!(!ls1.covers(&[0.0, 0.1]));
    assert!(ls1.with_tolerance(0.2, 0.0).covers(&[0.0, 0.1]));

    // Primitives cover themselves
    assert!(ls1.covers(&ls1));
    ```
     */
    fn covers<'a, T>(&self, other: T) -> bool
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>;

    /**
    Returns the intersections between two [`Primitive`]s in the stack-allocated
    [`PrimitiveIntersections`] enum.

    Two primitives can have at most two distinct intersection points. If one
    primitive [`covers`](Primitive::covers) another, their intersection is
    represented by the start and end points of the covered primitive (if
    available). By definition, a primitive does not self-intersect, but a
    primitive does intersect with an equal primitive.

    [`PrimitiveIntersections`] is stack-allocated and can represent up to two
    intersection points without heap allocation.

    By default, [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_RELATIVE`] are used for
    floating-point comparisons. For custom tolerances, use
    [`WithTolerance::with_tolerance`].

    # Examples
    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let ls1 = LineSegment::new([0.0, 0.0], [3.0, 0.0]).expect("valid inputs");
    let ls2 = LineSegment::new([1.0, 1.0], [1.0, -1.0]).expect("valid inputs");
    let ls3 = LineSegment::new([1.0, 0.0], [4.0, 0.0]).expect("valid inputs");
    let arc_segment = ArcSegment::from_center_radius_start_sweep_angle([1.0, 0.0], 1.0, 0.0, PI).expect("valid inputs");

    // Intersection between two perpendicular line segments
    assert_eq!(
        ls1.intersections_primitive(&ls2),
        PrimitiveIntersections::One([1.0, 0.0])
    );

    // Intersection between line and arc segments
    assert_eq!(
        ls1.intersections_primitive(&arc_segment),
        PrimitiveIntersections::Two([[0.0, 0.0], [2.0, 0.0]])
    );
    assert_eq!(
        arc_segment.intersections_primitive(&ls1),
        PrimitiveIntersections::Two([[0.0, 0.0], [2.0, 0.0]])
    );

    // ls1 and the point [0.0, 0.1] don't intersect with default tolerances, but
    // with sufficiently large custom tolerances, they do
    assert_eq!(ls1.intersections_primitive(&[0.0, 0.1]), PrimitiveIntersections::Zero);
    assert_eq!(
        ls1.with_tolerance(0.2, 0.0).intersections_primitive(&[0.0, 0.1]),
        PrimitiveIntersections::One([0.0, 0.1])
    );

    // If two segments overlap, the start and end point of the common parts are reported
    assert_eq!(
        ls1.intersections_primitive(&ls3),
        PrimitiveIntersections::Two([[3.0, 0.0], [1.0, 0.0]])
    );

    // Self-intersection
    assert_eq!(ls1.intersections_primitive(&ls1), PrimitiveIntersections::Zero);

    // Intersections with equal primitive
    let ls1_cloned = ls1.clone();
    assert_eq!(
        ls1.intersections_primitive(&ls1_cloned),
        PrimitiveIntersections::Two([[0.0, 0.0], [3.0, 0.0]])
    );
    ```
     */
    fn intersections_primitive<'a, T: Primitive>(&self, other: &'a T) -> PrimitiveIntersections
    where
        &'a T: Into<GeometryRef<'a>>,
        Self: Sized;

    /**
    Returns all intersections between `self` and another geometric entity.

    The intersections are returned in a [`Vec`]. For a lazy iterator, use one of
    the specialized `intersections_` methods provided by the
    [`Composite`](crate::composite::Composite) trait. If the second geometric
    entity is also a [`Primitive`], consider using [`intersections_primitive`](Primitive::intersections_primitive) instead.

    By definition, a geometric entity does not self-intersect, but it does
    intersect with an equal entity.

    By default, [`DEFAULT_EPSILON`] and [`DEFAULT_MAX_RELATIVE`] are used for
    floating-point comparisons. For custom tolerances, use
    [`WithTolerance::with_tolerance`].

    # Examples

    ```
    use planar_geo::prelude::*;

    let ls1 = LineSegment::new([0.0, 0.0], [3.0, 0.0]).expect("valid inputs");
    let ls2 = LineSegment::new([1.0, 1.0], [1.0, -1.0]).expect("valid inputs");

    // Same result as with intersections_primitive:
    let intersection = Intersection {
        point: [1.0, 0.0],
        left: SegmentKey::from_segment_idx(0),
        right: SegmentKey::from_segment_idx(0)
    };
    assert_eq!(ls1.intersections(&ls2), vec![intersection]);

    // Usage with a composite
    let polysegment = Polysegment::from_points(&[[1.0, 1.0], [1.0, -1.0], [2.0, -1.0], [2.0, 1.0]]);
    let i1 = Intersection {
        point: [1.0, 0.0],
        left: SegmentKey::from_segment_idx(0),
        right: SegmentKey::from_segment_idx(0)
    };
    let i2 = Intersection {
        point: [2.0, 0.0],
        left: SegmentKey::from_segment_idx(0),
        right: SegmentKey::from_segment_idx(2)
    };
    assert_eq!(ls1.intersections(&polysegment), vec![i1, i2]);

    // Using custom tolerances
    let intersection = Intersection {
        point: [0.0, 0.1],
        left: SegmentKey::from_segment_idx(0),
        right: SegmentKey::from_segment_idx(0)
    };
    assert_eq!(ls1.intersections(&[0.0, 0.1]), Vec::new());
    assert_eq!(ls1.with_tolerance(0.2, 0.0).intersections(&[0.0, 0.1]), vec![intersection]);
    ```
     */
    fn intersections<'a, T>(&self, other: T) -> Vec<crate::composite::Intersection>
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>;
}

pub(crate) trait PrimitiveWithTol: Primitive {
    fn covers_point_tol(&self, point: &[f64; 2], epsilon: f64, max_relative: f64) -> bool;

    fn covers_line_segment_tol(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> bool;

    fn covers_arc_segment_tol(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> bool;

    fn covers_line_tol(&self, line: &Line, epsilon: f64, max_relative: f64) -> bool;

    /**
    Returns if `self` covers `other` (which can be any geometric type).

    Internally, this function converts `self` and `other` to [`GeometryRef`]s
    and then matches them to select the specific `covers_` function for the
    pairing (e.g. [`Primitive::covers_line`] if `other` is a [`Line`]). To
    avoid matching for maximum performance, consider using the specific method
    directly.
     */
    fn covers_tol<'a, T>(&self, other: T, epsilon: f64, max_relative: f64) -> bool
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
                    return this.covers_point_tol(pt_other, epsilon, max_relative);
                }
                return false;
            }
            GeometryRef::BoundingBox(this) => {
                let bb_other = bounding_box::BoundingBox::from(&geo_other);
                return this.approx_covers(&bb_other, epsilon);
            }
            GeometryRef::ArcSegment(this) => match geo_other {
                GeometryRef::Point(o) => {
                    return this.covers_point_tol(o, epsilon, max_relative);
                }
                GeometryRef::ArcSegment(o) => {
                    return this.covers_arc_segment_tol(o, epsilon, max_relative);
                }
                GeometryRef::Segment(segment) => match segment {
                    Segment::LineSegment(_) => return false,
                    Segment::ArcSegment(o) => {
                        return this.covers_arc_segment_tol(o, epsilon, max_relative);
                    }
                },
                _ => return false,
            },
            GeometryRef::LineSegment(this) => match geo_other {
                GeometryRef::Point(o) => {
                    return this.covers_point_tol(o, epsilon, max_relative);
                }
                GeometryRef::LineSegment(o) => {
                    return this.covers_line_segment_tol(o, epsilon, max_relative);
                }
                GeometryRef::Segment(segment) => match segment {
                    Segment::LineSegment(o) => {
                        return this.covers_line_segment_tol(o, epsilon, max_relative);
                    }
                    Segment::ArcSegment(_) => return false,
                },
                _ => return false,
            },
            GeometryRef::Line(this) => match geo_other {
                GeometryRef::Point(o) => return this.covers_point_tol(o, epsilon, max_relative),
                GeometryRef::Line(o) => {
                    return this.covers_line_tol(o, epsilon, max_relative);
                }
                GeometryRef::LineSegment(o) => {
                    return this.covers_line_segment_tol(o, epsilon, max_relative);
                }
                GeometryRef::BoundingBox(o) => {
                    let c = crate::contour::Contour::from(o);
                    self.covers_tol(&c, epsilon, max_relative)
                }
                GeometryRef::Polysegment(o) => o.segments().all(|s| {
                    if let Segment::LineSegment(ls) = s {
                        this.covers_line_segment_tol(ls, epsilon, max_relative)
                    } else {
                        false
                    }
                }),
                GeometryRef::Contour(o) => self.covers_tol(o.polysegment(), epsilon, max_relative),
                _ => return false,
            },
            GeometryRef::Segment(segment) => match segment {
                Segment::LineSegment(this) => match geo_other {
                    GeometryRef::Point(o) => {
                        return this.covers_point_tol(o, epsilon, max_relative);
                    }
                    GeometryRef::LineSegment(ls_other) => {
                        return this.covers_line_segment_tol(ls_other, epsilon, max_relative);
                    }
                    GeometryRef::Segment(segment) => match segment {
                        Segment::LineSegment(ls_other) => {
                            return this.covers_line_segment_tol(ls_other, epsilon, max_relative);
                        }
                        Segment::ArcSegment(_) => return false,
                    },
                    _ => return false,
                },
                Segment::ArcSegment(this) => match geo_other {
                    GeometryRef::Point(o) => {
                        return this.covers_point_tol(o, epsilon, max_relative);
                    }
                    GeometryRef::ArcSegment(as_other) => {
                        return this.covers_arc_segment_tol(as_other, epsilon, max_relative);
                    }
                    GeometryRef::Segment(segment) => match segment {
                        Segment::LineSegment(_) => return false,
                        Segment::ArcSegment(as_other) => {
                            return this.covers_arc_segment_tol(as_other, epsilon, max_relative);
                        }
                    },
                    _ => return false,
                },
            },
            _ => false, // Won't be reached anyway, since self is a primitive
        }
    }

    fn intersections_point_tol(
        &self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        if self.covers_point_tol(point, epsilon, max_relative) {
            return PrimitiveIntersections::One(*point);
        }
        return PrimitiveIntersections::Zero;
    }

    fn intersections_line_tol(
        &self,
        line: &Line,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections;

    fn intersections_line_segment_tol(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections;

    fn intersections_arc_segment_tol(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections;

    fn intersections_segment_tol<'a, T>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections
    where
        T: Into<SegmentRef<'a>>,
    {
        let segment: SegmentRef = segment.into();
        match segment {
            SegmentRef::LineSegment(line_segment) => {
                self.intersections_line_segment_tol(line_segment, epsilon, max_relative)
            }
            SegmentRef::ArcSegment(arc_segment) => {
                self.intersections_arc_segment_tol(arc_segment, epsilon, max_relative)
            }
        }
    }
}

impl crate::private::Sealed for [f64; 2] {}
impl<'p> crate::private::Sealed for ToleranceContext<'p, [f64; 2]> {}

impl WithTolerance for [f64; 2] {}

impl Primitive for [f64; 2] {
    fn covers<'a, T>(&self, other: T) -> bool
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>,
    {
        self.with_tolerance(DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
            .covers(other)
    }

    fn intersections_primitive<'a, T: Primitive>(&self, other: &'a T) -> PrimitiveIntersections
    where
        &'a T: Into<GeometryRef<'a>>,
        Self: Sized,
    {
        self.with_tolerance(DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
            .intersections_primitive(other)
    }

    fn intersections<'a, T>(&self, other: T) -> Vec<crate::composite::Intersection>
    where
        Self: Sized,
        for<'b> &'b Self: Into<GeometryRef<'b>>,
        T: Into<GeometryRef<'a>>,
    {
        self.with_tolerance(DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
            .intersections(other)
    }
}

impl<'c> Primitive for ToleranceContext<'c, [f64; 2]> {
    fn covers<'a, T>(&self, other: T) -> bool
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>,
    {
        self.inner
            .covers_tol(other, self.epsilon, self.max_relative)
    }

    fn intersections_primitive<'a, T: Primitive>(&self, other: &'a T) -> PrimitiveIntersections
    where
        &'a T: Into<GeometryRef<'a>>,
        Self: Sized,
    {
        let geo_ref: GeometryRef = other.into();
        match geo_ref {
            GeometryRef::Point(pt) => {
                pt.intersections_point_tol(self.inner, self.epsilon, self.max_relative)
            }
            GeometryRef::ArcSegment(arc_segment) => {
                arc_segment.intersections_point_tol(self.inner, self.epsilon, self.max_relative)
            }
            GeometryRef::LineSegment(line_segment) => {
                line_segment.intersections_point_tol(self.inner, self.epsilon, self.max_relative)
            }
            GeometryRef::Line(line) => {
                line.intersections_point_tol(self.inner, self.epsilon, self.max_relative)
            }
            GeometryRef::Segment(segment) => {
                segment.intersections_point_tol(self.inner, self.epsilon, self.max_relative)
            }
            _ => unreachable!(
                "since other is a Primitive and Primitive is sealed, it cannot be another type"
            ),
        }
    }

    fn intersections<'a, T>(&self, other: T) -> Vec<crate::composite::Intersection>
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>,
    {
        let geo_ref: GeometryRef<'_> = self.inner.into();
        return geo_ref.intersections_tol::<true, _>(other, self.epsilon, self.max_relative);
    }
}

impl PrimitiveWithTol for [f64; 2] {
    fn covers_point_tol(&self, point: &[f64; 2], epsilon: f64, max_relative: f64) -> bool {
        return self.relative_eq(&point, epsilon, max_relative);
    }

    fn covers_arc_segment_tol(
        &self,
        _arc_segment: &ArcSegment,
        _epsilon: f64,
        _max_relative: f64,
    ) -> bool {
        return false;
    }

    fn covers_line_segment_tol(
        &self,
        _line_segment: &LineSegment,
        _epsilon: f64,
        _max_relative: f64,
    ) -> bool {
        return false;
    }

    fn covers_line_tol(&self, _line: &Line, _epsilon: f64, _max_relative: f64) -> bool {
        return false;
    }

    fn intersections_line_tol(
        &self,
        line: &Line,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        line.intersections_point_tol(self, epsilon, max_relative)
    }

    fn intersections_line_segment_tol(
        &self,
        line_segment: &LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        line_segment.intersections_point_tol(self, epsilon, max_relative)
    }

    fn intersections_arc_segment_tol(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        arc_segment.intersections_point_tol(self, epsilon, max_relative)
    }
}
