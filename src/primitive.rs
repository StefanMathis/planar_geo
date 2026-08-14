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
    DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE, ToleranceContext, WithTolerance,
    composite::CompositeWithTol,
    geometry::GeometryRef,
    line::Line,
    segment::{ArcSegment, LineSegment, Segment, SegmentRef},
};
use approx::RelativeEq;
use rayon::prelude::*;

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

A primitive "covers" another one if all points of the latter are part of the
former. The generic [`Primitive::covers`] method provides a unified interface to
check this. Specialized `covers_` methods for specific types are available as
well and are particularily useful when using [`Primitive`] to define a trait
object (which cannot use generic functions).

By definition, a primitive always covers itself.

```
use planar_geo::prelude::*;

let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();;
assert!(ls.covers(&[0.5, 0.0]));
assert!(!ls.covers(&[0.5, 0.1]));
assert!(ls.covers(&ls));
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

let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();;

// Self-intersection
assert_eq!(
    ls.intersections_primitive(&ls),
    PrimitiveIntersections::Zero
);

// Intersections with equal primitive
let ls_cloned = ls.clone();
assert_eq!(
    ls.intersections_primitive(&ls_cloned),
    PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]])
);
```
 */
pub trait Primitive: crate::private::Sealed + Sync {
    fn covers<'a, T>(&self, other: T) -> bool
    where
        Self: Sized,
        T: Into<GeometryRef<'a>>;

    fn intersections_primitive<'a, T: Primitive>(&self, other: &'a T) -> PrimitiveIntersections
    where
        &'a T: Into<GeometryRef<'a>>,
        Self: Sized;

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

    fn intersections_primitive_tol<T>(
        &self,
        primitive: &T,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections
    where
        T: PrimitiveWithTol;
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
        let geo_ref: GeometryRef<'_> = other.into();
        match geo_ref {
            GeometryRef::Point(pt) => pt
                .intersections_point_tol(self.inner, self.epsilon, self.max_relative)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::ArcSegment(arc_segment) => arc_segment
                .intersections_point_tol(self.inner, self.epsilon, self.max_relative)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::LineSegment(line_segment) => line_segment
                .intersections_point_tol(self.inner, self.epsilon, self.max_relative)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::Line(line) => line
                .intersections_point_tol(self.inner, self.epsilon, self.max_relative)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::Segment(segment) => segment
                .intersections_point_tol(self.inner, self.epsilon, self.max_relative)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::BoundingBox(bounding_box) => {
                let c = crate::contour::Contour::from(bounding_box);
                c.intersections_point_par_tol(self.inner, self.epsilon, self.max_relative)
                    .collect()
            }
            GeometryRef::Polysegment(polysegment) => polysegment
                .intersections_point_par_tol(self.inner, self.epsilon, self.max_relative)
                .collect(),
            GeometryRef::Contour(contour) => contour
                .intersections_point_par_tol(self.inner, self.epsilon, self.max_relative)
                .collect(),
            GeometryRef::Shape(shape) => shape
                .intersections_point_par_tol(self.inner, self.epsilon, self.max_relative)
                .collect(),
        }
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

    fn intersections_primitive_tol<T>(
        &self,
        primitive: &T,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections
    where
        T: PrimitiveWithTol,
    {
        primitive.intersections_point_tol(self, epsilon, max_relative)
    }
}
