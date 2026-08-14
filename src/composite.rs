/*!
This module contains the definition of the [`Composite`] trait as well as the
[`Intersection`] and [`SegmentKey`] structs.

The geometric types in this crate are either "primitives" (which implement the
[`Primitive`] trait) or "composites" (which are defined from multiple primitives
and implement the [`Composite`] trait). See the trait documentation for details.
 */

use rayon::prelude::*;

use crate::ToleranceContext;
use crate::contour::Contour;
use crate::geometry::GeometryRef;
use crate::line::Line;
use crate::polysegment::Polysegment;
use crate::primitive::PrimitiveWithTol;
use crate::segment::{Segment, SegmentPolygonizer, SegmentRef};
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
    /// the key is not used on a [`Shape`].
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
    /// Sealed trait for [`Composite`]
    pub trait Sealed {}
}

// =============================================================================

/// All points of `other` are contained within `self`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum Contained {
    /// `other` is contained in `self`.
    Inside,
}

impl std::fmt::Display for Contained {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "inside self")
    }
}

/// Not all points of `other` are contained within `self`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum NotContained {
    /// A segment of `other` intersects a boundary of `self` at the given
    /// intersection.
    Intersection(Intersection),
    /// `other` is within the nth hole of `self`, where n is given as the
    /// anonymous field of the variant. Can only be created if `self` is a
    /// [`Shape`].
    InsideHole(usize),
    /// The [`BoundingBox`](bounding_box::BoundingBox) of `other` is not
    /// contained inside the [`BoundingBox`](bounding_box::BoundingBox) of
    /// `self`, therefore `other` cannot be contained in `self`.
    OutsideBoundingBox,
    /// `other` is not contained inside the contour of `self`. This is a generic
    /// fallback variant if no more specific reason is available.
    OutsideContour,
    /// `self` is a [`Polysegment`] which has no surface area and therefore
    /// cannot contain anything.
    NoSurfaceArea,
    /// A [`Segment`] of `other` lies on the boundary of `self`. The segment can
    /// be retrieved from `other` via the given [`SegmentKey`].
    OnBoundary(SegmentKey),
    /// The given point of `other` is outside `self`.
    PointOutside([f64; 2]),
}

impl std::fmt::Display for NotContained {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NotContained::Intersection(intersection) => {
                write!(f, "self and other intersect at {:?}", intersection)
            }
            NotContained::InsideHole(i) => write!(f, "other inside hole {i} of the self"),
            NotContained::OutsideBoundingBox => write!(f, "outside bounding box of self"),
            NotContained::OutsideContour => write!(f, "outside contour of self"),
            NotContained::OnBoundary(segment_key) => {
                write!(f, "other on boundary segment {:?} of self", segment_key)
            }
            NotContained::PointOutside(pt) => write!(f, "point {:?} of other outside self", pt),
            NotContained::NoSurfaceArea => write!(
                f,
                "self has no surface area and can therefore not contain other"
            ),
        }
    }
}

impl std::error::Error for NotContained {}

/// All points of `other` are covered by `self`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum Covered {
    /// `other` is covered by `self`.
    Inside,
    /// `other` lies on a boundary [`Segment`] of `self`, the segment can be
    /// retrieved by the given [`SegmentKey`].
    OnBoundary(SegmentKey),
}

impl std::fmt::Display for Covered {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Covered::Inside => write!(f, "inside"),
            Covered::OnBoundary(segment_key) => {
                write!(f, "other on boundary segment {:?} of self", segment_key)
            }
        }
    }
}

/// Not all points of `other` are covered by `self`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum NotCovered {
    /// `other` is within the nth hole of `self`, where n is given as the
    /// anonymous field of the variant. Can only be created if `self` is a
    /// [`Shape`].
    InsideHole(usize),
    /// A segment of `other` intersects a boundary of `self` at the given
    /// intersection and is at least partially outside `self`.
    Intersection(Intersection),
    /// `other` is not covered by `self`. This is a generic/ fallback variant if
    /// no more specific reason is available.
    Outside,
    /// The specified segment of `other` is not covered by `self`.
    SegmentNotCovered(SegmentKey),
    /// The [`BoundingBox`](bounding_box::BoundingBox) of `other` is not
    /// covered by the [`BoundingBox`](bounding_box::BoundingBox) of `self`,
    /// therefore `other` cannot be covered by `self`.
    OutsideBoundingBox,
}

impl std::fmt::Display for NotCovered {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NotCovered::InsideHole(i) => write!(f, "other inside hole {i} of self"),
            NotCovered::Intersection(intersection) => write!(
                f,
                "other partially not covered by self, intersect at {:?}",
                intersection
            ),
            NotCovered::Outside => write!(f, "outside of self"),
            NotCovered::OutsideBoundingBox => write!(f, "outside bounding box"),
            NotCovered::SegmentNotCovered(segment_key) => write!(
                f,
                "segment {:?} of other is not covered by self",
                segment_key
            ),
        }
    }
}

impl std::error::Error for NotCovered {}

/// `self` contains at least one point of `other`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum Overlap {
    /// The given point of `other` is contained in `self`.
    Point([f64; 2]),
    /// `self` and `other` are identical (and therefore overlap by definition).
    Identical,
    /// The [`Segment`] specified by `key` of `self` overlaps `other`, of
    /// `key_of_self` is true. Otherwhise, the `key` specifies a segment of
    /// `other` which overlaps with `self`.
    Segment {
        /// Key of the [`Segment`] which overlaps.
        key: SegmentKey,
        /// Whether the `key` points to a [`Segment`] of `self` or `other`.
        key_of_self: bool,
    },
    /// Depending on the value of `self_covers_other`, either `self` covers
    /// `other` or `other` covers `self`.
    Covers {
        /// Whether `self` covers `other` or `other` covers `self`.
        self_covers_other: bool,
    },
}

impl std::fmt::Display for Overlap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Overlap::Point(pt) => write!(f, "point {:?} inside contour", pt),
            Overlap::Identical => write!(f, "self and other are identical"),
            Overlap::Segment { key, key_of_self } => {
                if *key_of_self {
                    write!(f, "segment {:?} of self overlaps with other", key)
                } else {
                    write!(f, "segment {:?} of other overlaps with self", key)
                }
            }
            Overlap::Covers { self_covers_other } => {
                if *self_covers_other {
                    write!(f, "self covers other")
                } else {
                    write!(f, "other covers self")
                }
            }
        }
    }
}

impl Overlap {
    pub(crate) fn switch(self) -> Self {
        match self {
            Overlap::Segment { key, key_of_self } => Overlap::Segment {
                key,
                key_of_self: !key_of_self,
            },
            Overlap::Covers { self_covers_other } => Overlap::Covers {
                self_covers_other: !self_covers_other,
            },
            _ => return self,
        }
    }
}

/// `self` contains not a single point of `other`. See docstring of
/// [`Composite`] for more information.
#[derive(Clone, Debug, PartialEq)]
pub enum NoOverlap {
    /// `self` is a [`Polysegment`] which has no surface area and therefore
    /// cannot contain anything.
    NoSurfaceArea,
    /// The [`BoundingBox`](bounding_box::BoundingBox)es of `self` and `other`
    /// do not overlap, hence `self` and `other` cannot overlap either.
    DisjointBoundingBoxes,
    /// No point of `other` is contained in `self`. This is a generic fallback
    /// in case no more specific reason has been found.
    NoPointContained,
    /// `other` is within the nth hole of `self`, where n is given as the
    /// anonymous field of the variant. Can only be created if `self` is a
    /// [`Shape`].
    InsideHole(usize),
}

impl std::fmt::Display for NoOverlap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NoOverlap::NoSurfaceArea => write!(
                f,
                "self has no surface area and can therefore not contain any point of other"
            ),
            NoOverlap::DisjointBoundingBoxes => write!(f, "bounding boxes do not overlap"),
            NoOverlap::NoPointContained => write!(f, "no point of other is contained within self"),
            NoOverlap::InsideHole(i) => write!(f, "other inside hole {i} of self"),
        }
    }
}

impl std::error::Error for NoOverlap {}

// =============================================================================

/**
A trait for types composed of multiple [`Segment`]s:
[`Polysegment`]s, [`Contour`]s and [`Shape`]s.

This trait provides methods for retrieving properties shared between all
composite types (like e.g. the number of segments), coverage, containment and
overlap checks and intersection calculation. It is not meant to be implemented
for other types, hence it is
[sealed](<https://rust-lang.github.io/api-guidelines/future-proofing.htm>).

# Properties

The following methods are available:
- [`Composite::centroid`] for finding the center of mass / centroid,
- [`Composite::num_segments`] for reading out the number of segments and
- [`Composite::segment`] for retrieving a reference to a particular segment with
a [`SegmentKey`].

# Containment and coverage

This trait provides a variety of methods to check if `self` contains, covers or
partially contains another geometric object such as a point, a [`Segment`] or
another [`Composite`]. These methods start with `contains_`, `contains_any` or
`covers_` followed by the type name of the second object (e.g.
[`Composite::contains_contour`]). They return a `Result<OkEnum, ErrEnum>` with
the `Ok` variant signalling that the check succeded (e.g. in case of
[`Composite::contains_contour`] the contour is in fact contained within `self`)
and `Err` obviously indicating that the check failed. The enums contain the
first reason the algorithm found which unambiguously proves the result. In fact,
there might often be multiple reasons. For example, when a [`Contour`] is not
contained in `self`, this might be due to the bounding boxes being disjoint and
due to a particular point not being contained. In this case, the algorithm would
first check the bounding boxes and return that reason
[`NotContained::OutsideBoundingBox`] wrappend in `Err`.

## Containment

**Return types**: [`Contained`]` and `[`NotContained`]`

**Function names**: `contains_` + type name

A composite "contains" another geometric entity if if all points of the latter
are within it (and in contrast to the concept of "covers", not on its
boundaries). This means that a [`Polysegment`] cannot contain anything, since
it has no surface area and only consists of its boundary.

## Partial containment

**Return types**: [`Overlap`]` and `[`NoOverlap`]`

**Function names**: `contains_any` + type name

A composite "contains any" part of another geometric entity if at least one
point of the latter lies within the composite's interior (excluding
boundaries).

## Coverage

**Return types**: [`Covered`]` and `[`NotCovered`]`

**Function names**: `covers_` + type name

Similar to the [`Primitive`] trait, a composite "covers" another geometric
entity if all points of the latter are within it or on its boundaries.

# Intersection

Different to primitives, there can be an arbitrary number of intersections
between two intersections, which is why the corresponding methods return
an iterator. Since intersection calculation between composites can be easily
parallelized, there is a serial and a parallel variant (the latter returning
a [`ParallelIterator`] and having a `_par` suffix).`

These iterators always return an [`Intersection`] object. The "left" side of the
object refers to the first argument `self`, whereas the "right" side
refers to the type of the other composite. See the docstring of [`Intersection`]
for more.

In contrast to primitives, composites can self-intersect. The self-intersection
points can be calculated by using `self` as the second argument `other`:

```
use planar_geo::prelude::*;

let sc = Polysegment::from_points(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]);

let mut iter = sc.intersections_polysegment(&sc);

// Intersection in the "cross" middle
assert_eq!(iter.next().unwrap().point, [0.5, 0.5]);

// The polysegment start and end point "touch", this is also counted as an intersection
assert_eq!(iter.next().unwrap().point, [0.0, 0.0]);
assert!(iter.next().is_none());

// By contrast, a contour is closed by default, hence start and end point are
// connected and therefore are not an intersection
let c = Contour::new(sc.clone());
let mut iter = c.intersections_contour(&c);

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
takes `epsilon` and `max_relative` as additional arguments to specify a certain
absolute and relative tolerance for intersection detection.

All intersection functions first check if the bounding boxes of the two
primitives overlap (short-circuiting the evaluation if they don't). Hence, it is
not necessary to check this before calling an intersection method.
 */
pub trait Composite: private::Sealed + Sync {
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
        assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9]).unwrap())));
    }

    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    // SegmentKey of Shape is SegmentKey. This accessor retrieves the third
    // segment of the second contour / first hole of the shape
    let segment = shape.segment(SegmentKey::new(1, 2));
    assert_eq!(segment, Some(&Segment::from(LineSegment::new([0.9, 0.9], [0.1, 0.9]).unwrap())));
    ```
     */
    fn segment(&self, key: SegmentKey) -> Option<&Segment>;

    /**
    Returns an iterator over all segments of `self` and their keys.

    If `self` is a [`Shape`], the segments of the outer contour are returned
    first, followed by those of the holes (in arbitrary order)
     */
    fn iter<'a>(&'a self) -> impl Iterator<Item = (SegmentKey, &'a Segment)>;

    /**
    Returns a parallel iterator over all segments of `self` and their keys in
    arbitrary order.
     */
    fn par_iter<'a>(&'a self) -> impl ParallelIterator<Item = (SegmentKey, &'a Segment)>;

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

    fn contains<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Contained, NotContained>;

    fn contains_any<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Overlap, NoOverlap>;

    fn covers<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Covered, NotCovered>;

    fn intersections_point<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl Iterator<Item = Intersection> + 'a;

    fn intersections_point_par<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    fn intersections_line<'a>(&'a self, line: &'a Line) -> impl Iterator<Item = Intersection> + 'a;

    fn intersections_line_par<'a>(
        &'a self,
        line: &'a Line,
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

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

    let ls = LineSegment::new([-1.0, 0.5], [2.0, 0.5]).expect("valid inputs");
    let mut intersections = polysegment.intersections_segment(&ls);

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
    fn intersections_segment<'a, T>(
        &'a self,
        segment: T,
    ) -> impl Iterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>;

    /**
    Returns a parallelized iterator over all intersections of `self` with the
    `primitive`.

    This is the parallelized variant of [`Composite::intersections_primitive`].
    See its docstring for more information and examples.
    TODO: [`Composite::intersections`] is more generic, but allocates eagerly!
     */
    fn intersections_segment_par<'a, T>(
        &'a self,
        segment: T,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>;

    /**
    Returns a iterator over all intersections of `self` with the `polysegment`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.
    TODO: [`Composite::intersections`] is more generic, but allocates eagerly!

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let left = Polysegment::from_points(vertices);

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let right = Polysegment::from_points(vertices);

    let mut intersections = left.intersections_polysegment(&right);

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns a iterator over all intersections of `self` with the `contour`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.
    TODO: [`Composite::intersections`] is more generic, but allocates eagerly!

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));

    let mut intersections = polysegment.intersections_contour(&contour);

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    /**
    Returns an iterator over all intersections of `self` with the `shape`.

    This method is mainly used to implement
    [`Composite::intersections_composite`], consider using that method instead.
    TODO

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

    let mut intersections = polysegment.intersections_shape(&shape);

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
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

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

    let intersections = polysegment.intersections(&shape);
    assert_eq!(intersections.len(), 4);
    ```
     */
    fn intersections<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized;

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

    let intersections = polysegment.intersections_par(&shape);
    assert_eq!(intersections.len(), 4);
    ```
     */
    fn intersections_par<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized;
}

/// This trait provides the actual implementation for [`Composite`] methods
/// where a tolerance is required for the underlying algorithm. For example, the
/// implementation of [`Composite::contains_point`] looks like this for
/// [`Polysegment`], [`Contour`] or [`Shape`]:
///
/// ```ignore
/// fn contains_point(&self, point: [f64; 2]) -> Result<Contained, NotContained> {
///     self.contains_point_tol(&point)
/// }
/// ```
/// If custom tolerances are needed, a [`ToleranceContext`] can be created
/// with `with_tolerance` and the specified tolerances then replace
/// [`DEFAULT_EPSILON`](crate::DEFAULT_EPSILON) and
/// [`DEFAULT_MAX_RELATIVE`](crate::DEFAULT_MAX_RELATIVE) in that operation:
///
/// `contour.with_tolerance(1e-9, 1e-8).contains_point(&point)
///
/// results in
///
/// `contour.contains_point_tol(&point, 1e-9, 1e-8)
pub(crate) trait CompositeWithTol: Composite {
    fn contains_point_tol(
        &self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained>;

    fn contains_segment_tol<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained>;

    fn contains_polysegment_tol(
        &self,
        polysegment: &Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        if let Some(err) = polysegment.segments_par().find_map_any(|s| {
            match self.contains_segment_tol(s, epsilon, max_relative) {
                Ok(_) => None,
                Err(err) => Some(err),
            }
        }) {
            return Err(err);
        }
        return Ok(Contained::Inside);
    }

    fn contains_contour_tol(
        &self,
        contour: &Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        return self.contains_polysegment_tol(contour.polysegment(), epsilon, max_relative);
    }

    fn contains_shape_tol(
        &self,
        shape: &Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        return self.contains_contour_tol(shape.contour(), epsilon, max_relative);
    }

    fn contains_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained>;

    fn contains_any_segment_tol<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap>;

    fn contains_any_polysegment_tol(
        &self,
        polysegment: &Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap> {
        let hit = polysegment
            .segments_par()
            .enumerate()
            .find_map_any(|(segment_idx, s)| {
                // Returns a segment of `polysegment` that overlaps with `self`.
                // Therefore, `segment_of_self` is false.
                match self.contains_any_segment_tol(s, epsilon, max_relative) {
                    Ok(_) => Some(Overlap::Segment {
                        key: SegmentKey::from_segment_idx(segment_idx),
                        key_of_self: false,
                    }),
                    Err(_) => None,
                }
            });

        return hit.ok_or(NoOverlap::NoPointContained);
    }

    fn contains_any_contour_tol(
        &self,
        contour: &Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap>;

    fn contains_any_shape_tol(
        &self,
        shape: &Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap>;

    fn contains_any_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap>;

    fn covers_point_tol(
        &self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered>;

    fn covers_segment_tol<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered>;

    fn covers_polysegment_tol(
        &self,
        polysegment: &Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        let first_seg = polysegment.front().ok_or(NotCovered::Outside)?;
        let first_cover = self.covers_segment_tol(first_seg, epsilon, max_relative)?;

        // Skip the first element, because we already tested it. If all other
        // segments are covered as well, return the result of the first
        // segment.
        match polysegment
            .segments_par()
            .enumerate()
            .skip(1)
            .find_map_any(
                |(i, s)| match self.covers_segment_tol(s, epsilon, max_relative) {
                    Ok(_) => None,
                    Err(_) => Some(i),
                },
            ) {
            Some(i) => {
                return Err(NotCovered::SegmentNotCovered(SegmentKey::from_segment_idx(
                    i,
                )));
            }
            None => return Ok(first_cover),
        }
    }

    fn covers_contour_tol(
        &self,
        contour: &Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        return self.covers_polysegment_tol(contour.polysegment(), epsilon, max_relative);
    }

    fn covers_shape_tol(
        &self,
        shape: &Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        return self.covers_contour_tol(shape.contour(), epsilon, max_relative);
    }

    fn covers_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered>;

    fn intersections_point_tol<'a>(
        &'a self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a {
        let point = point.clone();
        self.iter()
            .map(move |(k, s)| {
                s.intersections_point_tol(&point, epsilon, max_relative)
                    .into_iter()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_point_par_tol<'a>(
        &'a self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        let point = point.clone();
        self.par_iter()
            .map(move |(k, s)| {
                s.intersections_point_tol(&point, epsilon, max_relative)
                    .into_iter()
                    .par_bridge()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_line_tol<'a>(
        &'a self,
        line: &'a Line,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.iter()
            .map(move |(k, s)| {
                s.intersections_line_tol(line, epsilon, max_relative)
                    .into_iter()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_line_par_tol<'a>(
        &'a self,
        line: &'a Line,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.par_iter()
            .map(move |(k, s)| {
                s.intersections_line_tol(line, epsilon, max_relative)
                    .into_iter()
                    .par_bridge()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_segment_tol<'a, T>(
        &'a self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        let segment_ref: SegmentRef = segment.into();
        self.iter()
            .map(move |(k, s)| {
                s.intersections_segment_tol(segment_ref, epsilon, max_relative)
                    .into_iter()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_segment_par_tol<'a, T>(
        &'a self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        let segment_ref: SegmentRef = segment.into();
        self.par_iter()
            .map(move |(k, s)| {
                s.intersections_segment_tol(segment_ref, epsilon, max_relative)
                    .into_iter()
                    .par_bridge()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_polysegment_tol<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a;

    fn intersections_polysegment_par_tol<'a>(
        &'a self,
        polysegment: &'a Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    fn intersections_contour_tol<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_polysegment_tol(contour.polysegment(), epsilon, max_relative)
    }

    fn intersections_contour_par_tol<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_polysegment_par_tol(contour.polysegment(), epsilon, max_relative)
    }

    fn intersections_shape_tol<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a;

    fn intersections_shape_par_tol<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a;

    fn intersections_primitive_par_tol<'a, T>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        &'a T: Into<GeometryRef<'a>>,
        T: PrimitiveWithTol,
    {
        self.par_iter()
            .map(move |(k, s)| {
                s.intersections_primitive_tol(primitive, epsilon, max_relative)
                    .into_iter()
                    .par_bridge()
                    .map(move |point| Intersection {
                        point,
                        left: k,
                        right: Default::default(),
                    })
            })
            .flatten()
    }

    fn intersections_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Vec<Intersection>
    where
        Self: Sized;

    fn intersections_par_tol<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Vec<Intersection>
    where
        Self: Sized;
}

impl<'c, T: CompositeWithTol> private::Sealed for ToleranceContext<'c, T> {}

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
        line: SegmentPolygonizer,
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
    pub fn segment_polygonizer<'a, T>(&self, segment: T, index: usize) -> SegmentPolygonizer
    where
        T: Into<SegmentRef<'a>>,
    {
        let segment: SegmentRef = segment.into();
        match self {
            Polygonizer::PerType { arc, line } => match segment {
                SegmentRef::LineSegment(_) => return *line,
                SegmentRef::ArcSegment(_) => return *arc,
            },
            Polygonizer::Individual { default, map } => return *map.get(&index).unwrap_or(default),
        };
    }
}

impl Default for Polygonizer {
    fn default() -> Self {
        return Polygonizer::PerType {
            arc: SegmentPolygonizer::default(),
            line: SegmentPolygonizer::default(),
        };
    }
}

/**
An iterator over all points created with [`Polysegment::polygonize`] /
[`Contour::polygonize`].
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
