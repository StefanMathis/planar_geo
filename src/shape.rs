/*!
Defines the [`Shape`] type, a [`Composite`] representing the external boundaries
of an object in 2D space.

A [`Shape`] is made up of multiple [`Contour`]s, where one contour describes the
object extents ("outer contour"), whereas all other contours define "holes"
within the object.

Most users should interact with this module through the [`Shape`] type
itself; see its documentation for details on construction, invariants, and
usage.
 */

use crate::CentroidData;
use crate::ToleranceContext;
use crate::composite::*;
use crate::error::ShapeConstructorError;
use crate::geometry::GeometryRef;
use crate::polysegment::Polysegment;
use crate::segment::{Segment, SegmentRef};
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE};
use crate::{Transformation, contour::Contour};
use bounding_box::{BoundingBox, ToBoundingBox};
use rayon::prelude::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/**
A collection of [`Contour`]s which describe the boundaries / extents of a
geometric body in 2D space.

*/
#[doc = ""]
#[cfg_attr(feature = "doc-images", doc = "![Shape example][example_shape]")]
#[cfg_attr(
    feature = "doc-images",
    embed_doc_image::embed_doc_image("example_shape", "docs/img/example_shape.svg")
)]
#[cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images not enabled**. Compile docs with
    `cargo doc --features 'doc-images'` and Rust version >= 1.54."
)]
/**

A [`Shape`] always has at least one [`Contour`] which describes the outer
extents of the body represented by it. Any additional contours are interpreted
as "holes" within the "outer" contour. These contours need to fulfill the
following invariants:
- The outer contour contains all holes (see [`Contour::contains`]).
- No hole may contain another hole
- No contour intersects with itself or any other contour.
- No contour might be empty (i.e. contain no
[`Segment`](Segment)).

# Constructing and modifying a shape

The standard constructor which checks all the aforementioned invariants is
[`new`](Shape::new), which takes a [`Vec<Contour`](Contour), interprets the first
element as the outer contour and all other elements as holes. Alternatively, a
shape can be constructed from its outer contour via
[`from_outer`](Shape::from_outer) and holes can be added incrementally with
[`add_hole`](Shape::add_hole). Holes can also be removed with
[`remove_hole`](Shape::remove_hole).

# Access of individual contours and segments

Individual segments can be retrieved with [`Composite::segment`]. The outer
contour is available with [`contour`](Shape::contour) method, while a slice
of all holes can be accessed with [`holes`](Shape::holes). The
[`contours`](Shape::contours) method returns a slice where the first element is
the outer contour and all other elements are the holes. A shape can also be
converted into the underlying [`Vec<Contour>`](Contour) with its [`From`]
implementation.

# Serialization and deserialization

When the `serde` feature is enabled, a shape can be serialized and
deserialized using the [serde] crate. It uses the same serialized
representation as a [`Vec<Contour>`].
 */
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Shape(Vec<Contour>);

impl Shape {
    /**
    Creates a new [`Shape`] out of the given `contours`. The first element is
    interpreted as the outer contour of the shape, the other elements are
    interpreted as holes.

    The given vector of [`Contour`]s must fulfill the following conditions:
    - It must not be empty.
    - None of its elements (contours) must be empty.
    - None of the contours must intersect itself or each other.
    - The first ("outer") contour must contain all other contours / "holes".
    - No hole is contained within another hole.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Shape fulfills all conditions
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_ok());

    // Given vector is empty
    assert!(Shape::new(Vec::new()).is_err());

    // One of the contours is empty
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));
    let c2 = Contour::new(Polysegment::new());
    assert!(Shape::new(vec![c1, c2]).is_err());

    // The contours intersect
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [1.1, 0.1], [1.1, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_err());

    // Second contour is not inside the first contour
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[1.0, 0.1], [2.0, 0.1], [2.0, 0.9], [1.0, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_err());
    ```
     */
    pub fn new(contours: Vec<Contour>) -> Result<Self, ShapeConstructorError<Vec<Contour>>> {
        if contours.len() == 0 {
            return Err(ShapeConstructorError::EmptyVec);
        }
        for (idx, contour) in contours.iter().enumerate() {
            if contour.is_empty() {
                return Err(ShapeConstructorError::EmptyContour {
                    input: contours,
                    idx,
                });
            }
        }
        let this = Shape(contours);

        let outer = this.contour();
        for (first_hole_idx, first_hole) in this.holes().iter().enumerate() {
            if let Err(e) = outer.contains_contour(&first_hole) {
                return Err(ShapeConstructorError::HoleOutsideContour {
                    input: this.0,
                    reason: e,
                    idx: first_hole_idx + 1, // First element of this.0 is the outer contour.
                });
            }

            for (second_hole_idx, second_hole) in
                this.holes()[(first_hole_idx + 1)..].iter().enumerate()
            {
                if let Ok(reason) = first_hole.contains_contour(&second_hole) {
                    return Err(ShapeConstructorError::HoleInsideHole {
                        input: this.0,
                        reason,
                        outer_hole_idx: first_hole_idx,
                        inner_hole_idx: second_hole_idx,
                    });
                }

                if let Ok(reason) = second_hole.contains_contour(&first_hole) {
                    return Err(ShapeConstructorError::HoleInsideHole {
                        input: this.0,
                        reason,
                        outer_hole_idx: second_hole_idx,
                        inner_hole_idx: first_hole_idx,
                    });
                }
            }
        }

        // Check if the shape intersects itself
        if let Some(intersection) = this
            .intersections_shape_par(&this)
            .find_map_any(|v| Some(v))
        {
            return Err(ShapeConstructorError::Intersection {
                input: this.0,
                intersection,
            });
        }

        return Ok(this);
    }

    /**
    Creates a contour with no holes from its outer contour.

    The contour must not be empty and must not intersect itself.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Contour fulfill al criteria
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c = Contour::new(Polysegment::from_points(vertices));
    assert!(Shape::from_outer(c).is_ok());

    // Contour is empty
    let c = Contour::new(Polysegment::new());
    assert!(Shape::from_outer(c).is_err());

    // Contour intersects itself
    let vertices = &[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]];
    let c = Contour::new(Polysegment::from_points(vertices));
    assert!(Shape::from_outer(c).is_err());
    ```
     */
    pub fn from_outer(outer: Contour) -> Result<Self, ShapeConstructorError<Vec<Contour>>> {
        if outer.is_empty() {
            return Err(ShapeConstructorError::EmptyContour {
                input: vec![outer],
                idx: 0,
            });
        }
        if let Some(i) = outer
            .intersections_contour_par(&outer)
            .find_map_any(|v| Some(v))
        {
            return Err(ShapeConstructorError::Intersection {
                input: vec![outer],
                intersection: Intersection {
                    point: i.point,
                    left: i.left,
                    right: i.right,
                },
            });
        }

        return Ok(Shape(vec![outer]));
    }

    /**
    Returns a reference to all [`Contour`]s defining `self`. The first contour
    is the "outer" contour (see [`Shape::contour`]), all other contours are
    interpreted as "holes" (see [`Shape::holes`]).
     */
    pub fn contours(&self) -> &[Contour] {
        return &self.0;
    }

    /**
    Returns the "outer" [`Contour`] of `self`. This is the first element of the
    underlying [`Vec<Contour>`].
     */
    pub fn contour<'a>(&'a self) -> &'a Contour {
        /*
        Safety: This function is safe since it is guaranteed during the Shape
        element construction that it contains at least the contour (e.g. the
        field "lines" has a length of at least 1).
        */
        return unsafe { self.0.get_unchecked(0) };
    }

    /**
    Returns all "holes" of `self`. If the [`Shape`] has no holes, the returned
    slice is empty.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let c3 = Contour::new(Polysegment::from_points(vertices));

    let shape = Shape::new(vec![c1, c2, c3]).expect("valid input");
    assert_eq!(shape.holes().len(), 2);
    assert_eq!(shape.holes().len() + 1, shape.contours().len());
    ```
     */
    pub fn holes(&self) -> &[Contour] {
        return &self.0[1..self.0.len()]; // Skip first element since it is the contour of the shape
    }

    /**
    Returns the surface area of `self`.

    This function calculates the surface area of all contours (via
    [`Contour::area`]) and then subtracts the areas of all holes from that of
    the outer contour.

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));

    let shape = Shape::new(vec![contour.clone(), hole.clone()]).expect("valid input");
    assert_eq!(shape.area(), contour.area() - hole.area());
    ```
     */
    pub fn area(&self) -> f64 {
        let contour_area = self.contour().area();
        let holes: f64 = self
            .holes()
            .par_iter()
            .map(|polysegment| return polysegment.area())
            .sum();
        return contour_area - holes;
    }

    /**
    Tries to add a hole to `self`. The given contour must fulfill the following
    conditions:
    - It must be within the outer contour of `self`.
    - It must not be empty.
    - It must not intersect any of the contours making up the shape.
    - It must not be within another hole.

    If any of these conditions are not met, an error is returned which contains
    the input `hole`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Shaoe without an hole
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut shape = Shape::new(vec![Contour::new(Polysegment::from_points(vertices))]).expect("valid input");
    assert_eq!(shape.holes().len(), 0);

    // Add a hole to the shape - this works
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_ok());
    assert_eq!(shape.holes().len(), 1);

    // Adding this hole does not work because it intersects with the outer contour
    let vertices = &[[1.0, 0.1], [1.0, 0.1], [2.0, 0.9], [2.0, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_err());
    assert_eq!(shape.holes().len(), 1);
    ```
     */
    pub fn add_hole(&mut self, hole: Contour) -> Result<(), ShapeConstructorError<Contour>> {
        if hole.is_empty() {
            return Err(ShapeConstructorError::EmptyContour {
                input: hole,
                idx: 0,
            });
        }

        // Is the hole self-intersecting?
        if let Some(i) = hole
            .intersections_contour_par(&hole)
            .find_map_any(|v| Some(v))
        {
            return Err(ShapeConstructorError::Intersection {
                input: hole,
                intersection: Intersection {
                    point: i.point,
                    left: i.left,
                    right: i.right,
                },
            });
        }

        let outer = self.contour();

        if let Err(e) = outer.contains_contour(&hole) {
            return Err(ShapeConstructorError::HoleOutsideContour {
                input: hole,
                reason: e,
                idx: 0, // Dummy value
            });
        }

        for (shape_hole_idx, shape_hole) in self.holes().iter().enumerate() {
            if let Ok(reason) = hole.contains_contour(shape_hole) {
                return Err(ShapeConstructorError::HoleInsideHole {
                    input: hole,
                    reason,
                    outer_hole_idx: shape_hole_idx,
                    inner_hole_idx: 0,
                });
            }

            if let Ok(reason) = shape_hole.contains_contour(&hole) {
                return Err(ShapeConstructorError::HoleInsideHole {
                    input: hole,
                    reason,
                    outer_hole_idx: 0,
                    inner_hole_idx: shape_hole_idx,
                });
            }
        }

        if let Some(intersection) = self
            .intersections_contour_par(&hole)
            .find_map_any(|v| Some(v))
        {
            return Err(ShapeConstructorError::Intersection {
                input: hole,
                intersection,
            });
        }

        self.0.push(hole);

        return Ok(());
    }

    /**
    Removes the hole with the given `index` from `self` and returns the
    associated contour. If no hole exists for the given `index`, `None` is
    returned instead.

    This operation performs a [`Vec::swap_remove`] on the underlying vector.
    This means that the order of holes can change (but the first element of
    the vector which represents the contour will not change position).

    # Examples

    ```
    use planar_geo::prelude::*;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let c3 = Contour::new(Polysegment::from_points(vertices));

    let mut shape = Shape::new(vec![c1, c2, c3]).unwrap();

    assert_eq!(shape.holes().len(), 2);
    assert!(shape.remove_hole(2).is_none());

    assert!(shape.remove_hole(1).is_some());
    assert_eq!(shape.holes().len(), 1);

    assert!(shape.remove_hole(0).is_some());
    assert_eq!(shape.holes().len(), 0);

    assert!(shape.remove_hole(0).is_none());
    assert_eq!(shape.holes().len(), 0);
    ```
     */
    pub fn remove_hole(&mut self, index: usize) -> Option<Contour> {
        if self.holes().len() > index {
            return Some(self.0.swap_remove(index + 1));
        }
        return None;
    }

    /**
    Returns the convex hull of the polygonized outer [`contour`](Shape::contour)
    of `self`.

    This method is implemented as `self.contour().convex_hull(polygonizer)`, see
    [`Contour::convex_hull`] for details and examples.
     */
    #[cfg(feature = "convex_hull")]
    pub fn convex_hull(
        &self,
        polygonizer: crate::composite::Polygonizer,
    ) -> planar_convex_hull::ConvexHullIter {
        return self.contour().convex_hull(polygonizer);
    }

    /// Wraps `self` in a [`ToleranceContext`] using the specified `epsilon` and
    /// `max_relative` tolerances.
    ///
    /// The [`ToleranceContext`] applies these tolerances to floating-point
    /// comparisons performed by geometric operations, such as finding
    /// intersections. See [`ToleranceContext`] for details and examples.
    pub fn with_tolerance<'a>(
        &'a self,
        epsilon: f64,
        max_relative: f64,
    ) -> ToleranceContext<'a, Self> {
        ToleranceContext {
            inner: self,
            epsilon,
            max_relative,
        }
    }
}

impl crate::composite::private::Sealed for Shape {}

impl Composite for Shape {
    fn segment(&self, key: SegmentKey) -> Option<&crate::segment::Segment> {
        let contour = self.0.get(key.contour_idx)?;
        return contour.segment(key);
    }

    fn num_segments(&self) -> usize {
        return self.0.iter().map(Composite::num_segments).sum();
    }

    fn centroid(&self) -> [f64; 2] {
        return CentroidData::from(self).into();
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self.0.iter().enumerate().flat_map(|(i, c)| {
            c.iter().map(move |(mut k, s)| {
                k.contour_idx = i;
                return (k, s);
            })
        });
    }

    fn par_iter<'a>(
        &'a self,
    ) -> impl ParallelIterator<Item = (SegmentKey, &'a crate::segment::Segment)> {
        return self.0.par_iter().enumerate().flat_map(|(i, c)| {
            c.par_iter().map(move |(mut k, s)| {
                k.contour_idx = i;
                return (k, s);
            })
        });
    }

    fn contains_point(&self, point: &[f64; 2]) -> Result<Contained, NotContained> {
        self.contains_point_tol(point, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Contained, NotContained> {
        self.contains_segment_tol(segment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_polysegment(&self, polysegment: &Polysegment) -> Result<Contained, NotContained> {
        self.contains_polysegment_tol(polysegment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_contour(&self, contour: &Contour) -> Result<Contained, NotContained> {
        self.contains_contour_tol(contour, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_shape(&self, shape: &Shape) -> Result<Contained, NotContained> {
        self.contains_shape_tol(shape, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Contained, NotContained> {
        self.contains_tol(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_any_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Overlap, NoOverlap> {
        self.contains_any_segment_tol(segment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_any_polysegment(&self, polysegment: &Polysegment) -> Result<Overlap, NoOverlap> {
        self.contains_any_polysegment_tol(polysegment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_any_contour(&self, contour: &Contour) -> Result<Overlap, NoOverlap> {
        self.contains_any_contour_tol(contour, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_any_shape(&self, shape: &crate::prelude::Shape) -> Result<Overlap, NoOverlap> {
        self.contains_any_shape_tol(shape, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn contains_any<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Overlap, NoOverlap> {
        self.contains_any_tol(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers_point(&self, point: &[f64; 2]) -> Result<Covered, NotCovered> {
        self.covers_point_tol(point, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Covered, NotCovered> {
        self.covers_segment_tol(segment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers_polysegment(&self, polysegment: &Polysegment) -> Result<Covered, NotCovered> {
        self.covers_polysegment_tol(polysegment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers_contour(&self, contour: &Contour) -> Result<Covered, NotCovered> {
        self.covers_contour_tol(contour, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers_shape(&self, shape: &Shape) -> Result<Covered, NotCovered> {
        self.covers_shape_tol(shape, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn covers<'a, T: Into<GeometryRef<'a>>>(&'a self, other: T) -> Result<Covered, NotCovered> {
        self.covers_tol(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_polysegment<'a>(
        &'a self,
        polysegment: &'a Polysegment,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_polysegment_tol(polysegment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_polysegment_par<'a>(
        &'a self,
        polysegment: &'a Polysegment,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_polysegment_par_tol(polysegment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_contour_tol(contour, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a Contour,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_contour_par_tol(contour, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_shape_tol(shape, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_shape_par_tol(shape, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_point<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_point_tol(point, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_point_par<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_point_par_tol(point, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_line<'a>(
        &'a self,
        line: &'a crate::prelude::Line,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.intersections_line_tol(line, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_line_par<'a>(
        &'a self,
        line: &'a crate::prelude::Line,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.intersections_line_par_tol(line, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_segment<'a, T>(&'a self, segment: T) -> impl Iterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        self.intersections_segment_tol(segment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_segment_par<'a, T>(
        &'a self,
        segment: T,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        self.intersections_segment_par_tol(segment, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        self.intersections_tol(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }

    fn intersections_par<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        self.intersections_par_tol(other, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
    }
}

impl<'c> Composite for ToleranceContext<'c, Shape> {
    fn segment(&self, key: SegmentKey) -> Option<&Segment> {
        self.inner.segment(key)
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = (SegmentKey, &'a Segment)> {
        self.inner.iter()
    }

    fn par_iter<'a>(&'a self) -> impl ParallelIterator<Item = (SegmentKey, &'a Segment)> {
        self.inner.par_iter()
    }

    fn num_segments(&self) -> usize {
        self.inner.num_segments()
    }

    fn centroid(&self) -> [f64; 2] {
        self.inner.centroid()
    }

    fn contains_point(&self, point: &[f64; 2]) -> Result<Contained, NotContained> {
        self.inner
            .contains_point_tol(point, self.epsilon, self.max_relative)
    }

    fn contains_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Contained, NotContained> {
        self.inner
            .contains_segment_tol(segment, self.epsilon, self.max_relative)
    }

    fn contains_polysegment(&self, polysegment: &Polysegment) -> Result<Contained, NotContained> {
        self.inner
            .contains_polysegment_tol(polysegment, self.epsilon, self.max_relative)
    }

    fn contains_contour(&self, contour: &Contour) -> Result<Contained, NotContained> {
        self.inner
            .contains_contour_tol(contour, self.epsilon, self.max_relative)
    }

    fn contains_shape(&self, shape: &Shape) -> Result<Contained, NotContained> {
        self.inner
            .contains_shape_tol(shape, self.epsilon, self.max_relative)
    }

    fn contains<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Contained, NotContained> {
        self.inner
            .contains_tol(other, self.epsilon, self.max_relative)
    }

    fn contains_any_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Overlap, NoOverlap> {
        self.inner
            .contains_any_segment_tol(segment, self.epsilon, self.max_relative)
    }

    fn contains_any_polysegment(&self, polysegment: &Polysegment) -> Result<Overlap, NoOverlap> {
        self.inner
            .contains_any_polysegment_tol(polysegment, self.epsilon, self.max_relative)
    }

    fn contains_any_contour(&self, contour: &Contour) -> Result<Overlap, NoOverlap> {
        self.inner
            .contains_any_contour_tol(contour, self.epsilon, self.max_relative)
    }

    fn contains_any_shape(&self, shape: &crate::prelude::Shape) -> Result<Overlap, NoOverlap> {
        self.inner
            .contains_any_shape_tol(shape, self.epsilon, self.max_relative)
    }

    fn contains_any<'a, T: Into<GeometryRef<'a>>>(&self, other: T) -> Result<Overlap, NoOverlap> {
        self.inner
            .contains_any_tol(other, self.epsilon, self.max_relative)
    }

    fn covers_point(&self, point: &[f64; 2]) -> Result<Covered, NotCovered> {
        self.inner
            .covers_point_tol(point, self.epsilon, self.max_relative)
    }

    fn covers_segment<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
    ) -> Result<Covered, NotCovered> {
        self.inner
            .covers_segment_tol(segment, self.epsilon, self.max_relative)
    }

    fn covers_polysegment(&self, polysegment: &Polysegment) -> Result<Covered, NotCovered> {
        self.inner
            .covers_polysegment_tol(polysegment, self.epsilon, self.max_relative)
    }

    fn covers_contour(&self, contour: &Contour) -> Result<Covered, NotCovered> {
        self.inner
            .covers_contour_tol(contour, self.epsilon, self.max_relative)
    }

    fn covers_shape(&self, shape: &Shape) -> Result<Covered, NotCovered> {
        self.inner
            .covers_shape_tol(shape, self.epsilon, self.max_relative)
    }

    fn covers<'a, T: Into<GeometryRef<'a>>>(&'a self, other: T) -> Result<Covered, NotCovered> {
        self.inner
            .covers_tol(other, self.epsilon, self.max_relative)
    }

    fn intersections_polysegment<'a>(
        &'a self,
        polysegment: &'a Polysegment,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.inner
            .intersections_polysegment_tol(polysegment, self.epsilon, self.max_relative)
    }

    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.inner
            .intersections_contour_tol(contour, self.epsilon, self.max_relative)
    }

    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a Contour,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.inner
            .intersections_contour_par_tol(contour, self.epsilon, self.max_relative)
    }

    fn intersections_polysegment_par<'a>(
        &'a self,
        polysegment: &'a Polysegment,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.inner
            .intersections_polysegment_par_tol(polysegment, self.epsilon, self.max_relative)
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.inner
            .intersections_shape_tol(shape, self.epsilon, self.max_relative)
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a crate::prelude::Shape,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.inner
            .intersections_shape_par_tol(shape, self.epsilon, self.max_relative)
    }

    fn intersections_point<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.inner
            .intersections_point_tol(point, self.epsilon, self.max_relative)
    }

    fn intersections_point_par<'a>(
        &'a self,
        point: &[f64; 2],
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.inner
            .intersections_point_par_tol(point, self.epsilon, self.max_relative)
    }

    fn intersections_line<'a>(
        &'a self,
        line: &'a crate::prelude::Line,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.inner
            .intersections_line_tol(line, self.epsilon, self.max_relative)
    }

    fn intersections_line_par<'a>(
        &'a self,
        line: &'a crate::prelude::Line,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.inner
            .intersections_line_par_tol(line, self.epsilon, self.max_relative)
    }

    fn intersections_segment<'a, T>(&'a self, segment: T) -> impl Iterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        self.inner
            .intersections_segment_tol(segment, self.epsilon, self.max_relative)
    }

    fn intersections_segment_par<'a, T>(
        &'a self,
        segment: T,
    ) -> impl ParallelIterator<Item = Intersection> + 'a
    where
        T: Into<SegmentRef<'a>>,
    {
        self.inner
            .intersections_segment_par_tol(segment, self.epsilon, self.max_relative)
    }

    fn intersections<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        self.inner
            .intersections_tol(other, self.epsilon, self.max_relative)
    }

    fn intersections_par<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        self.inner
            .intersections_par_tol(other, self.epsilon, self.max_relative)
    }
}

impl CompositeWithTol for Shape {
    fn intersections_polysegment_tol<'a>(
        &'a self,
        polysegment: &'a crate::prelude::Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.contours()
            .iter()
            .enumerate()
            .flat_map(move |(c1, contour_1)| {
                contour_1
                    .polysegment()
                    .intersections_polysegment_tol(polysegment, epsilon, max_relative)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: SegmentKey {
                            contour_idx: c1,
                            segment_idx: i.left.segment_idx,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_polysegment_par_tol<'a>(
        &'a self,
        polysegment: &'a crate::prelude::Polysegment,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.contours()
            .par_iter()
            .enumerate()
            .flat_map(move |(c1, contour_1)| {
                contour_1
                    .polysegment()
                    .intersections_polysegment_par_tol(polysegment, epsilon, max_relative)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: SegmentKey {
                            contour_idx: c1,
                            segment_idx: i.left.segment_idx,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_shape_tol<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> impl Iterator<Item = Intersection> + 'a {
        self.contours()
            .iter()
            .enumerate()
            .flat_map(move |(idx_1, contour_1)| {
                shape
                    .contours()
                    .iter()
                    .enumerate()
                    .map(move |(idx_2, contour_2)| {
                        contour_1
                            .intersections_contour_tol(contour_2, epsilon, max_relative)
                            .map(move |i| Intersection {
                                point: i.point,
                                left: SegmentKey {
                                    contour_idx: idx_1,
                                    segment_idx: i.left.segment_idx,
                                },
                                right: SegmentKey {
                                    contour_idx: idx_2,
                                    segment_idx: i.right.segment_idx,
                                },
                            })
                    })
                    .flatten()
            })
    }

    fn intersections_shape_par_tol<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_relative: f64,
    ) -> impl ParallelIterator<Item = Intersection> + 'a {
        self.contours()
            .par_iter()
            .enumerate()
            .flat_map(move |(idx_1, contour_1)| {
                shape
                    .contours()
                    .par_iter()
                    .enumerate()
                    .map(move |(idx_2, contour_2)| {
                        contour_1
                            .intersections_contour_par_tol(contour_2, epsilon, max_relative)
                            .map(move |i| Intersection {
                                point: i.point,
                                left: SegmentKey {
                                    contour_idx: idx_1,
                                    segment_idx: i.left.segment_idx,
                                },
                                right: SegmentKey {
                                    contour_idx: idx_2,
                                    segment_idx: i.right.segment_idx,
                                },
                            })
                    })
                    .flatten()
            })
    }

    fn contains_point_tol(
        &self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        match self
            .contours()
            .par_iter()
            .enumerate()
            .find_map_any(|(i, c)| {
                if i == 0 {
                    // Special handling of outer contour
                    if let Err(noc) = c.contains_point_tol(point, epsilon, max_relative) {
                        return Some(noc);
                    }
                } else {
                    if let Ok(con) = c.covers_point_tol(point, epsilon, max_relative) {
                        match con {
                            Covered::Inside => return Some(NotContained::InsideHole(i)),
                            Covered::OnBoundary(segment_key) => {
                                return Some(NotContained::OnBoundary(SegmentKey {
                                    contour_idx: i,
                                    segment_idx: segment_key.segment_idx,
                                }));
                            }
                        }
                    }
                }
                return None;
            }) {
            Some(n) => Err(n),
            None => Ok(Contained::Inside),
        }
    }

    fn contains_segment_tol<'a, T: Into<crate::prelude::SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        let segment: SegmentRef = segment.into();
        match self
            .contours()
            .par_iter()
            .enumerate()
            .find_map_any(|(i, c)| {
                if i == 0 {
                    // Special handling of outer contour
                    if let Err(noc) = c.contains_segment_tol(segment.clone(), epsilon, max_relative)
                    {
                        return Some(noc);
                    }
                } else {
                    if let Ok(_) = c.covers_segment_tol(segment.clone(), epsilon, max_relative) {
                        return Some(NotContained::InsideHole(i));
                    }
                }
                return None;
            }) {
            Some(n) => Err(n),
            None => Ok(Contained::Inside),
        }
    }

    fn contains_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Contained, NotContained> {
        match other.into() {
            GeometryRef::Point(pt) => self.contains_point_tol(pt, epsilon, max_relative),
            GeometryRef::ArcSegment(arc_segment) => {
                self.contains_segment_tol(arc_segment, epsilon, max_relative)
            }
            GeometryRef::LineSegment(line_segment) => {
                self.contains_segment_tol(line_segment, epsilon, max_relative)
            }
            GeometryRef::Line(_) => Err(NotContained::OutsideBoundingBox),
            GeometryRef::Segment(segment) => {
                self.contains_segment_tol(segment, epsilon, max_relative)
            }
            GeometryRef::Polysegment(polysegment) => {
                self.contains_polysegment_tol(polysegment, epsilon, max_relative)
            }
            GeometryRef::Contour(contour) => {
                self.contains_contour_tol(contour, epsilon, max_relative)
            }
            GeometryRef::Shape(shape) => self.contains_shape_tol(shape, epsilon, max_relative),
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(bounding_box);
                self.contains_contour_tol(&contour, epsilon, max_relative)
            }
        }
    }

    fn covers_point_tol(
        &self,
        point: &[f64; 2],
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        match self
            .contours()
            .par_iter()
            .enumerate()
            .find_map_any(|(i, c)| {
                if i == 0 {
                    // Special handling of outer contour
                    if let Err(noc) = c.covers_point_tol(point, epsilon, max_relative) {
                        return Some(noc);
                    }
                } else {
                    if let Ok(_) = c.contains_point_tol(point, epsilon, max_relative) {
                        return Some(NotCovered::InsideHole(i));
                    }
                }
                return None;
            }) {
            Some(n) => Err(n),
            None => Ok(Covered::Inside),
        }
    }

    fn covers_segment_tol<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        let segment: SegmentRef = segment.into();

        let result = self
            .contour()
            .covers_segment_tol(segment.clone(), epsilon, max_relative)?;

        for (i, hole) in self.holes().iter().enumerate() {
            if hole
                .contains_any_segment_tol(segment.clone(), epsilon, max_relative)
                .is_ok()
            {
                return Err(NotCovered::InsideHole(i));
            }
        }
        return Ok(result);
    }

    fn contains_any_segment_tol<'a, T: Into<SegmentRef<'a>>>(
        &self,
        segment: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap> {
        let segment: SegmentRef = segment.into();
        /*
        The segment overlaps the shape if it overlaps the outer contour and is
        not contained in one of the shape's holes
         */
        let overlaps =
            self.contour()
                .contains_any_segment_tol(segment.clone(), epsilon, max_relative)?;

        if let Some(i) = self
            .holes()
            .par_iter()
            .enumerate()
            .find_map_any(|(i, hole)| {
                if hole
                    .covers_segment_tol(segment.clone(), epsilon, max_relative)
                    .is_ok()
                {
                    return Some(i);
                }
                return None;
            })
        {
            return Err(NoOverlap::InsideHole(i));
        }

        return Ok(overlaps);
    }

    fn contains_any_contour_tol(
        &self,
        contour: &Contour,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap> {
        // If other either not overlaps the outer contour of the shape or if it
        // is completely covered by one of the holes, then there is no overlap
        let overlaps = self
            .contour()
            .contains_any_contour_tol(contour, epsilon, max_relative)?;

        if let Some(i) = self
            .holes()
            .par_iter()
            .enumerate()
            .find_map_any(|(i, hole)| {
                if hole
                    .covers_contour_tol(contour, epsilon, max_relative)
                    .is_ok()
                {
                    return Some(i);
                }
                return None;
            })
        {
            return Err(NoOverlap::InsideHole(i));
        }

        return Ok(overlaps);
    }

    fn contains_any_shape_tol(
        &self,
        other: &Self,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap> {
        if std::ptr::eq(self, other) {
            // A Shape naturally overlaps itself -> Just return the key to the
            // first segment of self / other (any other segment would also work)
            return Ok(Overlap::Identical);
        }

        if other
            .contains_any_contour_tol(self.contour(), epsilon, max_relative)
            .is_ok()
        {
            return self.contains_any_contour_tol(other.contour(), epsilon, max_relative);
        }
        return Err(NoOverlap::NoPointContained);
    }

    fn contains_any_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Overlap, NoOverlap> {
        match other.into() {
            GeometryRef::Point(pt) => match self.contains_point_tol(pt, epsilon, max_relative) {
                Ok(_) => Ok(Overlap::Point(*pt)),
                Err(_) => Err(NoOverlap::DisjointBoundingBoxes),
            },
            GeometryRef::ArcSegment(arc_segment) => {
                self.contains_any_segment_tol(arc_segment, epsilon, max_relative)
            }
            GeometryRef::LineSegment(line_segment) => {
                self.contains_any_segment_tol(line_segment, epsilon, max_relative)
            }
            GeometryRef::Line(line) => match self.intersections_line(&line).next() {
                Some(pt) => Ok(Overlap::Point(pt.point)),
                None => Err(NoOverlap::DisjointBoundingBoxes),
            },
            GeometryRef::Segment(segment) => {
                self.contains_any_segment_tol(segment, epsilon, max_relative)
            }
            GeometryRef::Polysegment(polysegment) => {
                self.contains_any_polysegment_tol(polysegment, epsilon, max_relative)
            }
            GeometryRef::Contour(contour) => {
                self.contains_any_contour_tol(contour, epsilon, max_relative)
            }
            GeometryRef::Shape(shape) => self.contains_any_shape_tol(shape, epsilon, max_relative),
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(bounding_box);
                self.contains_any_contour_tol(&contour, epsilon, max_relative)
            }
        }
    }

    fn covers_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Result<Covered, NotCovered> {
        match other.into() {
            GeometryRef::Point(pt) => self.covers_point_tol(pt, epsilon, max_relative),
            GeometryRef::ArcSegment(arc_segment) => {
                self.covers_segment_tol(arc_segment, epsilon, max_relative)
            }
            GeometryRef::LineSegment(line_segment) => {
                self.covers_segment_tol(line_segment, epsilon, max_relative)
            }
            GeometryRef::Line(_) => Err(NotCovered::OutsideBoundingBox),
            GeometryRef::Segment(segment) => {
                self.covers_segment_tol(segment, epsilon, max_relative)
            }
            GeometryRef::Polysegment(polysegment) => {
                self.covers_polysegment_tol(polysegment, epsilon, max_relative)
            }
            GeometryRef::Contour(contour) => {
                self.covers_contour_tol(contour, epsilon, max_relative)
            }
            GeometryRef::Shape(shape) => self.covers_shape_tol(shape, epsilon, max_relative),
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(bounding_box);
                self.covers_contour_tol(&contour, epsilon, max_relative)
            }
        }
    }

    fn intersections_tol<'a, T: Into<GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        let geo_self: GeometryRef = self.into();
        let geo_other: GeometryRef = other.into();
        geo_self.intersections_tol::<false, _>(geo_other, epsilon, max_relative)
    }

    fn intersections_par_tol<'a, T: Into<crate::geometry::GeometryRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> Vec<Intersection>
    where
        Self: Sized,
    {
        let geo_self: GeometryRef = self.into();
        let geo_other: GeometryRef = other.into();
        geo_self.intersections_tol::<true, _>(geo_other, epsilon, max_relative)
    }
}

impl TryFrom<Contour> for Shape {
    type Error = ShapeConstructorError<Vec<Contour>>;

    fn try_from(value: Contour) -> Result<Self, Self::Error> {
        return Shape::new(vec![value]);
    }
}

impl TryFrom<Polysegment> for Shape {
    type Error = ShapeConstructorError<Vec<Contour>>;

    fn try_from(value: Polysegment) -> Result<Self, Self::Error> {
        return Shape::new(vec![value.into()]);
    }
}

impl Transformation for Shape {
    fn translate(&mut self, shift: [f64; 2]) -> () {
        self.0
            .as_mut_slice()
            .par_iter_mut()
            .for_each(|polygon| polygon.translate(shift));
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) -> () {
        self.0
            .as_mut_slice()
            .par_iter_mut()
            .for_each(|polygon| polygon.rotate(center, angle));
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        self.0.as_mut_slice().par_iter_mut().for_each(|polygon| {
            polygon.line_reflection(start, stop);
        })
    }

    fn scale(&mut self, factor: f64) -> () {
        self.0
            .as_mut_slice()
            .par_iter_mut()
            .for_each(|polysegment| {
                polysegment.scale(factor);
            })
    }
}

impl ToBoundingBox for Shape {
    fn bounding_box(&self) -> BoundingBox {
        return self.contour().bounding_box();
    }
}

impl From<&Shape> for CentroidData {
    fn from(value: &Shape) -> Self {
        let contour_centroid = CentroidData::from(value.contour());
        return value
            .holes()
            .iter()
            .map(|contour| CentroidData::from(contour))
            .fold(contour_centroid, |prev, curr| prev.subtract(&curr));
    }
}

impl From<Shape> for Vec<Contour> {
    fn from(value: Shape) -> Self {
        return value.0;
    }
}
