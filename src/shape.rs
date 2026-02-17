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

use crate::composite::Composite;
use crate::error::ShapeConstructorError;
use crate::primitive::Primitive;
use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS, composite::*};
use crate::{Transformation, contour::Contour};
use bounding_box::BoundingBox;
use num::Integer;
use rayon::prelude::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/**
A collection of [`Contour`]s which describe the boundaries / extents of a
geometric body in 2D space.

A [`Shape`] always has at least one [`Contour`] which describes the outer
extents of the body represented by it. Any additional contours are interpreted
as "holes" within the "outer" contour. These contours need to fulfill the
following invariants:
- The outer contour contains all holes (see [`Contour::contains`]).
- No hole may contain another hole
- No contour intersects with itself or any other contour.
- No contour might be empty (i.e. contain no
[`Segment`](crate::segment::Segment)).

# Constructing and modifying a shape

The standard constructor which checks all the aforementioned invariants is
[`new`](Shape::new), which takes a [`Vec<Contour`], interprets the first
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
converted into the underlying [`Vec<Contour>`] with the [`From`] implementation.

# Serialization and deserialization

When the `serde` feature is enabled, a shape can be serialized and
deserialized using the [`serde`] crate. It uses the same serialized
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

    */
    #[cfg_attr(
        docsrs,
        doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/example_shape.svg \"Shape\")"
    )]
    #[cfg_attr(
        not(docsrs),
        doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/example_shape.svg)"
    )]
    /**

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
    let c1 = Contour::new(SegmentChain::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_ok());

    // Given vector is empty
    assert!(Shape::new(Vec::new()).is_err());

    // One of the contours is empty
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(SegmentChain::from_points(vertices));
    let c2 = Contour::new(SegmentChain::new());
    assert!(Shape::new(vec![c1, c2]).is_err());

    // The contours intersect
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(SegmentChain::from_points(vertices));
    let vertices = &[[0.1, 0.1], [1.1, 0.1], [1.1, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_err());

    // Second contour is not inside the first contour
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(SegmentChain::from_points(vertices));
    let vertices = &[[1.0, 0.1], [2.0, 0.1], [2.0, 0.9], [1.0, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(vertices));
    assert!(Shape::new(vec![c1, c2]).is_err());
    ```
     */
    pub fn new(
        contours: Vec<Contour>,
    ) -> Result<Self, ShapeConstructorError<Vec<Contour>, ShapeIdx>> {
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
            if !outer.contains(&first_hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                return Err(ShapeConstructorError::HoleOutsideContour {
                    input: this.0,
                    idx: first_hole_idx + 1, // First element of this.0 is the outer contour.
                });
            }

            for (second_hole_idx, second_hole) in
                this.holes()[(first_hole_idx + 1)..].iter().enumerate()
            {
                if first_hole.contains(second_hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                    return Err(ShapeConstructorError::HoleInsideHole {
                        input: this.0,
                        outer_hole_idx: first_hole_idx,
                        inner_hole_idx: second_hole_idx,
                    });
                }
                if second_hole.contains(first_hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                    return Err(ShapeConstructorError::HoleInsideHole {
                        input: this.0,
                        outer_hole_idx: second_hole_idx,
                        inner_hole_idx: first_hole_idx,
                    });
                }
            }
        }

        if let Some(intersection) = this
            .intersections_shape_par(&this, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
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
    let c = Contour::new(SegmentChain::from_points(vertices));
    assert!(Shape::from_outer(c).is_ok());

    // Contour is empty
    let c = Contour::new(SegmentChain::new());
    assert!(Shape::from_outer(c).is_err());

    // Contour intersects itself
    let vertices = &[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]];
    let c = Contour::new(SegmentChain::from_points(vertices));
    assert!(Shape::from_outer(c).is_err());
    ```
     */
    pub fn from_outer(outer: Contour) -> Result<Self, ShapeConstructorError<Contour, SegmentIdx>> {
        if outer.is_empty() {
            return Err(ShapeConstructorError::EmptyContour {
                input: outer,
                idx: 0,
            });
        }
        if let Some(i) = outer
            .intersections_contour_par(&outer, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .find_map_any(|v| Some(v))
        {
            return Err(ShapeConstructorError::Intersection {
                input: outer,
                intersection: Intersection {
                    point: i.point,
                    left: ShapeIdx::new(0usize, i.left),
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
    let c1 = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let c3 = Contour::new(SegmentChain::from_points(vertices));

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
    let contour = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));

    let shape = Shape::new(vec![contour.clone(), hole.clone()]).expect("valid input");
    assert_eq!(shape.area(), contour.area() - hole.area());
    ```
     */
    pub fn area(&self) -> f64 {
        let contour_area = self.contour().area();
        let holes: f64 = self
            .holes()
            .par_iter()
            .map(|segment_chain| return segment_chain.area())
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
    let mut shape = Shape::new(vec![Contour::new(SegmentChain::from_points(vertices))]).expect("valid input");
    assert_eq!(shape.holes().len(), 0);

    // Add a hole to the shape - this works
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));
    assert!(shape.add_hole(hole).is_ok());
    assert_eq!(shape.holes().len(), 1);

    // Adding this hole does not work because it intersects with the outer contour
    let vertices = &[[1.0, 0.1], [1.0, 0.1], [2.0, 0.9], [2.0, 0.9]];
    let hole = Contour::new(SegmentChain::from_points(vertices));
    assert!(shape.add_hole(hole).is_err());
    assert_eq!(shape.holes().len(), 1);
    ```
     */
    pub fn add_hole(
        &mut self,
        hole: Contour,
    ) -> Result<(), ShapeConstructorError<Contour, SegmentIdx>> {
        if hole.is_empty() {
            return Err(ShapeConstructorError::EmptyContour {
                input: hole,
                idx: 0,
            });
        }

        let outer = self.contour();
        if !outer.contains(&hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
            return Err(ShapeConstructorError::HoleOutsideContour {
                input: hole,
                idx: 0, // Dummy value
            });
        }

        for (shape_hole_idx, shape_hole) in self.holes().iter().enumerate() {
            if hole.contains(shape_hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                return Err(ShapeConstructorError::HoleInsideHole {
                    input: hole,
                    outer_hole_idx: shape_hole_idx,
                    inner_hole_idx: 0,
                });
            }
            if shape_hole.contains(&hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
                return Err(ShapeConstructorError::HoleInsideHole {
                    input: hole,
                    outer_hole_idx: 0,
                    inner_hole_idx: shape_hole_idx,
                });
            }
        }

        if let Some(intersection) = self
            .intersections_contour_par(&hole, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
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
    let c1 = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(vertices));

    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let c3 = Contour::new(SegmentChain::from_points(vertices));

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
    Returns an iterator over all [`Contours`]s of `self`.
     */
    pub fn iter(&self) -> std::slice::Iter<'_, Contour> {
        return self.0.iter();
    }

    /**
    Returns a parallel iterator over all [`Contours`]s of `self`.
     */
    pub fn par_iter(&self) -> rayon::slice::Iter<'_, Contour> {
        return self.0.par_iter();
    }
}

impl Composite for Shape {
    type SegmentKey = ShapeIdx;

    fn segment(&self, key: Self::SegmentKey) -> Option<&crate::segment::Segment> {
        let contour = self.0.get(key.contour_idx)?;
        return contour.segment(key.segment_idx);
    }

    fn num_segments(&self) -> usize {
        return self.0.iter().map(Composite::num_segments).sum();
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
        self.contours()
            .iter()
            .enumerate()
            .flat_map(move |(idx_1, contour_1)| {
                contour_1
                    .segment_chain()
                    .intersections_primitive(primitive, epsilon, max_ulps)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: ShapeIdx {
                            contour_idx: idx_1,
                            segment_idx: i.left,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_primitive_par<'a, T: Primitive + std::marker::Sync>(
        &'a self,
        primitive: &'a T,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ()>> + 'a {
        self.contours()
            .par_iter()
            .enumerate()
            .flat_map(move |(idx_1, contour_1)| {
                contour_1
                    .segment_chain()
                    .intersections_primitive_par(primitive, epsilon, max_ulps)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: ShapeIdx {
                            contour_idx: idx_1,
                            segment_idx: i.left,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_segment_chain<'a>(
        &'a self,
        segment_chain: &'a crate::prelude::SegmentChain,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        self.contours()
            .iter()
            .enumerate()
            .flat_map(move |(c1, contour_1)| {
                contour_1
                    .segment_chain()
                    .intersections_segment_chain(segment_chain, epsilon, max_ulps)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: ShapeIdx {
                            contour_idx: c1,
                            segment_idx: i.left,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_segment_chain_par<'a>(
        &'a self,
        segment_chain: &'a crate::prelude::SegmentChain,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        self.contours()
            .par_iter()
            .enumerate()
            .flat_map(move |(c1, contour_1)| {
                contour_1
                    .segment_chain()
                    .intersections_segment_chain_par(segment_chain, epsilon, max_ulps)
                    .map(move |i| Intersection {
                        point: i.point,
                        left: ShapeIdx {
                            contour_idx: c1,
                            segment_idx: i.left,
                        },
                        right: i.right,
                    })
            })
    }

    fn intersections_contour<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        self.intersections_segment_chain(contour.segment_chain(), epsilon, max_ulps)
    }

    fn intersections_contour_par<'a>(
        &'a self,
        contour: &'a Contour,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, SegmentIdx>> + 'a {
        self.intersections_segment_chain_par(contour.segment_chain(), epsilon, max_ulps)
    }

    fn intersections_shape<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl Iterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> + 'a {
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
                            .intersections_contour(contour_2, epsilon, max_ulps)
                            .map(move |i| Intersection {
                                point: i.point,
                                left: ShapeIdx {
                                    contour_idx: idx_1,
                                    segment_idx: i.left,
                                },
                                right: ShapeIdx {
                                    contour_idx: idx_2,
                                    segment_idx: i.right,
                                },
                            })
                    })
                    .flatten()
            })
    }

    fn intersections_shape_par<'a>(
        &'a self,
        shape: &'a Shape,
        epsilon: f64,
        max_ulps: u32,
    ) -> impl ParallelIterator<Item = Intersection<Self::SegmentKey, ShapeIdx>> + 'a {
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
                            .intersections_contour_par(contour_2, epsilon, max_ulps)
                            .map(move |i| Intersection {
                                point: i.point,
                                left: ShapeIdx {
                                    contour_idx: idx_1,
                                    segment_idx: i.left,
                                },
                                right: ShapeIdx {
                                    contour_idx: idx_2,
                                    segment_idx: i.right,
                                },
                            })
                    })
                    .flatten()
            })
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
            .intersections_shape(self, epsilon, max_ulps)
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
            .intersections_shape_par(self, epsilon, max_ulps)
            .map(Intersection::switch);
    }

    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        // First coarse, but fast test: Check if the point is inside the polygon
        let bb = BoundingBox::from(self);
        if !bb.approx_contains_point(point, epsilon, max_ulps) {
            return false;
        }

        // Check if the point is located on the edge of one of the contours
        if self
            .0
            .par_iter()
            .find_any(|contour| {
                contour
                    .segment_chain()
                    .contains_point(point, epsilon, max_ulps)
            })
            .is_some()
        {
            return true;
        }

        // Use the ray casting algorithm
        let ray = crate::line::Line::from_point_angle(point, 0.0);

        let mut counter = 0;
        for contour in &self.0 {
            for segment in contour.segment_chain().iter() {
                counter += segment
                    .intersections_primitive(&ray, epsilon, max_ulps)
                    .len();
            }
        }

        return counter.is_odd();
    }
}

impl From<Shape> for Contour {
    fn from(shape: Shape) -> Contour {
        return shape
            .0
            .into_iter()
            .next()
            .expect("the shape contains at least one line.");
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
            .for_each(|segment_chain| {
                segment_chain.scale(factor);
            })
    }
}

impl From<&Shape> for BoundingBox {
    fn from(value: &Shape) -> BoundingBox {
        return value.contour().into();
    }
}

impl From<&Shape> for crate::CentroidData {
    fn from(value: &Shape) -> Self {
        let contour_centroid = crate::CentroidData::from(value.contour());
        return value
            .holes()
            .iter()
            .map(|contour| crate::CentroidData::from(contour))
            .fold(contour_centroid, |prev, curr| prev.subtract(&curr));
    }
}

impl From<Shape> for Vec<Contour> {
    fn from(value: Shape) -> Self {
        return value.0;
    }
}
