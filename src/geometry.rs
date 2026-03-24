/*!
Container enums for geometric types.

Sometimes, it can be useful to store different geometric types (e.g. a
[`Segment`] and a [`Contour`]) in a collection (e.g. a [`Vec`]). For this purpose,
this module offers the following container enums:
- [`Geometry`] for owned types
- [`GeometryRef`] for borrowed types
- [`GeometryCow`] for types wrapped in a [`Cow`] smart pointer.

A geometric type such as a [`Segment`] can be wrapped into its container via
[`From`]:

```
use bounding_box::BoundingBox;
use planar_geo::prelude::*;

let segment: Segment = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).expect("points not identical").into();
let contour = Contour::from(BoundingBox::new(0.0, 1.0, 2.0, 3.0));
let geometries: &[Geometry] = &[segment.into(), contour.into()];
```

It is also possible to convert between the different containers: For example, a
shared reference to a [`Geometry`] can be converted into a [`GeometryRef`]:

```
use planar_geo::prelude::*;

let segment: Segment = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).expect("points not identical").into();
let geo = Geometry::from(segment);
let geo_ref = GeometryRef::from(&geo);
```

These containers enable the generic `intersections` interface which avoids the
need to use the specialized intersection methods from [`Primitive`] or
[`Composite`] (it defers to those under the hood). Its main disadvantage is that
it needs to eagerly allocate a vector to hold the [`Intersection`]s, because the
specialized methods return different iterator types. If that is not an issues,
using the `intersections` methods greatly simplifies finding the intersections
between any two geometric types

```
use planar_geo::prelude::*;

let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let contour = Contour::new(Polysegment::from_points(vertices));
let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
let hole = Contour::new(Polysegment::from_points(vertices));
let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
let polysegment = Polysegment::from_points(vertices);

let line = Line::from_point_angle([0.0, 0.5], 0.0);

// Generic interface
let intersections: Vec<Intersection> = polysegment.intersections(&shape, e, m);
assert_eq!(intersections.len(), 4);
let intersections: Vec<Intersection> = shape.intersections(&line, e, m);
assert_eq!(intersections.len(), 4);

// Specialized interface (note the explicit allocation from the iterators)
let intersections: Vec<Intersection> = polysegment.intersections_composite(&shape, e, m).collect();
assert_eq!(intersections.len(), 4);
let intersections: Vec<Intersection> = shape.intersections_primitive(&line, e, m).collect();
assert_eq!(intersections.len(), 4);
```

Another use for these containers is to find the common bounding box of multiple
different geometric types (the containers implement
[`bounding_box::ToBoundingBox`]):

```
use bounding_box::BoundingBox;
use planar_geo::prelude::*;

let segment: Segment = LineSegment::new([0.0, 0.0], [2.0, 0.0], 0.0, 0).expect("points not identical").into();
let contour = Contour::from(BoundingBox::new(0.0, 1.0, 2.0, 3.0));
let geometries: &[Geometry] = &[segment.into(), contour.into()];

let bb = BoundingBox::from_bounded_entities(geometries.iter()).expect("at least one item");
assert_eq!(bb.xmin(), 0.0);
assert_eq!(bb.xmax(), 2.0);
assert_eq!(bb.ymin(), 0.0);
assert_eq!(bb.ymax(), 3.0);
```
If the `cairo` feature is activated, it is also possible to use the container
types for drawing onto a cairo canvas.
*/
#![cfg_attr(
    feature = "cairo",
    cfg_attr(
        all(),
        doc = r#"See the module documentation of [drawable](crate::draw::drawable) for more."#,
    )
)]

use std::borrow::Cow;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use bounding_box::BoundingBox;

use rayon::prelude::*;

use crate::Transformation;
use crate::composite::{Composite, Intersection, SegmentKey};
use crate::contour::Contour;
use crate::line::Line;
use crate::polysegment::Polysegment;
use crate::primitive::Primitive;
use crate::segment::{ArcSegment, LineSegment, Segment};
use crate::shape::Shape;

/**
A container enum for owned geometric types. See the
[module-level documentation](crate::geometry) for more.
 */
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Geometry {
    /// An owned point (`[f64; 2]`).
    Point([f64; 2]),
    /// An owned [`BoundingBox`].
    BoundingBox(BoundingBox),
    /// An owned [`ArcSegment`].
    ArcSegment(ArcSegment),
    /// An owned [`LineSegment`].
    LineSegment(LineSegment),
    /// An owned [`Line`].
    Line(Line),
    /// An owned [`Segment`].
    Segment(Segment),
    /// An owned [`Polysegment`].
    Polysegment(Polysegment),
    /// An owned [`Contour`].
    Contour(Contour),
    /// An owned [`Shape`].
    Shape(Shape),
}

impl Geometry {
    /**
    If `self` holds a [`Composite`], returns a reference to the [`Segment`]
    specified by the given [`SegmentKey`]. Otherwise, returns `None`.
     */
    pub fn segment(&self, key: SegmentKey) -> Option<&Segment> {
        match self {
            Geometry::Point(_) => None,
            Geometry::BoundingBox(_) => None,
            Geometry::ArcSegment(_) => None,
            Geometry::LineSegment(_) => None,
            Geometry::Line(_) => None,
            Geometry::Segment(_) => None,
            Geometry::Polysegment(polysegment) => polysegment.segment(key),
            Geometry::Contour(contour) => contour.segment(key),
            Geometry::Shape(shape) => shape.segment(key),
        }
    }

    /**
    Returns all intersections between `self` and `other`.

    See [`GeometryRef::intersections`] for details and examples.
     */
    pub fn intersections<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let this: GeometryRef = self.into();
        this.intersections(other, epsilon, max_ulps)
    }

    /**
    Returns all intersections between `self` and `other`.

    This is a parallelized version of [`Geometry::intersections`], see its
    docstring for details. It uses parallel variants of the specialized
    intersection algorithms, if available.
     */
    pub fn intersections_par<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let this: GeometryRef = self.into();
        this.intersections_par(other, epsilon, max_ulps)
    }
}

impl Transformation for Geometry {
    fn translate(&mut self, shift: [f64; 2]) {
        match self {
            Geometry::Point(elem) => elem.translate(shift),
            Geometry::BoundingBox(elem) => elem.translate(shift),
            Geometry::ArcSegment(elem) => elem.translate(shift),
            Geometry::LineSegment(elem) => elem.translate(shift),
            Geometry::Line(elem) => elem.translate(shift),
            Geometry::Segment(elem) => elem.translate(shift),
            Geometry::Polysegment(elem) => elem.translate(shift),
            Geometry::Contour(elem) => elem.translate(shift),
            Geometry::Shape(elem) => elem.translate(shift),
        }
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        match self {
            Geometry::Point(elem) => elem.rotate(center, angle),
            Geometry::BoundingBox(elem) => elem.rotate(center, angle),
            Geometry::ArcSegment(elem) => elem.rotate(center, angle),
            Geometry::LineSegment(elem) => elem.rotate(center, angle),
            Geometry::Line(elem) => elem.rotate(center, angle),
            Geometry::Segment(elem) => elem.rotate(center, angle),
            Geometry::Polysegment(elem) => elem.rotate(center, angle),
            Geometry::Contour(elem) => elem.rotate(center, angle),
            Geometry::Shape(elem) => elem.rotate(center, angle),
        }
    }

    fn scale(&mut self, factor: f64) {
        match self {
            Geometry::Point(elem) => elem.scale(factor),
            Geometry::BoundingBox(elem) => elem.scale(factor),
            Geometry::ArcSegment(elem) => elem.scale(factor),
            Geometry::LineSegment(elem) => elem.scale(factor),
            Geometry::Line(elem) => elem.scale(factor),
            Geometry::Segment(elem) => elem.scale(factor),
            Geometry::Polysegment(elem) => elem.scale(factor),
            Geometry::Contour(elem) => elem.scale(factor),
            Geometry::Shape(elem) => elem.scale(factor),
        }
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        match self {
            Geometry::Point(elem) => elem.line_reflection(start, stop),
            Geometry::BoundingBox(elem) => elem.line_reflection(start, stop),
            Geometry::ArcSegment(elem) => elem.line_reflection(start, stop),
            Geometry::LineSegment(elem) => elem.line_reflection(start, stop),
            Geometry::Line(elem) => elem.line_reflection(start, stop),
            Geometry::Segment(elem) => elem.line_reflection(start, stop),
            Geometry::Polysegment(elem) => elem.line_reflection(start, stop),
            Geometry::Contour(elem) => elem.line_reflection(start, stop),
            Geometry::Shape(elem) => elem.line_reflection(start, stop),
        }
    }
}

impl From<&Geometry> for BoundingBox {
    fn from(value: &Geometry) -> Self {
        let geom_ref = GeometryRef::from(value);
        return BoundingBox::from(geom_ref);
    }
}

impl From<[f64; 2]> for Geometry {
    fn from(value: [f64; 2]) -> Self {
        Geometry::Point(value)
    }
}

impl From<BoundingBox> for Geometry {
    fn from(value: BoundingBox) -> Self {
        Geometry::BoundingBox(value)
    }
}

impl From<ArcSegment> for Geometry {
    fn from(value: ArcSegment) -> Self {
        Geometry::ArcSegment(value)
    }
}

impl From<LineSegment> for Geometry {
    fn from(value: LineSegment) -> Self {
        Geometry::LineSegment(value)
    }
}

impl From<Line> for Geometry {
    fn from(value: Line) -> Self {
        Geometry::Line(value)
    }
}

impl From<Segment> for Geometry {
    fn from(value: Segment) -> Self {
        Geometry::Segment(value)
    }
}

impl From<Polysegment> for Geometry {
    fn from(value: Polysegment) -> Self {
        Geometry::Polysegment(value)
    }
}

impl From<Contour> for Geometry {
    fn from(value: Contour) -> Self {
        Geometry::Contour(value)
    }
}

impl From<Shape> for Geometry {
    fn from(value: Shape) -> Self {
        Geometry::Shape(value)
    }
}

/**
A container enum for borrowed geometric types. See the
[module-level documentation](crate::geometry) for more.
 */
#[derive(Debug, Clone)]
pub enum GeometryRef<'a> {
    /// A borrowed point (`[f64; 2]`).
    Point(&'a [f64; 2]),
    /// A borrowed [`BoundingBox`].
    BoundingBox(&'a BoundingBox),
    /// A borrowed [`ArcSegment`].
    ArcSegment(&'a ArcSegment),
    /// A borrowed [`LineSegment`].
    LineSegment(&'a LineSegment),
    /// A borrowed [`Line`].
    Line(&'a Line),
    /// A borrowed [`Segment`].
    Segment(&'a Segment),
    /// A borrowed [`Polysegment`].
    Polysegment(&'a Polysegment),
    /// A borrowed [`Contour`].
    Contour(&'a Contour),
    /// A borrowed [`Shape`].
    Shape(&'a Shape),
}

impl<'a> GeometryRef<'a> {
    /**
    If `self` holds a [`Composite`], returns a reference to the [`Segment`]
    specified by the given [`SegmentKey`]. Otherwise, returns `None`.
     */
    pub fn segment(&self, key: SegmentKey) -> Option<&Segment> {
        match self {
            GeometryRef::Point(_) => None,
            GeometryRef::BoundingBox(_) => None,
            GeometryRef::ArcSegment(_) => None,
            GeometryRef::LineSegment(_) => None,
            GeometryRef::Line(_) => None,
            GeometryRef::Segment(_) => None,
            GeometryRef::Polysegment(polysegment) => polysegment.segment(key),
            GeometryRef::Contour(contour) => contour.segment(key),
            GeometryRef::Shape(shape) => shape.segment(key),
        }
    }

    /**
    Returns all intersections between `self` and `other`.

    This function uses the specialized intersection functions provided by the
    [`Primitive`] and [`Composite`] traits, depending on the underlying types of
    `self` and `other`:
    - Both are [`Primitive`]s: [`Primitive::intersections_primitive`].
    - One of them is a [`Primitive`], the other one is a [`Composite`]:
    [`Composite::intersections_primitive`].
    - Both are [`Composite`]s: [`Composite::intersections_composite`].
    - One of them is a point [`f64; 2`]: [`Primitive::contains_point`] or
    [`Composite::contains_point`].
    - A bounding box is converted to a [`Contour`] ([`Composite`]), then one of
    the functions listed above is used.

    Since all of these underlying functions return different types (iterators,
    enums, booleans), this function unifies the outputs into an allocated
    vector. If allocation is undesirable, consider using the specialized
    functions instead.

    # Examples

    ```
    use planar_geo::prelude::*;

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));
    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

    let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
    let polysegment = Polysegment::from_points(vertices);

    let line = Line::from_point_angle([0.0, 0.5], 0.0);

    // Using GeometryRef
    let intersections: Vec<Intersection> = GeometryRef::from(&polysegment).intersections(&shape, e, m);
    assert_eq!(intersections.len(), 4);

    // Since all of the geometric types in this crate provide an intersection
    // method themselves, it is not necesssary to convert them into GeometryRef
    // explicitly (the conversion happens under the hood)
    let intersections: Vec<Intersection> = shape.intersections(&line, e, m);
    assert_eq!(intersections.len(), 4);

    // Parallelized variants
    let intersections: Vec<Intersection> = polysegment.intersections_par(&shape, e, m);
    assert_eq!(intersections.len(), 4);
    let intersections: Vec<Intersection> = shape.intersections_par(&line, e, m);
    assert_eq!(intersections.len(), 4);
    ```
     */
    pub fn intersections<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let geo_ref: GeometryRef = other.into();
        match self {
            GeometryRef::Point(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::BoundingBox(bounding_box) => {
                geo_ref.intersections_composite(&Contour::from(*bounding_box), epsilon, max_ulps)
            }
            GeometryRef::ArcSegment(elem) => {
                geo_ref.intersections_primitive(*elem, epsilon, max_ulps)
            }
            GeometryRef::LineSegment(elem) => {
                geo_ref.intersections_primitive(*elem, epsilon, max_ulps)
            }
            GeometryRef::Line(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::Segment(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::Polysegment(elem) => {
                geo_ref.intersections_composite(*elem, epsilon, max_ulps)
            }
            GeometryRef::Contour(elem) => geo_ref.intersections_composite(*elem, epsilon, max_ulps),
            GeometryRef::Shape(elem) => geo_ref.intersections_composite(*elem, epsilon, max_ulps),
        }
    }

    /**
    Returns all intersections between `self` and `other`.

    This is a parallelized version of [`GeometryRef::intersections`], see its
    docstring for details. It uses parallel variants of the specialized
    intersection algorithms, if available.
     */
    pub fn intersections_par<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let geo_ref: GeometryRef = other.into();
        match self {
            GeometryRef::Point(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::BoundingBox(bounding_box) => geo_ref.intersections_composite_par(
                &Contour::from(*bounding_box),
                epsilon,
                max_ulps,
            ),
            GeometryRef::ArcSegment(elem) => {
                geo_ref.intersections_primitive(*elem, epsilon, max_ulps)
            }
            GeometryRef::LineSegment(elem) => {
                geo_ref.intersections_primitive(*elem, epsilon, max_ulps)
            }
            GeometryRef::Line(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::Segment(elem) => geo_ref.intersections_primitive(*elem, epsilon, max_ulps),
            GeometryRef::Polysegment(elem) => {
                geo_ref.intersections_composite_par(*elem, epsilon, max_ulps)
            }
            GeometryRef::Contour(elem) => {
                geo_ref.intersections_composite_par(*elem, epsilon, max_ulps)
            }
            GeometryRef::Shape(elem) => {
                geo_ref.intersections_composite_par(*elem, epsilon, max_ulps)
            }
        }
    }

    pub(crate) fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        match self {
            GeometryRef::Point(point) => {
                if other.contains_point((*point).clone(), epsilon, max_ulps) {
                    return vec![(**point).into()];
                } else {
                    return Vec::new();
                }
            }
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(*bounding_box);
                contour
                    .intersections_primitive(other, epsilon, max_ulps)
                    .collect()
            }
            GeometryRef::ArcSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::LineSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::Line(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::Segment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .into_iter()
                .map(From::from)
                .collect(),
            GeometryRef::Polysegment(elem) => elem
                .intersections_primitive(other, epsilon, max_ulps)
                .collect(),
            GeometryRef::Contour(elem) => elem
                .intersections_primitive(other, epsilon, max_ulps)
                .collect(),
            GeometryRef::Shape(elem) => elem
                .intersections_primitive(other, epsilon, max_ulps)
                .collect(),
        }
    }

    pub(crate) fn intersections_composite<T: Composite>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        match self {
            GeometryRef::Point(point) => {
                if other.contains_point((*point).clone(), epsilon, max_ulps) {
                    return vec![(**point).into()];
                } else {
                    return Vec::new();
                }
            }
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(*bounding_box);
                contour
                    .intersections_composite(other, epsilon, max_ulps)
                    .collect()
            }
            GeometryRef::ArcSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::LineSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Line(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Segment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Polysegment(elem) => other
                .intersections_composite(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Contour(elem) => other
                .intersections_composite(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Shape(elem) => other
                .intersections_composite(*elem, epsilon, max_ulps)
                .collect(),
        }
    }

    pub(crate) fn intersections_composite_par<T: Composite>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        match self {
            GeometryRef::Point(point) => {
                if other.contains_point((*point).clone(), epsilon, max_ulps) {
                    return vec![(**point).into()];
                } else {
                    return Vec::new();
                }
            }
            GeometryRef::BoundingBox(bounding_box) => {
                let contour = Contour::from(*bounding_box);
                contour
                    .intersections_composite_par(other, epsilon, max_ulps)
                    .collect()
            }
            GeometryRef::ArcSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::LineSegment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Line(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Segment(elem) => other
                .intersections_primitive(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Polysegment(elem) => other
                .intersections_composite_par(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Contour(elem) => other
                .intersections_composite_par(*elem, epsilon, max_ulps)
                .collect(),
            GeometryRef::Shape(elem) => other
                .intersections_composite_par(*elem, epsilon, max_ulps)
                .collect(),
        }
    }
}

impl<'a> From<&'a Geometry> for GeometryRef<'a> {
    fn from(value: &'a Geometry) -> Self {
        match value {
            Geometry::Point(elem) => GeometryRef::Point(elem),
            Geometry::BoundingBox(elem) => GeometryRef::BoundingBox(elem),
            Geometry::ArcSegment(elem) => GeometryRef::ArcSegment(elem),
            Geometry::LineSegment(elem) => GeometryRef::LineSegment(elem),
            Geometry::Line(elem) => GeometryRef::Line(elem),
            Geometry::Segment(elem) => GeometryRef::Segment(elem),
            Geometry::Polysegment(elem) => GeometryRef::Polysegment(elem),
            Geometry::Contour(elem) => GeometryRef::Contour(elem),
            Geometry::Shape(elem) => GeometryRef::Shape(elem),
        }
    }
}

impl<'a> From<GeometryRef<'a>> for BoundingBox {
    fn from(value: GeometryRef<'a>) -> Self {
        match value {
            GeometryRef::Point(elem) => elem.into(),
            GeometryRef::BoundingBox(elem) => elem.clone(),
            GeometryRef::ArcSegment(elem) => elem.into(),
            GeometryRef::LineSegment(elem) => elem.into(),
            GeometryRef::Line(elem) => elem.into(),
            GeometryRef::Segment(elem) => elem.into(),
            GeometryRef::Polysegment(elem) => elem.into(),
            GeometryRef::Contour(elem) => elem.into(),
            GeometryRef::Shape(elem) => elem.into(),
        }
    }
}

impl<'a> From<&'a [f64; 2]> for GeometryRef<'a> {
    fn from(value: &'a [f64; 2]) -> Self {
        GeometryRef::Point(value)
    }
}

impl<'a> From<&'a BoundingBox> for GeometryRef<'a> {
    fn from(value: &'a BoundingBox) -> Self {
        GeometryRef::BoundingBox(value)
    }
}

impl<'a> From<&'a ArcSegment> for GeometryRef<'a> {
    fn from(value: &'a ArcSegment) -> Self {
        GeometryRef::ArcSegment(value)
    }
}

impl<'a> From<&'a LineSegment> for GeometryRef<'a> {
    fn from(value: &'a LineSegment) -> Self {
        GeometryRef::LineSegment(value)
    }
}

impl<'a> From<&'a Line> for GeometryRef<'a> {
    fn from(value: &'a Line) -> Self {
        GeometryRef::Line(value)
    }
}

impl<'a> From<&'a Segment> for GeometryRef<'a> {
    fn from(value: &'a Segment) -> Self {
        GeometryRef::Segment(value)
    }
}

impl<'a> From<&'a Polysegment> for GeometryRef<'a> {
    fn from(value: &'a Polysegment) -> Self {
        GeometryRef::Polysegment(value)
    }
}

impl<'a> From<&'a Contour> for GeometryRef<'a> {
    fn from(value: &'a Contour) -> Self {
        GeometryRef::Contour(value)
    }
}

impl<'a> From<&'a Shape> for GeometryRef<'a> {
    fn from(value: &'a Shape) -> Self {
        GeometryRef::Shape(value)
    }
}

/**
A container enum for geometric types wrapped in a [`Cow`] smart pointer. See
the [module-level documentation](crate::geometry) for more.
 */
#[derive(Debug, Clone)]
pub enum GeometryCow<'a> {
    /// A [`Cow`]-wrapped point (`[f64; 2]`).
    Point(Cow<'a, [f64; 2]>),
    /// A [`Cow`]-wrapped [`BoundingBox`].
    BoundingBox(Cow<'a, BoundingBox>),
    /// A [`Cow`]-wrapped [`ArcSegment`].
    ArcSegment(Cow<'a, ArcSegment>),
    /// A [`Cow`]-wrapped [`LineSegment`].
    LineSegment(Cow<'a, LineSegment>),
    /// A [`Cow`]-wrapped [`Line`].
    Line(Cow<'a, Line>),
    /// A [`Cow`]-wrapped [`Segment`].
    Segment(Cow<'a, Segment>),
    /// A [`Cow`]-wrapped [`Polysegment`].
    Polysegment(Cow<'a, Polysegment>),
    /// A [`Cow`]-wrapped [`Contour`].
    Contour(Cow<'a, Contour>),
    /// A [`Cow`]-wrapped [`Shape`].
    Shape(Cow<'a, Shape>),
}

impl<'a> GeometryCow<'a> {
    /**
    If `self` holds a [`Composite`], returns a reference to the [`Segment`]
    specified by the given [`SegmentKey`]. Otherwise, returns `None`.
     */
    pub fn segment(&self, key: SegmentKey) -> Option<&Segment> {
        match self {
            GeometryCow::Point(_) => None,
            GeometryCow::BoundingBox(_) => None,
            GeometryCow::ArcSegment(_) => None,
            GeometryCow::LineSegment(_) => None,
            GeometryCow::Line(_) => None,
            GeometryCow::Segment(_) => None,
            GeometryCow::Polysegment(polysegment) => polysegment.segment(key),
            GeometryCow::Contour(contour) => contour.segment(key),
            GeometryCow::Shape(shape) => shape.segment(key),
        }
    }

    /**
    Returns all intersections between `self` and `other`.

    See [`GeometryRef::intersections`] for details and examples.
     */
    pub fn intersections<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let this: GeometryRef = self.into();
        this.intersections(other, epsilon, max_ulps)
    }

    /**
    Returns all intersections between `self` and `other`.

    This is a parallelized version of [`GeometryCow::intersections`], see its
    docstring for details. It uses parallel variants of the specialized
    intersection algorithms, if available.
     */
    pub fn intersections_par<'b, T: Into<GeometryRef<'b>>>(
        &self,
        other: T,
        epsilon: f64,
        max_ulps: u32,
    ) -> Vec<Intersection> {
        let this: GeometryRef = self.into();
        this.intersections_par(other, epsilon, max_ulps)
    }
}

impl<'a> From<Geometry> for GeometryCow<'a> {
    fn from(value: Geometry) -> Self {
        match value {
            Geometry::Point(elem) => GeometryCow::Point(Cow::Owned(elem)),
            Geometry::BoundingBox(elem) => GeometryCow::BoundingBox(Cow::Owned(elem)),
            Geometry::ArcSegment(elem) => GeometryCow::ArcSegment(Cow::Owned(elem)),
            Geometry::LineSegment(elem) => GeometryCow::LineSegment(Cow::Owned(elem)),
            Geometry::Line(elem) => GeometryCow::Line(Cow::Owned(elem)),
            Geometry::Segment(elem) => GeometryCow::Segment(Cow::Owned(elem)),
            Geometry::Polysegment(elem) => GeometryCow::Polysegment(Cow::Owned(elem)),
            Geometry::Contour(elem) => GeometryCow::Contour(Cow::Owned(elem)),
            Geometry::Shape(elem) => GeometryCow::Shape(Cow::Owned(elem)),
        }
    }
}

impl<'a> From<&'a Geometry> for GeometryCow<'a> {
    fn from(value: &'a Geometry) -> Self {
        match value {
            Geometry::Point(elem) => GeometryCow::Point(Cow::Borrowed(elem)),
            Geometry::BoundingBox(elem) => GeometryCow::BoundingBox(Cow::Borrowed(elem)),
            Geometry::ArcSegment(elem) => GeometryCow::ArcSegment(Cow::Borrowed(elem)),
            Geometry::LineSegment(elem) => GeometryCow::LineSegment(Cow::Borrowed(elem)),
            Geometry::Line(elem) => GeometryCow::Line(Cow::Borrowed(elem)),
            Geometry::Segment(elem) => GeometryCow::Segment(Cow::Borrowed(elem)),
            Geometry::Polysegment(elem) => GeometryCow::Polysegment(Cow::Borrowed(elem)),
            Geometry::Contour(elem) => GeometryCow::Contour(Cow::Borrowed(elem)),
            Geometry::Shape(elem) => GeometryCow::Shape(Cow::Borrowed(elem)),
        }
    }
}

impl<'a> From<&'a GeometryCow<'a>> for GeometryRef<'a> {
    fn from(value: &'a GeometryCow<'a>) -> Self {
        match value {
            GeometryCow::Point(elem) => GeometryRef::Point(elem.as_ref()),
            GeometryCow::BoundingBox(elem) => GeometryRef::BoundingBox(elem.as_ref()),
            GeometryCow::ArcSegment(elem) => GeometryRef::ArcSegment(elem.as_ref()),
            GeometryCow::LineSegment(elem) => GeometryRef::LineSegment(elem.as_ref()),
            GeometryCow::Line(elem) => GeometryRef::Line(elem.as_ref()),
            GeometryCow::Segment(elem) => GeometryRef::Segment(elem.as_ref()),
            GeometryCow::Polysegment(elem) => GeometryRef::Polysegment(elem.as_ref()),
            GeometryCow::Contour(elem) => GeometryRef::Contour(elem.as_ref()),
            GeometryCow::Shape(elem) => GeometryRef::Shape(elem.as_ref()),
        }
    }
}

impl<'a> From<GeometryCow<'a>> for BoundingBox {
    fn from(value: GeometryCow<'a>) -> Self {
        (&value).into()
    }
}

impl<'a> From<&'a GeometryCow<'a>> for BoundingBox {
    fn from(value: &'a GeometryCow<'a>) -> Self {
        GeometryRef::from(value).into()
    }
}

impl From<[f64; 2]> for GeometryCow<'_> {
    fn from(value: [f64; 2]) -> Self {
        GeometryCow::Point(Cow::Owned(value))
    }
}

impl From<BoundingBox> for GeometryCow<'_> {
    fn from(value: BoundingBox) -> Self {
        GeometryCow::BoundingBox(Cow::Owned(value))
    }
}

impl From<ArcSegment> for GeometryCow<'_> {
    fn from(value: ArcSegment) -> Self {
        GeometryCow::ArcSegment(Cow::Owned(value))
    }
}

impl From<LineSegment> for GeometryCow<'_> {
    fn from(value: LineSegment) -> Self {
        GeometryCow::LineSegment(Cow::Owned(value))
    }
}

impl From<Line> for GeometryCow<'_> {
    fn from(value: Line) -> Self {
        GeometryCow::Line(Cow::Owned(value))
    }
}

impl From<Segment> for GeometryCow<'_> {
    fn from(value: Segment) -> Self {
        GeometryCow::Segment(Cow::Owned(value))
    }
}

impl From<Polysegment> for GeometryCow<'_> {
    fn from(value: Polysegment) -> Self {
        GeometryCow::Polysegment(Cow::Owned(value))
    }
}

impl From<Contour> for GeometryCow<'_> {
    fn from(value: Contour) -> Self {
        GeometryCow::Contour(Cow::Owned(value))
    }
}

impl From<Shape> for GeometryCow<'_> {
    fn from(value: Shape) -> Self {
        GeometryCow::Shape(Cow::Owned(value))
    }
}

impl<'a> From<&'a [f64; 2]> for GeometryCow<'a> {
    fn from(value: &'a [f64; 2]) -> Self {
        GeometryCow::Point(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a BoundingBox> for GeometryCow<'a> {
    fn from(value: &'a BoundingBox) -> Self {
        GeometryCow::BoundingBox(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a ArcSegment> for GeometryCow<'a> {
    fn from(value: &'a ArcSegment) -> Self {
        GeometryCow::ArcSegment(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a LineSegment> for GeometryCow<'a> {
    fn from(value: &'a LineSegment) -> Self {
        GeometryCow::LineSegment(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a Line> for GeometryCow<'a> {
    fn from(value: &'a Line) -> Self {
        GeometryCow::Line(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a Segment> for GeometryCow<'a> {
    fn from(value: &'a Segment) -> Self {
        GeometryCow::Segment(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a Polysegment> for GeometryCow<'a> {
    fn from(value: &'a Polysegment) -> Self {
        GeometryCow::Polysegment(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a Contour> for GeometryCow<'a> {
    fn from(value: &'a Contour) -> Self {
        GeometryCow::Contour(Cow::Borrowed(value))
    }
}

impl<'a> From<&'a Shape> for GeometryCow<'a> {
    fn from(value: &'a Shape) -> Self {
        GeometryCow::Shape(Cow::Borrowed(value))
    }
}
