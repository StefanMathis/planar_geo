/*!
This module contains the structs [`Drawable`], [`DrawableRef`] and
[`DrawableCow`] - helper structs for drawing a [`Geometry`], [`GeometryRef`] or
[`GeometryCow`] respectively.

Besides the aforementioned enums, each of these structs holds a [`Style`] and
implements a `draw` method: [`Drawable::draw`], [`DrawableRef::draw`] and
[`DrawableCow::draw`]. The `draw` method calls the `draw` method of the
contained geometric type (e.g.
[`ArcSegment::draw`](crate::segment::ArcSegment::draw)) with its [`Style`] field
as the second argument. The helper structs are therefore useful to group a
geometric type together with its [`Style`] - for example, to pass them around.

There are two geometric types which cannot be drawn:
- The "point" type `[f64; 2]` (because it has no size)
- The [`Line`](crate::line::Line) (because it has an infinite length and cannot
be drawn onto a finite canvas)
For these types, the `draw` method simply returns `Ok(())`
 */

use bounding_box::{BoundingBox, ToBoundingBox};

use super::Style;
use crate::{Transformation, geometry::*};

/**
A wrapper around a [`Geometry`] and its corresponding [`Style`].

It provides the [`Drawable::draw`] method for drawing itself onto a
[`cairo::Context`]. See the [module level documentation](crate::draw::drawable)
for more.

A special feature of [`Geometry`] compared to [`GeometryRef`] and
[`GeometryCow`] is the fact that it owns its [`Drawable::geometry`]. Therefore,
it implements [`Transformation`].
 */
#[derive(Clone)]
pub struct Drawable {
    /// Underlying geometry object.
    pub geometry: Geometry,
    /// [`Style`] for the [`Drawable::geometry`], used in [`Drawable::draw`].
    pub style: Style,
}

impl Drawable {
    /**
    Returns a new instance of a [`Drawable`] from any [`Geometry`] and its
    [`Style`].
     */
    pub fn new<G: Into<Geometry>>(geo: G, style: Style) -> Self {
        Drawable {
            geometry: geo.into(),
            style,
        }
    }

    /**
    Draws the wrapped [`Geometry`] onto the `context` using the [`Style`]
    provided by the [`Drawable::style`] field.

    See the [module level documentation](crate::draw::drawable) for more.
     */
    pub fn draw(&self, context: &cairo::Context) -> Result<(), cairo::Error> {
        DrawableRef::from(self).draw(context)
    }
}

impl Transformation for Drawable {
    fn translate(&mut self, shift: [f64; 2]) {
        self.geometry.translate(shift);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        self.geometry.rotate(center, angle);
    }

    fn scale(&mut self, factor: f64) {
        self.geometry.scale(factor);
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        self.geometry.line_reflection(start, stop);
    }
}

impl ToBoundingBox for Drawable {
    fn bounding_box(&self) -> BoundingBox {
        (&self.geometry).bounding_box()
    }
}

impl<I> From<(I, Style)> for Drawable
where
    I: Into<Geometry>,
{
    fn from(value: (I, Style)) -> Self {
        return Drawable {
            geometry: value.0.into(),
            style: value.1,
        };
    }
}

/**
A wrapper around a [`GeometryRef`] and its corresponding [`Style`].

It provides the [`Drawable::draw`] method for drawing itself onto a
[`cairo::Context`]. See the [module level documentation](crate::draw::drawable)
for more.
 */
#[derive(Clone)]
pub struct DrawableRef<'a> {
    /// Underlying geometry object.
    pub geometry: GeometryRef<'a>,
    /// [`Style`] for the [`DrawableRef::geometry`], used in [`Drawable::draw`].
    pub style: Style,
}

impl<'a> DrawableRef<'a> {
    /**
    Returns a new instance of a [`DrawableRef`] from any [`GeometryRef`] and its
    [`Style`].
     */
    pub fn new<G: Into<GeometryRef<'a>>>(geo: G, style: Style) -> Self {
        DrawableRef {
            geometry: geo.into(),
            style,
        }
    }

    /**
    Draws the wrapped [`GeometryRef`] onto the `context` using the [`Style`]
    provided by the [`DrawableRef::style`] field.

    See the [module level documentation](crate::draw::drawable) for more.
     */
    pub fn draw(&self, context: &cairo::Context) -> Result<(), cairo::Error> {
        match self.geometry {
            GeometryRef::Point(_) => Ok(()), // Points cannot be drawn
            GeometryRef::BoundingBox(elem) => {
                crate::contour::Contour::from(elem).draw(&self.style, context)
            }
            GeometryRef::ArcSegment(elem) => elem.draw(&self.style, context),
            GeometryRef::LineSegment(elem) => elem.draw(&self.style, context),
            GeometryRef::Line(_) => Ok(()), // Lines cannot be drawn
            GeometryRef::Segment(elem) => elem.draw(&self.style, context),
            GeometryRef::Polysegment(elem) => elem.draw(&self.style, context),
            GeometryRef::Contour(elem) => elem.draw(&self.style, context),
            GeometryRef::Shape(elem) => elem.draw(&self.style, context),
        }
    }
}

impl<'a> ToBoundingBox for DrawableRef<'a> {
    fn bounding_box(&self) -> BoundingBox {
        return self.geometry.bounding_box();
    }
}

impl<'a> From<DrawableRef<'a>> for Drawable {
    fn from(value: DrawableRef<'a>) -> Self {
        return Drawable {
            geometry: value.geometry.into(),
            style: value.style,
        };
    }
}

impl<'a> From<&'a Drawable> for DrawableRef<'a> {
    fn from(value: &'a Drawable) -> Self {
        return DrawableRef {
            geometry: (&value.geometry).into(),
            style: value.style.clone(),
        };
    }
}

impl<'a, I> From<(I, Style)> for DrawableRef<'a>
where
    I: Into<GeometryRef<'a>>,
{
    fn from(value: (I, Style)) -> Self {
        return DrawableRef {
            geometry: value.0.into(),
            style: value.1,
        };
    }
}

/**
A wrapper around a [`GeometryCow`] and its corresponding [`Style`].

It provides the [`Drawable::draw`] method for drawing itself onto a
[`cairo::Context`]. See the [module level documentation](crate::draw::drawable)
for more.
 */
#[derive(Clone)]
pub struct DrawableCow<'a> {
    /// Underlying geometry object.
    pub geometry: GeometryCow<'a>,
    /// [`Style`] for the [`DrawableCow::geometry`], used in [`Drawable::draw`].
    pub style: Style,
}

impl<'a> DrawableCow<'a> {
    /**
    Returns a new instance of a [`GeometryCow`] from any [`GeometryCow`] and its
    [`Style`].
     */
    pub fn new<G: Into<GeometryCow<'a>>>(geo: G, style: Style) -> Self {
        DrawableCow {
            geometry: geo.into(),
            style,
        }
    }

    /**
    Draws the wrapped [`GeometryCow`] onto the `context` using the [`Style`]
    provided by the [`DrawableCow::style`] field.

    See the [module level documentation](crate::draw::drawable) for more.
     */
    pub fn draw(&self, context: &cairo::Context) -> Result<(), cairo::Error> {
        DrawableRef::from(self).draw(context)
    }
}

impl<'a> ToBoundingBox for DrawableCow<'a> {
    fn bounding_box(&self) -> BoundingBox {
        return self.geometry.bounding_box();
    }
}

impl<'a> From<&'a Drawable> for DrawableCow<'a> {
    fn from(value: &'a Drawable) -> Self {
        return DrawableCow {
            geometry: (&value.geometry).into(),
            style: value.style.clone(),
        };
    }
}

impl<'a> From<&'a DrawableCow<'a>> for DrawableRef<'a> {
    fn from(value: &'a DrawableCow<'a>) -> Self {
        return DrawableRef {
            geometry: (&value.geometry).into(),
            style: value.style.clone(),
        };
    }
}

impl<'a> From<DrawableCow<'a>> for Drawable {
    fn from(value: DrawableCow<'a>) -> Self {
        return Drawable {
            geometry: value.geometry.into(),
            style: value.style,
        };
    }
}

impl<'a, I> From<(I, Style)> for DrawableCow<'a>
where
    I: Into<GeometryCow<'a>>,
{
    fn from(value: (I, Style)) -> Self {
        return DrawableCow {
            geometry: value.0.into(),
            style: value.1,
        };
    }
}
