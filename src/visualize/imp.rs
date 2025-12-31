/*!
Implementations of visualization functions for geometric types.

If the `visualize` feature is enabled, all geometric types ([`Segment`],
[`SegmentChain`], [`Contour`] and [`Shape`]) implement a drawing function with
the signature `draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error>`.
This function can be used to draw the type onto a [`cairo::Context`], using the
[`Style`] struct to define the appearance. The vertices / points of the type
are directly used in the coordinate system of the [`cairo::Context`]. For
example, if the [cairo transformation matrix](https://www.cairographics.org/manual/cairo-transformation.html)
matrix does not specify any scaling, a line from (0, 0) to (100, 0) is exactly
100 pixels long when the context is written to a .png file.

To refine the [`Style`] of a type, the following helper structs and enum are
provided by this module:
- [`Style`]: Defines the visual appearance of a type.
- [`LineStyle`]: Defines the stroking pattern for a line, e.g. solid, dashed or
dotted. This enum is used as a field of [`Style`].
- [`Text`]: Defines a text which is placed relative to the geometric type. This
struct is an optional field of [`Style`], but can also be used to display a
"standalone" text via its [`draw`](Text::draw) method (in that case, the text is
placed relative to the origin of the [`Context`](cairo::Context)).
- [`Anchor`]: Defines the orientation of a [`Text`] next to its geometric type
(e.g. on top of it).
*/

use bounding_box::BoundingBox;

use crate::Rotation2;
use crate::composite::Composite;
use crate::contour::Contour;
use crate::segment::{ArcSegment, LineSegment, Segment};
use crate::segment_chain::SegmentChain;
use crate::shape::Shape;

impl Segment {
    /**
    Draws the [`Segment`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        match self {
            Segment::LineSegment(s) => s.draw(style, context),
            Segment::ArcSegment(s) => s.draw(style, context),
        }
    }
}

impl ArcSegment {
    /**
    Draws the [`ArcSegment`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        let vertex = self.start();
        context.move_to(vertex[0], vertex[1]);
        let radius = self.radius();
        let center = self.center();
        let start_angle = self.start_angle();
        let stop_angle = self.stop_angle();

        // Counter-clockwise arc
        if self.is_positive() {
            context.arc(center[0], center[1], radius, start_angle, stop_angle);
        } else {
            context.arc_negative(center[0], center[1], radius, start_angle, stop_angle);
        }

        // Code shared between line and arc segments
        draw_common(style, context)?;

        if let Some(txt) = &style.text {
            let bb = if txt.anchor == crate::visualize::Anchor::Centroid {
                let c = self.centroid();
                BoundingBox::new(c[0], c[0], c[1], c[1])
            } else {
                BoundingBox::from(self)
            };
            txt.draw_with_bounding_box::<true>(context, &bb)?;
        }
        return Ok(());
    }
}

impl LineSegment {
    /**
    Draws the [`LineSegment`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        context.move_to(self.start()[0], self.start()[1]);
        context.line_to(self.stop()[0], self.stop()[1]);

        // Code shared between line and arc segments
        draw_common(style, context)?;

        if let Some(txt) = &style.text {
            let bb = if txt.anchor == crate::visualize::Anchor::Centroid {
                let c = self.centroid();
                BoundingBox::new(c[0], c[0], c[1], c[1])
            } else {
                BoundingBox::from(self)
            };
            txt.draw_with_bounding_box::<true>(context, &bb)?;
        }
        return Ok(());
    }
}

fn draw_common(style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
    // Set the style
    match &style.line_style {
        LineStyle::None => {}
        LineStyle::Dotted => {
            context.set_dash([1.0, 3.0].as_slice(), 0.0);
        }
        LineStyle::Dashed { pattern, offset } => {
            context.set_dash(pattern.as_slice(), *offset);
        }
        LineStyle::Solid => {
            context.set_dash([].as_slice(), 0.0);
        }
    };

    if let LineStyle::None = style.line_style {
        // Do nothing
    } else {
        let lc = &style.line_color;
        context.set_line_width(style.line_width);
        context.set_source_rgba(lc.r.into(), lc.g.into(), lc.b.into(), lc.a.into());
        context.set_line_cap(style.line_cap);

        // Reset the transformation to use the original line width,
        // then stroke, then use the transformation again.
        context.save()?;
        context.identity_matrix();
        context.stroke()?;
        context.restore()?;
    }

    return Ok(());
}

impl SegmentChain {
    /**
    Draws the [`SegmentChain`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        return self.draw_stroking::<false>(style, context);
    }

    fn draw_stroking<const CLOSE: bool>(
        &self,
        style: &Style,
        context: &cairo::Context,
    ) -> Result<(), cairo::Error> {
        self.draw_without_stroking::<CLOSE>(&style, context)?;
        if let LineStyle::None = style.line_style {
            // Do nothing
        } else {
            let lc = &style.line_color;
            context.set_line_width(style.line_width);
            context.set_source_rgba(lc.r.into(), lc.g.into(), lc.b.into(), lc.a.into());
            context.set_line_cap(style.line_cap);
            context.set_line_join(style.line_join);

            // Reset the transformation to use the original line width,
            // then stroke, then use the transformation again.
            context.save()?;
            context.identity_matrix();
            context.stroke()?;
            context.restore()?;
        }

        if let Some(txt) = &style.text {
            let bb = if txt.anchor == crate::visualize::Anchor::Centroid {
                let c = self.centroid();
                BoundingBox::new(c[0], c[0], c[1], c[1])
            } else {
                BoundingBox::from(self)
            };
            txt.draw_with_bounding_box::<true>(context, &bb)?;
        }
        return Ok(());
    }

    fn draw_without_stroking<const CLOSE: bool>(
        &self,
        style: &Style,
        context: &cairo::Context,
    ) -> Result<(), cairo::Error> {
        // Initialize the Cairo path
        if let Some(first_segment) = self.front() {
            match &style.line_style {
                LineStyle::None => {}
                LineStyle::Dotted => {
                    context.set_dash([1.0, 3.0].as_slice(), 0.0);
                }
                LineStyle::Dashed { pattern, offset } => {
                    context.set_dash(pattern.as_slice(), *offset);
                }
                LineStyle::Solid => {
                    context.set_dash([].as_slice(), 0.0);
                }
            };

            let vertex = first_segment.start();

            context.move_to(vertex[0], vertex[1]);

            for segment in self.iter() {
                match segment {
                    Segment::LineSegment(ls) => context.line_to(ls.stop()[0], ls.stop()[1]),
                    Segment::ArcSegment(arc) => {
                        let radius = arc.radius();
                        let center = arc.center();
                        let start_angle = arc.start_angle();
                        let stop_angle = arc.stop_angle();

                        // Counter-clockwise arc
                        if arc.is_positive() {
                            context.arc(center[0], center[1], radius, start_angle, stop_angle);
                        } else {
                            context.arc_negative(
                                center[0],
                                center[1],
                                radius,
                                start_angle,
                                stop_angle,
                            );
                        }
                    }
                }
            }
            if CLOSE {
                context.close_path();
            }
        }
        return Ok(());
    }
}

impl Contour {
    /**
    Draws the [`Contour`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        return self.segment_chain().draw_stroking::<true>(style, context);
    }
}

impl Shape {
    /**
    Draws the [`Shape`] onto the [`cairo::Context`] with the given [`Style`].
    See the [module level documentation](crate::visualize).
     */
    pub fn draw(&self, style: &Style, context: &cairo::Context) -> Result<(), cairo::Error> {
        let draw_contour = if let LineStyle::None = style.line_style {
            false
        } else {
            true
        };
        let fill_contour = style.background_color.a != 0.0;

        // Skip the entire evaluation in case no border and no background color has been given
        if draw_contour || fill_contour {
            self.contour()
                .segment_chain()
                .draw_without_stroking::<true>(style, context)?;
            // Fill path while substracting inner paths
            for hole in self.holes() {
                context.new_sub_path();
                hole.segment_chain()
                    .draw_without_stroking::<true>(style, context)?;
            }

            // Substract the holes regardless of the direction of the segment_chain (clockwise or counter-clockwise)
            if fill_contour {
                context.set_fill_rule(cairo::FillRule::EvenOdd);
                let fc = &style.background_color;
                context.set_source_rgba(fc.r.into(), fc.g.into(), fc.b.into(), fc.a.into());
                context.fill_preserve()?;
            }

            if draw_contour {
                let lc = &style.line_color;
                context.set_line_width(style.line_width);
                context.set_source_rgba(lc.r.into(), lc.g.into(), lc.b.into(), lc.a.into());
                context.set_line_cap(style.line_cap);
                context.set_line_join(style.line_join);

                // Reset the transformation to use the original line width,
                // then stroke, then use the transformation again.
                context.save()?;
                context.identity_matrix();
                context.stroke()?;
                context.restore()?;
            }
        }

        if let Some(txt) = &style.text {
            let bb = if txt.anchor == crate::visualize::Anchor::Centroid {
                let c = self.centroid();
                BoundingBox::new(c[0], c[0], c[1], c[1])
            } else {
                BoundingBox::from(self)
            };
            txt.draw_with_bounding_box::<true>(context, &bb)?;
        }
        return Ok(());
    }
}

/**
A representation of a color in RGB and A(lpha) values.

When used to draw an object, all values are clamped to [0, 1] by
[cairo](https://www.cairographics.org/manual/cairo-cairo-t.html#cairo-set-source-rgba).

For the color components, 0 means no component and 1 means maximum component.
If e.g. applied to the RGB8 color space, 0 corresponds to 0 and 1 corresponds
to 255.

The alpha value defines the translucency of the color. 0 means that the color
is invisible, while 1 means the color is opaque.
 */
#[derive(Copy, Clone, Debug)]
pub struct Color {
    /// Red component.
    pub r: f32,
    /// Green component.
    pub g: f32,
    /// Blue component.
    pub b: f32,
    /// Alpha value of the color.
    pub a: f32,
}

impl Color {
    /**
    Creates a new instance of [`Color`] from its components.
     */
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        return Self { r, g, b, a };
    }

    /**
    Creates a new instance of [`Color`] from RGB8 values (integers from 0 to 255).
     */
    pub fn from_rgba8(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self {
            r: r as f32 / 255.0,
            g: g as f32 / 255.0,
            b: b as f32 / 255.0,
            a: a as f32 / 255.0,
        }
    }
}

/**
This struct defines the visual appearance of a geometric object - i.e. a
[`SegmentChain`], a [`Contour`] or a [`Shape`]. These properties are derived
from the CSS standard, see links in field descriptions.
*/
#[derive(Debug, Clone)]
pub struct Style {
    /// Color of all lines of the geometric object. Corresponds to "border-color"
    /// as defined in <https://www.w3.org/TR/css-backgrounds-3/#border-color>.
    pub line_color: Color,
    /// Fill color for closed surface, e.g. as they appear in a [`Shape`].
    /// For other geometric objects, this property is meaningless.
    /// Corresponds to "background-color" as defined in
    /// <https://www.w3.org/TR/css-backgrounds-3/#background-color>.
    pub background_color: Color,
    /// Width of all lines of the geometric object in . Corresponds to "border-width"
    /// as defined in <https://www.w3.org/TR/css-backgrounds-3/#the-border-width>.
    pub line_width: f64,
    /// Style of all lines of the geometric object. See docstring of
    /// [`LineStyle`] for an example. Corresponds to "border-style" as defined
    /// in <https://www.w3.org/TR/css-backgrounds-3/#the-border-style>.
    pub line_style: LineStyle,
    /// How the end points of lines are rendered. See
    /// <https://www.cairographics.org/manual/cairo-cairo-t.html#cairo-line-cap-t>.
    pub line_cap: cairo::LineCap,
    /// How the junction points of individual [`Segment`]s are rendered. See
    /// <https://www.cairographics.org/manual/cairo-cairo-t.html#cairo-line-cap-t>.
    pub line_join: cairo::LineJoin,
    /// Specifies a [`Text`] which is displayed next to the geometric object.
    /// The location is defined by [`Text::anchor`], [`Text::fixed_anchor_offset`]
    /// and [`Text::scaled_anchor_offset`]. See the docstring of [`Text`].
    /// The [`Text`] is boxed to minimize the size of [`Style`] in the common
    /// case that no text is used in the visual representation of a geometric
    /// type.
    pub text: Option<Box<Text>>,
}

impl Style {
    /**
    Creates a new instance of [`Style`] from its components.
     */
    pub fn new(
        line_color: Color,
        background_color: Color,
        line_width: f64,
        line_style: LineStyle,
        line_cap: cairo::LineCap,
        line_join: cairo::LineJoin,
        text: Option<Box<Text>>,
    ) -> Self {
        return Style {
            line_color,
            background_color,
            line_width,
            line_style,
            line_cap,
            line_join,
            text,
        };
    }
}

impl Default for Style {
    fn default() -> Self {
        let black = Color::new(0.0, 0.0, 0.0, 1.0);
        let white = Color::new(1.0, 1.0, 1.0, 1.0);
        return Style {
            line_color: black,
            background_color: white,
            line_width: 0.5,
            line_style: LineStyle::Solid,
            line_cap: cairo::LineCap::Round,
            line_join: cairo::LineJoin::Miter,
            text: None,
        };
    }
}

/**
Definition of the line style.
*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/line_styles.svg \"Line style comparison\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/line_styles.svg)"
)]
/**
```rust
use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use cairo;

fn create_fn(cr: &cairo::Context) -> Result<(), cairo::Error> {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    // Set the background to white
    cr.set_source_rgb(1.0, 1.0, 1.0);
    cr.paint()?;

    let mut txt = Text {
        text: "".into(),
        anchor: Anchor::Bottom,
        fixed_anchor_offset: [0.0, -5.0],
        scaled_anchor_offset: [0.0, 0.0],
        color: Color::new(0.0, 0.0, 0.0, 1.0),
        font_size: 16.0,
        angle: 0.0,
    };

    let mut style = Style::default();
    style.line_width = 2.0;
    style.line_cap = cairo::LineCap::Butt;

    let segment = LineSegment::new([0.1, 0.5], [0.7, 0.5], e, m).unwrap();
    txt.text = "Solid".into();
    style.line_style = LineStyle::Solid;
    style.text = Some(Box::new(txt.clone()));
    segment.draw(&style, cr)?;

    let segment = LineSegment::new([0.8, 0.5], [1.4, 0.5], e, m).unwrap();
    txt.text = "Dotted".into();
    style.line_style = LineStyle::Dotted;
    style.text = Some(Box::new(txt.clone()));
    segment.draw(&style, cr)?;

    let segment = LineSegment::new([1.5, 0.5], [2.1, 0.5], e, m).unwrap();
    txt.text = "Dashed".into();
    style.line_style = LineStyle::default_dashed();
    style.text = Some(Box::new(txt.clone()));
    segment.draw(&style, cr)?;

    let segment = LineSegment::new([2.2, 0.5], [2.8, 0.5], e, m).unwrap();
    txt.text = "None".into();
    style.line_style = LineStyle::None;
    style.text = Some(Box::new(txt.clone()));
    segment.draw(&style, cr)?;

    return Ok(());
}

let view = Viewport::from_bounding_box(
    &BoundingBox::new(0.0, 2.9, 0.35, 0.55),
    SideLength::Long(600),
);

// Comment this in to create the image
// view.write_to_file("docs/line_styles.svg", create_fn).unwrap();
#
# assert!(
#     view.compare_or_create(std::path::Path::new("tests/img/line_styles.png"), create_fn, 0.98)
#         .is_ok()
# );
```
 */
#[derive(Debug, Clone)]
pub enum LineStyle {
    /// No line is drawn
    None,
    /// Draws a dotted border. Equal to [`LineStyle::Dashed`] with a pattern of
    /// `(1.0, 3.0)` and an offset of 0.
    Dotted,
    /// Draws a dashed border with the defined pattern. See
    /// <https://www.cairographics.org/manual/cairo-cairo-t.html#cairo-set-dash>.
    Dashed {
        /// Each value of the vector provides the length of alternate "on" and
        /// "off" portions of the dash pattern stroke. If any value is negative
        /// or if all values are 0, the underlying [cairo] library will return
        /// an error.
        pattern: Vec<f64>,
        /// An "off" portion offset which specifies an offset into the `pattern`
        /// at which the stroke begins.
        offset: f64,
    },
    /// Draws a solid line.
    Solid,
}

impl LineStyle {
    /**
    Returns [`LineStyle::Dashed{pattern: vec![4.0, 4.0], offset: 0.0}`](`LineStyle::Dashed`).
    This is a default dash pattern which is symmetric (non-solid and solid
    parts are equal), see [`LineStyle`] docstring.
     */
    pub fn default_dashed() -> Self {
        return LineStyle::Dashed {
            pattern: vec![4.0, 4.0],
            offset: 0.0,
        };
    }
}

/**
A displayable string. Besides the string itself, it also has a couple of
properties which define the placement and depiction of the string on the canvas.

Two properties which warrant further explanation are `fixed_anchor_offset` and
`scaled_anchor_offset`. The two images below show the same two
[`Contour`](crate::prelude::Contour) objects, but the
[`Viewport`](cairo_viewport::Viewport) of the second image has been scaled by
2. On the left side, `scaled_anchor_offset` is set to `(0.3, 0.3)` and on the
right side, `fixed_anchor_offset` is set to `(10, 10)`. Comparing the two images
show that the offset of the text from the top left corner is scaled with the
[`Viewport`](cairo_viewport::Viewport) for the left
[`Contour`](crate::prelude::Contour) , but not for the right
[`Contour`](crate::prelude::Contour). The font size is also not scaled.

If a [`Text`] is used as part of a [`Style`] for a geometric object (such as a
[`Shape`] or [`LineSegment`]), the offsets are applied to the [`Anchor`] of the
object (for example, in the images they are applied to the top left corner of
the contour). If a [`Text`] is drawn on its own, the offset is applied to the
origin of the [`cairo::Context`] coordinates. See docstring of [`Anchor`] for
an example.
*/
#[cfg_attr(
    docsrs,
    doc = "\n\nNormal image:\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/anchor_offset_scale_1.svg)"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\nNormal image:\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/anchor_offset_scale_1.svg)"
)]
#[cfg_attr(
    docsrs,
    doc = "\n\nZoomed in view:\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/anchor_offset_scale_2.svg)"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\nZoomed in view:\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/anchor_offset_scale_2.svg)"
)]
/**

These images were created with the code below:

```rust
use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;

fn draw_fn() -> Box<dyn FnOnce(&cairo::Context) -> Result<(), cairo::Error>> {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let mut style = Style::default();
    style.line_color = Color::new(1.0, 0.5, 0.5, 1.0);
    style.line_width = 2.0;

    let mut contour = Contour::new(SegmentChain::from_points(&[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]]));

    let drawing_fn = move |cr: &cairo::Context| {
        // Set the background to white
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        // Left box
        let txt = Text {
            text: "Scaled offset".into(),
            anchor: Anchor::TopLeft,
            fixed_anchor_offset: [0.0, 0.0],
            scaled_anchor_offset: [0.3, 0.3], // Values in contour scale
            color: Color::new(0.0, 0.0, 0.0, 1.0),
            font_size: 16.0,
            angle: 0.0,
        };
        style.text = Some(Box::new(txt));
        contour.draw(&style, cr)?;

        // Right box
        let txt = Text {
            text: "Fixed offset".into(),
            anchor: Anchor::TopLeft,
            fixed_anchor_offset: [10.0, 10.0], // Values in image scale (pixels)
            scaled_anchor_offset: [0.0, 0.0],
            color: Color::new(0.0, 0.0, 0.0, 1.0),
            font_size: 16.0,
            angle: 0.0,
        };
        style.line_color = Color::new(0.5, 0.8, 0.5, 1.0);
        style.text = Some(Box::new(txt));
        contour.translate([1.0, 0.0]);
        contour.draw(&style, cr)?;

        return Ok(());
    };
    return Box::new(drawing_fn);
}

// Prepare the view for the first image
let mut view =
    Viewport::from_bounding_box(&BoundingBox::new(0.0, 3.0, 0.0, 1.0), SideLength::Long(600));

// Comment this in to create the first image
// view.write_to_file("docs/anchor_offset_scale_1.svg", draw_fn()).unwrap();
#
# assert!(
#     view.compare_or_create(
#         std::path::Path::new("tests/img/anchor_offset_scale_1.png"),
#         draw_fn(),
#         0.98
#     )
#     .is_ok()
# );

// Zoom into the image
view.scale *= 2.0;

// Comment this in to create the second image
// view.write_to_file("docs/anchor_offset_scale_2.svg", draw_fn()).unwrap();
#
# assert!(
#    view.compare_or_create(
#         std::path::Path::new("tests/img/anchor_offset_scale_2.png"),
#         draw_fn(),
#         0.98
#     )
#    .is_ok()
# );
```
 */
#[derive(Debug, Clone)]
pub struct Text {
    /// A String containing the displayed text.
    pub text: String,
    /// Specifies the text position relative to a reference bounding box /
    /// point. See the docstring of [`Anchor`] for more.
    pub anchor: Anchor,
    /// Defines an offset from the reference bounding box / origin which
    /// is not scaled with `context.matrix().xx()`. See [`Text`] docstring.
    pub fixed_anchor_offset: [f64; 2],
    /// Defines an offset from the reference bounding box / origin which
    /// is scaled with `context.matrix().xx()`. See [`Text`] docstring.
    pub scaled_anchor_offset: [f64; 2],
    /// Defines the text color.
    pub color: Color,
    /// Defines the text font size in points (pt). If set to a negative value,
    /// the text is not shown. A pt is roughly 0.3528 mm on a screen. The
    /// font size is not scaled with the canvas, see example images in the
    /// [`Text`] docstring.
    pub font_size: f64,
    /// Rotational angle of the text in rad.
    pub angle: f64,
}

impl Text {
    /**
    Creates a new instance of [`Text`] from its components.
     */
    pub fn new(
        text: String,
        anchor: Anchor,
        fixed_anchor_offset: [f64; 2],
        scaled_anchor_offset: [f64; 2],
        color: Color,
        font_size: f64,
        angle: f64,
    ) -> Self {
        return Text {
            text,
            anchor,
            fixed_anchor_offset,
            scaled_anchor_offset,
            color,
            font_size,
            angle,
        };
    }

    /**
    Draw the [`Text`] onto the given [`cairo::Context`].

    The placement of the text is defined by the following fields:
    [`Text::anchor`], [`Text::fixed_anchor_offset`] and
    [`Text::scaled_anchor_offset`]. See the docstrings of [`Text`] and
    [`Anchor`] for details and examples.
     */
    pub fn draw(&self, context: &cairo::Context) -> Result<(), cairo::Error> {
        // This bounding box is a point. Together with INSIDE_BB = false, this
        // results in the text being placed around the point as defined by
        // the anchor.
        let bb = BoundingBox::new(0.0, 0.0, 0.0, 0.0);
        return self.draw_with_bounding_box::<false>(context, &bb);
    }

    /**
    This function performs the actual text drawing. It can be used either to
    draw a text "standalone" (INSIDE_BB = false) or within a bounding box of a
    geometric object (INSIDE_BB = true). In the latter mode, it is used
    throughout the drawing functions of different objects.
     */
    pub(crate) fn draw_with_bounding_box<const INSIDE_BB: bool>(
        &self,
        context: &cairo::Context,
        bounding_box: &BoundingBox,
    ) -> Result<(), cairo::Error> {
        let bb = bounding_box;
        let xmean = (bb.xmin() + bb.xmax()) / 2.0;
        let ymean = (bb.ymin() + bb.ymax()) / 2.0;
        context.set_source_rgba(
            self.color.r.into(),
            self.color.g.into(),
            self.color.b.into(),
            self.color.a.into(),
        );
        let font_size = self.font_size / context.matrix().xx(); // Is equivalent to context.matrix().yy
        context.set_font_size(font_size);

        let extents = context.text_extents(&self.text)?;
        let x: f64;
        let y: f64;

        let anchor = if INSIDE_BB {
            self.anchor
        } else {
            self.anchor.opposite()
        };

        match anchor {
            Anchor::Center | Anchor::Centroid => {
                x = xmean - (extents.width() * 0.5 + extents.x_bearing());
                y = ymean - (extents.height() * 0.5 + extents.y_bearing());
            }
            Anchor::Top => {
                x = xmean - (extents.width() * 0.5 + extents.x_bearing());
                y = bb.ymin() - extents.y_bearing();
            }
            Anchor::TopRight => {
                x = bb.xmax() - (extents.width() + extents.x_bearing());
                y = bb.ymin() - extents.y_bearing();
            }
            Anchor::Right => {
                x = bb.xmax() - (extents.width() + extents.x_bearing());
                y = ymean - (extents.height() * 0.5 + extents.y_bearing());
            }
            Anchor::BottomRight => {
                x = bb.xmax() - (extents.width() + extents.x_bearing());
                y = bb.ymax();
            }
            Anchor::Bottom => {
                x = xmean - (extents.width() * 0.5 + extents.x_bearing());
                y = bb.ymax();
            }
            Anchor::BottomLeft => {
                x = bb.xmin();
                y = bb.ymax();
            }
            Anchor::Left => {
                x = bb.xmin();
                y = ymean - (extents.height() * 0.5 + extents.y_bearing());
            }
            Anchor::TopLeft => {
                x = bb.xmin();
                y = bb.ymin() - extents.y_bearing();
            }
        }

        // Offset to anchor position y/y. Since the user coordinates are already
        // transformed, the fixed offset is transformed back into screen coordinates
        let offset_x =
            self.fixed_anchor_offset[0] / context.matrix().xx() + self.scaled_anchor_offset[0];
        let offset_y =
            self.fixed_anchor_offset[1] / context.matrix().yy() + self.scaled_anchor_offset[1];

        // Short-cut if the text is not rotated
        if self.angle == 0.0 {
            context.new_path();
            context.move_to(x + offset_x, y + offset_y);
            context.show_text(&self.text)?;
        } else {
            /*
            The text needs to be rotated around its center. To rotate a 2d entity (e.g. a point P)
            around a rotation center C, the entity first needs to be translated into a new
            coordinate system where the rotation center equals the origin.
            Then, the points of the entity are rotated by an rotational matrix R and then
            the entity is translated back into its original coordinate system (P').

            P' = R * (P - C) + C

            Since the "rotate" command of Cairo rotates any subsequent commands around the upper left
            corner, the center of the text box needs to be "counter-rotated" before the
            text can be shown.
            */
            // Get the context origin (distance to the upper left corner)
            let x0 = context.matrix().x0() / context.matrix().xx();
            let y0 = context.matrix().y0() / context.matrix().yy();

            // Lower left corner of the text box relative to the upper left context corner
            // (in Cairo coordinates, e.g. with the y-axis pointing down)
            let text_box_ll = [x + offset_x + x0, y + offset_y + y0];

            // Center of the text box
            let center_text_box = [
                text_box_ll[0] + 0.5 * extents.width(),
                text_box_ll[1] - 0.5 * extents.height(),
            ];

            let delta_ll_to_center = [-0.5 * extents.width(), 0.5 * extents.height()];
            let center_text_box_rot = Rotation2::new(-self.angle) * center_text_box; // Counter rotation

            context.save()?;
            context.translate(-x0, -y0); // Reset the Cairo context origin
            context.rotate(self.angle); // Rotate the text

            /*
            center_text_box_rot is rotated by Cairo => compensated by counter rotation,
            lower left corner of the text box appears unchanged in user coordinates.
            Adding delta_ll_to_center => Center of the text box appears unchanged in user coordinates.
            */
            context.move_to(
                center_text_box_rot[0] + delta_ll_to_center[0],
                center_text_box_rot[1] + delta_ll_to_center[1],
            );
            context.show_text(&self.text)?;
            context.restore()?;
        }
        return Ok(());
    }
}

/**
An [`Anchor`] describes the position of a [`Text`] relative to a reference point.

When a [`Text`] object is used as part of a [`Style`] struct to define the
visualization of a geometric object, the [`Anchor`] determines how the text is
placed relative to the bounding box of the geometric object. However, a [`Text`]
can also be shown "standalone" using [`Text::draw`]. In that case, the
[`Anchor`] defines the text placement relative to a point defined by
[`Text::fixed_anchor_offset`] and [`Text::scaled_anchor_offset`]. See docstring
of [`Text`].

The image below shows both cases: On the left side, the red bounding box of an
geometric object and the associated text placement options are shown. On the
right side, the [`Text`] has been drawn "standalone" (the text placement point
is visualized as a red cross).
*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/text_placement.svg \"Text placement with and without associated geometric object\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/text_placement.svg)"
)]
/**

This image was created with the following code:
```
use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;

let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

let mut style = Style::default();
style.line_color = Color::new(1.0, 0.5, 0.5, 1.0);
style.line_width = 2.0;
let mut txt = Text {
    text: "Placeholder".into(),
    anchor: Anchor::Center,
    fixed_anchor_offset: [0.0, 0.0],
    scaled_anchor_offset: [0.0, 0.0],
    color: Color::new(0.0, 0.0, 0.0, 1.0),
    font_size: 16.0,
    angle: 0.0,
};

// This contour is equivalent to its bounding box, since it is just a rectangle.
let contour = Contour::new(SegmentChain::from_points(&[[0.1, 0.1], [1.9, 0.1], [1.9, 0.9], [0.1, 0.9]]));

// A cross on the right of the contour which shows placement of standalone texts
let cc = [2.5, 0.5];
let hori = LineSegment::new([cc[0] - 0.2, cc[1]], [cc[0] + 0.2, cc[1]], e, m)
    .expect("points not identical");
let vert = LineSegment::new([cc[0], cc[1] - 0.2], [cc[0], cc[1] + 0.2], e, m)
    .expect("points not identical");

let view =
    Viewport::from_bounding_box(&BoundingBox::new(0.0, 3.0, 0.0, 1.0), SideLength::Long(600));

// This function can be used to draw the image
let draw_fn = |cr: &cairo::Context| {
    // Set the background to white
    cr.set_source_rgb(1.0, 1.0, 1.0);
    cr.paint()?;

    // Draw the contour first
    contour.draw(&style, cr)?;

    // Draw the cross
    hori.draw(&style, cr)?;
    vert.draw(&style, cr)?;

    // Then draw the contour again with an invisible line color to demonstrate the text placement
    style.line_color = Color::new(0.0, 0.0, 0.0, 0.0);

    let anchors_str = &[
        "Center",
        "TopLeft",
        "Top",
        "TopRight",
        "Right",
        "BottomRight",
        "Bottom",
        "BottomLeft",
        "Left",
    ];
    let anchors = &[
        Anchor::Center,
        Anchor::TopLeft,
        Anchor::Top,
        Anchor::TopRight,
        Anchor::Right,
        Anchor::BottomRight,
        Anchor::Bottom,
        Anchor::BottomLeft,
        Anchor::Left,
    ];
    for (anchor_str, anchor) in anchors_str.iter().zip(anchors.iter()) {
        txt.anchor = *anchor;
        txt.text = anchor_str.to_string();
        style.text = Some(Box::new(txt.clone()));
        contour.draw(&style, cr)?;
    }

    // Show text placement without an associated geometric type
    txt.scaled_anchor_offset = cc;
    let anchors_str = &["TopLeft", "TopRight", "BottomRight", "BottomLeft"];
    let anchors = &[
        Anchor::TopLeft,
        Anchor::TopRight,
        Anchor::BottomRight,
        Anchor::BottomLeft,
    ];
    for (anchor_str, anchor) in anchors_str.iter().zip(anchors.iter()) {
        txt.anchor = *anchor;
        txt.text = anchor_str.to_string();
        txt.draw(cr)?;
    }

    return Ok(());
};

// Comment this in to actually create the shown image
// view.write_to_file("docs/text_placement.svg", draw_fn).expect("image could not be created");
# assert!(
#     view.compare_or_create(
#         std::path::Path::new("tests/img/text_placement.png"),
#         draw_fn,
#         0.98
#     )
#     .is_ok()
# );
```
*/
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    /// If used to style a geometric type, the text is placed at its centroid
    /// (center of mass). If used to visualize a [`Text`] object directly, this
    /// variant behaves identical to [`Anchor::Center`].
    Centroid,
    /// Center of text extents is placed at the bounding box center / at the
    /// reference point.
    Center,
    /// Top left corner of text extents is placed in the top left corner of the
    /// bounding box / reference point.
    TopLeft,
    /// Top middle of the text extents is placed in the middle of the top edge
    /// of the bounding box / at the top of the reference point.
    Top,
    /// Top right corner of text extents is placed in the top right corner of the
    /// bounding box / reference point.
    TopRight,
    /// Right middle of the text extents is placed in the middle of the right edge
    /// of the bounding box / right to the reference point.
    Right,
    /// Bottom right corner of text extents is placed in the bottom right corner of the
    /// bounding box / reference point.
    BottomRight,
    /// Top middle of the text extents is placed in the middle of the bottom edge
    /// of the bounding box / at the bottom of the reference point.
    Bottom,
    /// Bottom left corner of text extents is placed in the bottom left corner of the
    /// bounding box / reference point.
    BottomLeft,
    /// Left middle of the text extents is placed in the middle of the left edge
    /// of the bounding box / left to the reference point.
    Left,
}

impl Anchor {
    /**
    Returns the "opposite" [`Anchor`]. See the examples below.

    ```
    use planar_geo::visualize::Anchor;

    // Centroid and center are not changed
    assert_eq!(Anchor::Centroid.opposite(), Anchor::Centroid);
    assert_eq!(Anchor::Center.opposite(), Anchor::Center);

    // All directional anchors are inverted
    assert_eq!(Anchor::TopLeft.opposite(), Anchor::BottomRight);
    assert_eq!(Anchor::Top.opposite(), Anchor::Bottom);
    assert_eq!(Anchor::TopRight.opposite(), Anchor::BottomLeft);
    assert_eq!(Anchor::Right.opposite(), Anchor::Left);
    assert_eq!(Anchor::BottomRight.opposite(), Anchor::TopLeft);
    assert_eq!(Anchor::Bottom.opposite(), Anchor::Top);
    assert_eq!(Anchor::BottomLeft.opposite(), Anchor::TopRight);
    assert_eq!(Anchor::Left.opposite(), Anchor::Right);
    ```
     */
    pub fn opposite(&self) -> Self {
        match self {
            Anchor::Centroid => Anchor::Centroid,
            Anchor::Center => Anchor::Center,
            Anchor::TopLeft => Anchor::BottomRight,
            Anchor::Top => Anchor::Bottom,
            Anchor::TopRight => Anchor::BottomLeft,
            Anchor::Right => Anchor::Left,
            Anchor::BottomRight => Anchor::TopLeft,
            Anchor::Bottom => Anchor::Top,
            Anchor::BottomLeft => Anchor::TopRight,
            Anchor::Left => Anchor::Right,
        }
    }
}
