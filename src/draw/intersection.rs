/*!
Module containing the [`Intersection::draw`] implementation. See its docstring
for more.
 */

use crate::composite::{Intersection, SegmentKey};
use crate::draw::{Color, DrawableRef, Style};
use crate::geometry::GeometryRef;
use crate::prelude::Polysegment;
use crate::shape::Shape;

/**
[`Intersection`]s are drawn as crosses, whose properties are defined by
this struct. See the [module documentation](crate::draw::intersection) for
an example.
 */
pub struct IntersectionStyle {
    /// Color of all lines of the geometric object. Corresponds to
    /// "border-color" as defined in <https://www.w3.org/TR/css-backgrounds-3/#border-color>.
    pub line_color: Color,
    /// Width of all lines of the geometric object. Corresponds to
    /// "border-color" as defined in <https://www.w3.org/TR/css-backgrounds-3/#the-border-width>.
    pub line_width: f64,
    /// How the end points of lines are rendered. See
    /// <https://www.cairographics.org/manual/cairo-cairo-t.html#cairo-line-cap-t>.
    pub line_cap: cairo::LineCap,
    /// Height / width of the intersection cross in points.
    pub cross_size: f64,
}

impl Default for IntersectionStyle {
    fn default() -> Self {
        Self {
            line_color: Color::new(1.0, 0.0, 0.0, 1.0),
            line_width: 5.0,
            line_cap: cairo::LineCap::Butt,
            cross_size: 10.0,
        }
    }
}

fn draw_segment_of_polysegment(
    polysegment: &Polysegment,
    style: &Style,
    segment_key: SegmentKey,
    context: &cairo::Context,
) -> Result<(), cairo::Error> {
    if let Some(segment) = polysegment.get(segment_key.segment_idx) {
        segment.draw(style, context)?;
    }
    return Ok(());
}

fn draw_segment_of_shape(
    shape: &Shape,
    style: &Style,
    segment_key: SegmentKey,
    context: &cairo::Context,
) -> Result<(), cairo::Error> {
    if let Some(contour) = shape.contours().get(segment_key.contour_idx) {
        draw_segment_of_polysegment(contour.polysegment(), style, segment_key, context)?;
    }
    return Ok(());
}

fn draw_segment_of_drawable(
    drawable: DrawableRef<'_>,
    segment_key: SegmentKey,
    context: &cairo::Context,
) -> Result<(), cairo::Error> {
    match drawable.geometry {
        GeometryRef::Point(_) => Ok(()),
        GeometryRef::BoundingBox(_) => Ok(()),
        GeometryRef::ArcSegment(arc_segment) => arc_segment.draw(&drawable.style, context),
        GeometryRef::LineSegment(line_segment) => line_segment.draw(&drawable.style, context),
        GeometryRef::Line(_) => Ok(()),
        GeometryRef::Segment(segment) => segment.draw(&drawable.style, context),
        GeometryRef::Polysegment(polysegment) => {
            draw_segment_of_polysegment(polysegment, &drawable.style, segment_key, context)
        }
        GeometryRef::Contour(contour) => draw_segment_of_polysegment(
            contour.polysegment(),
            &drawable.style,
            segment_key,
            context,
        ),
        GeometryRef::Shape(shape) => {
            draw_segment_of_shape(shape, &drawable.style, segment_key, context)
        }
    }
}

impl Intersection {
    /// Draws the intersection between two [`DrawableRef`]s onto the given
    /// [`cairo::Context`] using the given [`IntersectionStyle`].
    ///
    /// Since an intersection is just a point, it is represented as a "X"
    /// when drawing it. The visual properties of that "+" are provided by
    /// the given [`IntersectionStyle`]. If `left` or `right` are provided
    /// and are [`Composite`](crate::composite::Composite)s, the
    /// [`Intersection::left`] or [`Intersection::right`] key are used to
    /// access the respective segment. The segment is then drawn using
    /// [`DrawableRef::style`].
    ///
    /// The image below shows a [`Shape`] and an intersecting [`Polysegment`].
    /// The intersecting segments have been colored yellow.
    #[doc = ""]
    #[cfg_attr(feature = "doc-images", doc = "![Intersection][intersection_example]")]
    #[cfg_attr(
        feature = "doc-images",
        embed_doc_image::embed_doc_image(
            "intersection_example",
            "docs/img/intersection_example.svg"
        )
    )]
    #[cfg_attr(
        not(feature = "doc-images"),
        doc = "**Doc images not enabled**. Compile docs with `cargo doc --features 'doc-images'` and Rust version >= 1.54."
    )]
    ///
    /// This image from the [lib](crate) module documentation was created using
    /// the following code:
    /// ```
    /// use cairo_viewport::{SideLength, Viewport};
    /// use planar_geo::prelude::*;
    ///
    /// let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    /// let contour = Contour::new(Polysegment::from_points(vertices));
    ///
    /// let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    /// let hole = Contour::new(Polysegment::from_points(vertices));
    ///
    /// let shape = Shape::new(vec![contour.clone(), hole]).expect("valid inputs");
    ///
    /// let polysegment = Polysegment::from_points(&[[-1.0, 1.0], [-0.5, 0.5], [1.5, 0.5], [2.0, 1.0]]);
    ///
    /// let view = Viewport::from_bounding_box(
    ///     &BoundingBox::new(-1.2, 2.2, -0.1, 1.1),
    ///     SideLength::Long(500),
    /// );
    ///
    /// let mut style = Style::default();
    /// style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    /// style.line_width = 2.0;
    /// style.line_style = LineStyle::Solid;
    /// style.background_color = Color::from_rgba8(144, 213, 255, 255);
    ///
    /// // Yellow line color for the intersecting segments.
    /// let mut intersected_segments_style = Style::default();
    /// intersected_segments_style.line_color = Color::new(1.0, 1.0, 0.0, 1.0);
    /// intersected_segments_style.line_width = 3.0;
    /// intersected_segments_style.line_style = LineStyle::Solid;
    ///
    /// let intersection_style = IntersectionStyle::default();
    ///
    /// let draw_fn = |cr: &cairo::Context| {
    ///     // Set the background to white
    ///     cr.set_source_rgb(1.0, 1.0, 1.0);
    ///     cr.paint()?;
    ///
    ///     shape.draw(&style, cr)?;
    ///     polysegment.draw(&style, cr)?;
    ///     for i in shape.intersections_polysegment(&polysegment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
    ///         i.draw(
    ///             &intersection_style,
    ///             Some(DrawableRef::new(&shape, intersected_segments_style.clone())),
    ///             Some(DrawableRef::new(&polysegment, intersected_segments_style.clone())),
    ///             cr,
    ///         )?;
    ///     }
    ///     return Ok(());
    /// };
    ///
    /// // Comment this in to actually create the shown image
    /// // view.write_to_file("docs/img/intersection_example.svg", draw_fn).expect("image could not be created");
    /// # assert!(
    /// #     view.compare_or_create(
    /// #         std::path::Path::new("tests/img/intersection_polysegment_shape.png"),
    /// #         draw_fn, 0.99
    /// #     )
    /// #     .is_ok()
    /// # );
    /// ```
    pub fn draw(
        &self,
        style: &IntersectionStyle,
        left: Option<DrawableRef<'_>>,
        right: Option<DrawableRef<'_>>,
        context: &cairo::Context,
    ) -> Result<(), cairo::Error> {
        if let Some(s) = left {
            draw_segment_of_drawable(s, self.left, context)?;
        }

        if let Some(s) = right {
            draw_segment_of_drawable(s, self.right, context)?;
        }

        context.set_dash([].as_slice(), 0.0);
        context.set_line_width(style.line_width);
        let lc = &style.line_color;
        context.set_source_rgba(lc.r.into(), lc.g.into(), lc.b.into(), lc.a.into());

        let [x, y] = self.point;

        // First line
        let offset = 0.5 * style.cross_size / context.matrix().xx();
        context.move_to(x - offset, y - offset);
        context.line_to(x + offset, y + offset);
        context.save()?;
        context.identity_matrix();
        context.stroke()?;
        context.restore()?;

        // Second line
        context.move_to(x + offset, y - offset);
        context.line_to(x - offset, y + offset);
        context.save()?;
        context.identity_matrix();
        context.stroke()?;
        context.restore()?;
        return Ok(());
    }
}
