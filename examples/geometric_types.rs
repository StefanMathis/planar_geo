use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::f64::consts::PI;
use std::path::PathBuf;

pub const UNSCALED_FONT_SIZE: f64 = 16.0;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/{filename}.svg"));
}

fn main() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    // Line segment
    // =========================================================================
    let line_segment =
        LineSegment::new([0.0, 0.0], [2.0, -1.0], e, m).expect("end points not identical");

    let mut style = Style::default();
    style.background_color = Color::from_rgba8(144, 213, 255, 255);
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 2.0;

    let mut bb = BoundingBox::from(&line_segment);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_line_segment");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        return line_segment.draw(&style, cr);
    })
    .expect("image creation failed");

    // Arc segment
    // =========================================================================
    let arc_segment =
        ArcSegment::from_center_radius_start_offset_angle([4.0, 0.0], 1.0, PI, 0.75 * PI, e, m)
            .unwrap();

    let mut bb = BoundingBox::from(&arc_segment);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_arc_segment");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        return arc_segment.draw(&style, cr);
    })
    .expect("image creation failed");

    // Segments
    // =========================================================================

    let line_segment = Segment::LineSegment(line_segment);
    let arc_segment = Segment::ArcSegment(arc_segment);
    let mut bb = BoundingBox::from_bounded_entities([&line_segment, &arc_segment].into_iter())
        .expect("slice has elements");
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_segments");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        line_segment.draw(&style, cr)?;
        arc_segment.draw(&style, cr)?;
        return Ok(());
    })
    .expect("image creation failed");

    // Segment chain
    // =========================================================================
    let mut segment_chain = SegmentChain::new();

    segment_chain.push_back(
        ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, -0.5 * PI, e, m)
            .unwrap()
            .into(),
    );
    segment_chain.extend_back([-1.0, -1.0]);
    segment_chain.extend_back([-1.0, 0.0]);

    let mut bb = BoundingBox::from(&segment_chain);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_segment_chain");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        segment_chain.draw(&style, cr)?;
        return Ok(());
    })
    .expect("image creation failed");

    // Contour
    // =========================================================================
    let contour = Contour::new(segment_chain);

    let fp = file_path("example_contour");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        contour.draw(&style, cr)?;
        return Ok(());
    })
    .expect("image creation failed");

    // Shape
    // =========================================================================
    let hole = Contour::rectangle([-0.9, -0.9], [-0.1, -0.1]);

    let shape = Shape::new(vec![contour, hole]).expect("valid input");

    let fp = file_path("example_shape");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        shape.draw(&style, cr)?;
        return Ok(());
    })
    .expect("image creation failed");
}
