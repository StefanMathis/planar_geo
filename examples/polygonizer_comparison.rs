use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::f64::consts::PI;
use std::path::PathBuf;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/{filename}.svg"));
}

fn main() {
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        1.0,
        1.25 * PI,
        0.75 * PI,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();
    let s = Segment::from(arc);

    // Approximation by one segment
    let iter = s.polygonize(SegmentPolygonizer::MaximumSegmentLength(10.0));
    let pts: Vec<_> = iter.collect();
    let poly_one_seg = SegmentChain::from_points(&pts);

    // Approximation by two segments
    let iter = s.polygonize(SegmentPolygonizer::MaximumAngle(0.5 * PI));
    let pts: Vec<_> = iter.collect();
    let poly_two_seg = SegmentChain::from_points(&pts);

    // Approximation by three segments
    let iter = s.polygonize(SegmentPolygonizer::NumberSegments(3));
    let pts: Vec<_> = iter.collect();
    let poly_three_seg = SegmentChain::from_points(&pts);

    let mut bb = BoundingBox::from(&s);
    bb.scale(1.05);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let mut style_arc = Style::default();
    style_arc.line_width = 2.0;

    let mut style_poly = Style::default();
    style_poly.line_width = 1.0;
    style_poly.line_style = LineStyle::default_dashed();

    let fp = file_path("polygonized_arc");
    view.write_to_file(&fp, move |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        s.draw(&style_arc, cr)?;

        style_poly.line_color = Color {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        poly_one_seg.draw(&style_poly, cr)?;

        style_poly.line_color = Color {
            r: 0.0,
            g: 1.0,
            b: 0.0,
            a: 1.0,
        };
        poly_two_seg.draw(&style_poly, cr)?;

        style_poly.line_color = Color {
            r: 0.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        };
        poly_three_seg.draw(&style_poly, cr)?;

        return Ok(());
    })
    .expect("image creation failed");
}
