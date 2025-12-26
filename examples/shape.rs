use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};
use std::path::PathBuf;

pub const UNSCALED_FONT_SIZE: f64 = 16.0;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/{filename}.svg"));
}

fn main() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let first_arc =
        ArcSegment::from_center_radius_start_offset_angle([1.5, 0.0], 1.5, PI, -FRAC_PI_2, e, m)
            .expect("radius is positive and offset angle is not zero");
    let line = LineSegment::new([1.5, 1.5], [3.5, 1.5], e, m).expect("segment length is not zero");
    let second_arc = ArcSegment::from_center_radius_start_offset_angle(
        [3.5, 0.0],
        1.5,
        FRAC_PI_2,
        -FRAC_PI_2,
        e,
        m,
    )
    .expect("radius is positive and offset angle is not zero");

    // Connect those three segments to a segment chain
    let mut chain = SegmentChain::new();
    chain.push_back(first_arc.into());
    chain.push_back(line.into());
    chain.push_back(second_arc.into());

    // Create a contour out of the chain. If start and end of the chain are not
    // identical, a line segment is automatically added
    let outer_contour = Contour::new(chain);

    // Create a second contour via a simplified constructor
    let rect_hole = Contour::new(SegmentChain::from_points(&[
        [1.5, 0.2],
        [3.5, 0.2],
        [3.5, 1.3],
        [1.5, 1.3],
    ]));

    // First element of the vector is interpreted as outer contour, all further
    // elements are holes
    let shape = Shape::new(vec![outer_contour, rect_hole]).unwrap();

    // Drawing
    // =========================================================================

    let mut shape_draw = shape;
    shape_draw.line_reflection([0.0, 0.0], [1.0, 0.0]);

    let mut style = Style::default();
    style.background_color = Color::from_rgba8(144, 213, 255, 255);
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 2.0;

    let mut bb = BoundingBox::from(&shape_draw);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_shape");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        return shape_draw.draw(&style, cr);
    })
    .expect("image creation failed");
}
