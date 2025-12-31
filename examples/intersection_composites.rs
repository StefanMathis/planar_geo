use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;

fn main() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(SegmentChain::from_points(&vertices));

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(SegmentChain::from_points(&vertices));

    let shape = Shape::new(vec![c1, c2]).unwrap();

    let chain = SegmentChain::from_points(&[[-1.0, 1.0], [-0.5, 0.5], [1.5, 0.5], [2.0, 1.0]]);

    let view = Viewport::from_bounding_box(
        &BoundingBox::new(-1.2, 2.2, -0.1, 1.1),
        SideLength::Long(500),
    );

    let mut style = Style::default();
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 2.0;
    style.line_style = LineStyle::Solid;
    style.background_color = Color::from_rgba8(144, 213, 255, 255);

    let mut intersected_segments_style = Style::default();
    intersected_segments_style.line_color = Color::new(1.0, 1.0, 0.0, 1.0);
    intersected_segments_style.line_width = 3.0;
    intersected_segments_style.line_style = LineStyle::Solid;

    let intersection_style = IntersectionStyle::default();

    let draw_fn = |cr: &cairo::Context| {
        // Set the background to white
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        shape.draw(&style, cr)?;
        chain.draw(&style, cr)?;
        for i in shape.intersections_segment_chain(&chain, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
            i.draw(
                &intersection_style,
                Some((&shape, &intersected_segments_style)),
                Some((&chain, &intersected_segments_style)),
                cr,
            )?;
        }

        return Ok(());
    };

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(&format!("docs/intersection_composites.svg"));

    view.write_to_file(path, draw_fn)
        .expect("image could not be created");
}
