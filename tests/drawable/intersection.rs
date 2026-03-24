use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;

#[test]
fn test_intersection_visualization() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut contours = vec![Contour::new(Polysegment::from_points(&vertices))];

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    contours.push(Contour::new(Polysegment::from_points(&vertices)));

    let shape = Shape::new(contours).unwrap();

    let polysegment = Polysegment::from_points(&[
        [-1.0, 1.0],
        [-1.0, 0.5],
        [2.0, 0.5],
        [0.5, -1.0],
        [0.5, 2.0],
        [0.0, 2.0],
    ]);

    let view = Viewport::from_bounding_box(
        &BoundingBox::new(-1.2, 2.2, -1.2, 2.2),
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
        polysegment.draw(&style, cr)?;
        for i in shape.intersections_polysegment(&polysegment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS) {
            i.draw(
                &intersection_style,
                Some(DrawableRef::new(&shape, intersected_segments_style.clone())),
                Some(DrawableRef::new(
                    &polysegment,
                    intersected_segments_style.clone(),
                )),
                cr,
            )?;
        }

        return Ok(());
    };

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/intersection_shape_polysegment.png"),
            draw_fn,
            0.99
        )
        .is_ok()
    );
}
