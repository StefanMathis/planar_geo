use std::f64::consts::PI;

use cairo_viewport::{SideLength, Viewport};
use planar_geo::{
    contour::ArrowHeadSize,
    prelude::*,
    visualize::{Color, LineStyle, Style},
};

#[test]
fn test_horizontal_arrow() {
    let path =
        Contour::arrow_from_head_length_angle([0.0, 0.0], 1.0, 0.0, ArrowHeadSize::Height(0.5))
            .unwrap();
    let mut style = Style::default();
    style.line_style = LineStyle::Solid;
    style.background_color = Color::new(0.0, 0.0, 0.0, 1.0);
    let arrow = Shape::new(vec![path]);

    let view = Viewport::from_bounding_box(&(&arrow).into(), SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/horizontal_arrow.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                arrow.draw(&style, cr)
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_arrow_with_dashed_tail() {
    let angle = 20.0 / 180.0 * PI;
    let path =
        Contour::arrow_from_head_length_angle([0.0, 0.0], 1.0, angle, ArrowHeadSize::Height(0.1))
            .unwrap();
    let mut style = Style::default();
    style.line_style = LineStyle::default_dashed();
    let arrow = Shape::new(vec![path]);

    let segment_chain = Contour::rectangle([0.0, 0.0], [1.0, 1.0]);
    let rectangle = Shape::new(vec![segment_chain]);

    let view = Viewport::from_bounded_entities(
        [BoundingBox::from(&arrow), BoundingBox::from(&rectangle)].into_iter(),
        SideLength::Long(500),
    )
    .unwrap();

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/arrow_with_dashed_tail.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                arrow.draw(&style, cr)?;
                rectangle.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}
