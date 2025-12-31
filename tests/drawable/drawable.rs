use std::f64::consts::{FRAC_PI_2, TAU};

use bounding_box::BoundingBox;
use cairo_viewport::*;
use planar_geo::{
    DEFAULT_EPSILON, DEFAULT_MAX_ULPS,
    contour::Contour,
    segment::{ArcSegment, LineSegment, Segment},
    segment_chain::SegmentChain,
    shape::Shape,
    visualize::{Anchor, Color, LineStyle, Style, Text},
};

#[test]
fn test_segment_line() {
    let line: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
        .unwrap()
        .into();
    let mut style = Style::default();
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 3.0;
    style.line_style = LineStyle::Solid;

    let view = Viewport::from_bounding_box(
        &BoundingBox::new(-0.1, 1.1, -0.5, 0.5),
        SideLength::Long(500),
    );

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/line.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                line.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_segment_arc() {
    let line: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        1.0,
        0.0,
        FRAC_PI_2,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap()
    .into();

    let mut style = Style::default();
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 3.0;
    style.line_style = LineStyle::Solid;

    let view = Viewport::from_bounding_box(
        &BoundingBox::new(-0.1, 1.1, -0.1, 1.1),
        SideLength::Long(500),
    );

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/arc.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                line.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_text() {
    let mut text = Text::new(
        "Test".into(),
        Anchor::Center,
        [250.0, 250.0],
        [0.0, 0.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        100.0,
        0.0,
    );

    let view = Viewport::new([0.0, 0.0], 1.0, 500, 500);

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/text_corners.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                text.anchor = Anchor::TopLeft;
                text.draw(cr)?;
                text.anchor = Anchor::BottomLeft;
                text.draw(cr)?;
                text.anchor = Anchor::TopRight;
                text.draw(cr)?;
                text.anchor = Anchor::BottomRight;
                text.draw(cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    for anchor in [Anchor::Center, Anchor::Centroid].into_iter() {
        text.anchor = anchor;
        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/text_center.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;
                    text.draw(cr)?;
                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );
    }

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/text_edges.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                text.anchor = Anchor::Top;
                text.draw(cr)?;
                text.anchor = Anchor::Left;
                text.draw(cr)?;
                text.anchor = Anchor::Right;
                text.draw(cr)?;
                text.anchor = Anchor::Bottom;
                text.draw(cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_scaled_text() {
    let mut text = Text::new(
        "Test".into(),
        Anchor::Center,
        [250.0, 250.0],
        [0.0, 0.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        30.0,
        0.0,
    );

    let view = Viewport::new([0.0, 0.0], 1.0, 500, 500);
    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/scale_1_text_corners.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                text.anchor = Anchor::TopLeft;
                text.draw(cr)?;
                text.anchor = Anchor::BottomLeft;
                text.draw(cr)?;
                text.anchor = Anchor::TopRight;
                text.draw(cr)?;
                text.anchor = Anchor::BottomRight;
                text.draw(cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    let view = Viewport::new([0.0, 0.0], 2.0, 500, 500);
    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/scale_2_text_corners.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                text.anchor = Anchor::TopLeft;
                text.draw(cr)?;
                text.anchor = Anchor::BottomLeft;
                text.draw(cr)?;
                text.anchor = Anchor::TopRight;
                text.draw(cr)?;
                text.anchor = Anchor::BottomRight;
                text.draw(cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

/**
This is a regression test due to a bug which was found in the WindingExplorer app
 */
#[test]
fn test_shape_and_text_separated() {
    let vertices = vec![[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]];
    let radii = vec![50.0, 0.0, 10.0];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, radii, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );

    let shape = Shape::new(vec![contour]);
    let style = Style::new(
        Color::new(0.0, 0.0, 0.0, 1.0),
        Color::new(0.0, 0.0, 1.0, 0.0),
        1.0,
        LineStyle::Solid,
        cairo::LineCap::Round,
        cairo::LineJoin::Miter,
        None,
    );

    let text = Text::new(
        "Test".into(),
        Anchor::Center,
        [0.0, 0.0],
        [0.0, 0.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        20.0,
        0.0,
    );

    let view = Viewport::from_bounded_entity(&shape, SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/shape_and_text_separated.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                text.draw(cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_shape_scaled_text() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let vertices = vec![[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]];
    let radii = vec![50.0, 0.0, 10.0];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, radii, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );
    let mut style = Style::new(
        Color::new(0.0, 0.0, 1.0, 1.0),
        Color::new(0.0, 0.0, 1.0, 1.0),
        1.0,
        LineStyle::Solid,
        cairo::LineCap::Round,
        cairo::LineJoin::Miter,
        None,
    );
    let shape = Shape::new(vec![contour]);
    let view = Viewport::from_bounded_entity(&shape, SideLength::Long(500));

    // Plot without text
    // ====================================================================================

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_no_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Plot with text
    // ====================================================================================

    // Centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Center,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_centered_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Centroid text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Centroid,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_centroid_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Top-right text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopRight,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_right_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Right text with offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Right,
        [0.0, 0.0],
        [-2.0, 2.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_right_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Left text with vertical offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Left,
        [0.0, 0.0],
        [0.0, 2.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_left_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Top left text with vertical offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopLeft,
        [0.0, 0.0],
        [5.0, 5.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_left_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Bottom text with vertical offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_bottom_text.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

fn setup_rotated_text() -> (Shape, Style, Viewport) {
    let vertices = vec![[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]];
    let radii = vec![50.0, 0.0, 0.0];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, radii, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );
    let style = Style::new(
        Color::new(0.0, 0.0, 1.0, 1.0),
        Color::new(0.0, 0.0, 1.0, 1.0),
        1.0,
        LineStyle::Solid,
        cairo::LineCap::Round,
        cairo::LineJoin::Miter,
        None,
    );
    let shape = Shape::new(vec![contour]);
    let config = Viewport::from_bounded_entity(&shape, SideLength::Long(500));
    return (shape, style, config);
}

#[test]
fn test_rotated_text_center() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let (shape, mut style, view) = setup_rotated_text();

    // Unrotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Center,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_center_00_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Rotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Center,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::PI / 4.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_center_45_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Rotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Center,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::FRAC_PI_2,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_center_90_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_rotated_text_top_left() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let (shape, mut style, view) = setup_rotated_text();

    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopLeft,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_left_0_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Rotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopLeft,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::PI / 4.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_left_45_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_rotated_text_top_right() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let (shape, mut style, view) = setup_rotated_text();

    // Rotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopRight,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_right_0_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Rotated centered text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopRight,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::PI / 4.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_right_45_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_rotated_text_right() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let (shape, mut style, view) = setup_rotated_text();

    // With offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Right,
        [0.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::FRAC_PI_2,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_right_90_deg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // With offset
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::Right,
        [30.0, 0.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        std::f64::consts::FRAC_PI_2,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_right_90_deg_with_offset.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_unscaled_text_and_anchor() {
    let black = Color::new(0.0, 0.0, 0.0, 1.0);
    let vertices = vec![[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]];
    let radii = vec![50.0, 0.0, 10.0];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, radii, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );
    let mut style = Style::new(
        Color::new(0.0, 0.0, 1.0, 1.0),
        Color::new(0.0, 0.0, 1.0, 1.0),
        1.0,
        LineStyle::Solid,
        cairo::LineCap::Round,
        cairo::LineJoin::Miter,
        None,
    );
    let shape = Shape::new(vec![contour]);
    let view = Viewport::from_bounded_entity(&shape, SideLength::Long(500));

    // Top-right text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopRight,
        [-10.0, 10.0],
        [0.0, 0.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_right_text_unscaled.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Top-right text
    style.text = Some(Box::new(Text::new(
        "Test".into(),
        Anchor::TopRight,
        [0.0, 0.0],
        [-10.0, 10.0],
        black.clone(),
        20.0,
        0.0,
    )));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/triangle_top_right_text_unscaled_offset_scaled.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_circle_in_square() {
    let verts = vec![[1.0, 1.0], [1.0, -1.0], [-1.0, -1.0], [-1.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(
        &verts,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    ));

    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        0.5,
        0.0,
        TAU,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap()
    .into();
    let hole = Contour::new(SegmentChain::new(
        vec![segment],
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    ));
    let shape = Shape::from_contour_and_holes(contour, Some(vec![hole]));

    // Set the background_color to blue
    let mut style = Style::default();
    style.background_color = Color::new(0.0, 0.0, 1.0, 1.0);
    let view = Viewport::from_bounded_entity(&shape, SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/circle_in_square.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_two_shapes() {
    let verts = vec![[1.0, 0.0], [1.0, -1.0], [-1.0, -1.0], [-1.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(
        &verts,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    ));
    let shape_outer = Shape::from_contour_and_holes(contour, None);

    // Set the background_color to blue
    let mut style_inner = Style::default();
    style_inner.background_color = Color::new(0.0, 0.0, 1.0, 1.0);

    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        0.5,
        0.0,
        TAU,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap()
    .into();
    let contour = Contour::new(SegmentChain::new(
        vec![segment],
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    ));
    let shape_inner = Shape::from_contour_and_holes(contour, None);

    // Set the background_color to green
    let mut style_outer = Style::default();
    style_outer.background_color = Color::new(0.0, 1.0, 0.0, 1.0);

    let view = Viewport::from_bounded_entities(
        [
            BoundingBox::from(&shape_outer),
            BoundingBox::from(&shape_inner),
        ]
        .into_iter(),
        SideLength::Long(500),
    )
    .unwrap();

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/circle_in_square_different_color.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape_outer.draw(&style_outer, cr)?;
                shape_inner.draw(&style_inner, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_quarter_arc() {
    // Positive circle
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        0.5,
        FRAC_PI_2,
        FRAC_PI_2,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap()
    .into();
    let segment_chain = SegmentChain::new(vec![segment], DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    let mut style = Style::default();
    style.line_width = 3.0;

    let view =
        Viewport::from_bounding_box(&BoundingBox::from(&segment_chain), SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/quarter_arc_pos.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                segment_chain.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Negative circle
    let segment: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        0.5,
        -FRAC_PI_2,
        -FRAC_PI_2,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap()
    .into();
    let contour = Contour::new(SegmentChain::new(
        vec![segment],
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    ));
    let shape = Shape::from_contour_and_holes(contour, None);
    style.line_width = 3.0;

    let view = Viewport::from_bounded_entity(&shape, SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/quarter_arc_neg.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_block_with_fillets() {
    let vertices = vec![[10f64, 0.0], [10f64, 10f64], [-10f64, 10f64], [-10f64, 0.0]];
    let fillets = vec![1f64, 1f64, 1f64, 1f64];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, fillets, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );

    let shape = Shape::from_contour_and_holes(contour, None);

    let mut style = Style::default();
    style.background_color = Color::new(1.0, 0.0, 0.0, 1.0);

    let mut bb = BoundingBox::from(&shape);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    approx::assert_abs_diff_eq!(view.scale, 22.7272727, epsilon = 0.00001);
    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/block_with_fillets_1.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Now with another scale
    let vertices = vec![[10e3, 0.0], [10e3, 10e3], [-10e3, 10e3], [-10e3, 0.0]];
    let fillets = vec![1e3, 1e3, 1e3, 1e3];
    let contour = Contour::new(
        SegmentChain::from_fillets(vertices, fillets, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap(),
    );

    let shape = Shape::from_contour_and_holes(contour, None);
    style.background_color = Color::new(1.0, 0.0, 0.0, 1.0);

    let mut bb = BoundingBox::from(&shape);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    approx::assert_abs_diff_eq!(view.scale, 0.0227272727, epsilon = 0.00001);
    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/block_with_fillets_2.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    // Now with another scale
    let vertices = vec![[10e-3, 0.0], [10e-3, 10e-3], [-10e-3, 10e-3], [-10e-3, 0.0]];
    let fillets = vec![1e-3, 1e-3, 1e-3, 1e-3];
    let contour =
        SegmentChain::from_fillets(vertices, fillets, true, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();

    let shape = Shape::from_contour_and_holes(contour, None);
    style.background_color = Color::new(1.0, 0.0, 0.0, 1.0);

    let mut bb = BoundingBox::from(&shape);
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    approx::assert_abs_diff_eq!(view.scale, 22727.272727, epsilon = 0.00001);
    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/block_with_fillets_3.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;
                shape.draw(&style, cr)?;
                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_line_style() {
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

    view.write_to_file("docs/line_styles.svg", create_fn)
        .unwrap();

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/line_styles.png"),
            create_fn,
            0.99
        )
        .is_ok()
    );
}
