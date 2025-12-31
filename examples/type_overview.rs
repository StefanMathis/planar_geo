use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::{f64::consts::PI, path::PathBuf};

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/{filename}.svg"));
}

fn main() {
    // Abbreviations for tolerances to keep the subsequent function calls clear.
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let line = LineSegment::new([0.0, -5.0], [2.0, -5.0], e, m).unwrap();

    let line_text = Text::new(
        "LineSegment".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [1.0, -4.5],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    let arc = ArcSegment::from_center_radius_start_offset_angle(
        [1.0, -2.0],
        1.4,
        -0.25 * PI,
        -0.5 * PI,
        e,
        m,
    )
    .unwrap();

    let arc_text = Text::new(
        "ArcSegment".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [1.0, -2.5],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // =========================================================================

    let mut segment_chain = SegmentChain::new();
    segment_chain.push_back(
        LineSegment::new([0.0, -4.0], [-1.0, -4.0], e, m)
            .unwrap()
            .into(),
    );
    segment_chain.push_back(
        ArcSegment::from_start_center_angle([-1.0, -4.0], [0.0, -4.0], 0.5 * PI, e, m)
            .unwrap()
            .into(),
    );
    segment_chain.push_back(
        LineSegment::new([1.0, -4.0], [0.0, -3.0], e, m)
            .unwrap()
            .into(),
    );
    segment_chain.translate([5.0, 0.0]);

    let segment_chain_text = Text::new(
        "SegmentChain".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [5.0, -2.5],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    let mut c_chain = SegmentChain::new();
    c_chain.push_back(
        LineSegment::new([0.0, -4.0], [0.0, -5.2], e, m)
            .unwrap()
            .into(),
    );
    c_chain.push_back(
        ArcSegment::from_start_middle_stop([0.0, -5.2], [1.3, -4.0], [-1.3, -4.0], e, m)
            .unwrap()
            .into(),
    );
    let mut contour = Contour::new(c_chain);
    contour.translate([9.2, 0.0]);

    let contour_text = Text::new(
        "Contour".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [9.2, -2.5],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // Define the shape
    let mut s_chain = SegmentChain::new();
    s_chain.push_back(
        ArcSegment::from_start_center_angle([1.0, -4.0], [0.0, -4.0], 0.5 * PI, e, m)
            .unwrap()
            .into(),
    );
    s_chain.extend_back([-1.0, -3.0]);
    s_chain.extend_back([-1.0, -5.0]);
    s_chain.extend_back([1.0, -5.0]);

    let outer: Contour = s_chain.into();
    let hole: Contour =
        SegmentChain::from_points(&[[-0.5, -3.5], [-0.5, -4.5], [0.5, -4.5]]).into();
    let mut shape = Shape::new(vec![outer, hole]).expect("valid inputs");
    shape.translate([13.0, 0.0]);

    let shape_text = Text::new(
        "Shape".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [13.0, -2.5],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    let mut style = Style::default();
    style.line_width = 2.0;
    style.background_color = Color::new(1.0, 0.2, 0.2, 1.0);

    let bb = BoundingBox::new(-0.5, 14.1, -5.5, -2.0);
    let view = Viewport::from_bounding_box(&bb, SideLength::Width(600));

    let draw_fn = |cr: &cairo::Context| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        line.draw(&style, cr)?;
        line_text.draw(cr)?;

        arc.draw(&style, cr)?;
        arc_text.draw(cr)?;

        segment_chain.draw(&style, cr)?;
        segment_chain_text.draw(cr)?;

        contour.draw(&style, cr)?;
        contour_text.draw(cr)?;

        shape.draw(&style, cr)?;
        shape_text.draw(cr)?;
        return Ok(());
    };

    let fp = file_path("type_overview");
    view.write_to_file(fp, draw_fn)
        .expect("image could not be created");
}
