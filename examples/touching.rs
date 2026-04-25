use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::{f64::consts::FRAC_PI_2, path::PathBuf};

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/img/{filename}.svg"));
}

fn main() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let l1_not_touching = LineSegment::new([0.0, 0.0], [2.0, 0.0], e, m).unwrap();
    let l2_not_touching = LineSegment::new([0.0, 1.0], [2.0, 0.3], e, m).unwrap();

    let text_not_touching = Text::new(
        "Not touching".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [1.0, -1.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // =========================================================================

    let l1_touching = LineSegment::new([3.0, 0.0], [5.0, 0.0], e, m).unwrap();
    let l2_touching = LineSegment::new([3.0, 1.0], [4.0, 0.0], e, m).unwrap();

    let text_touching = Text::new(
        "Touching".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [4.0, -1.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // =========================================================================

    let l1_dividing = LineSegment::new([6.0, 0.0], [8.0, 0.0], e, m).unwrap();
    let l2_dividing = LineSegment::new([6.0, 1.0], [8.0, -0.25], e, m).unwrap();

    let text_dividing = Text::new(
        "Dividing".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [7.0, -1.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // =========================================================================

    let l_tangent = LineSegment::new([9.0, 0.0], [11.0, 0.0], e, m).unwrap();
    let a_tangent = ArcSegment::from_center_radius_start_offset_angle(
        [10.0, 1.5],
        1.5,
        -0.5 * FRAC_PI_2,
        -FRAC_PI_2,
        e,
        m,
    )
    .unwrap();

    let text_tangent = Text::new(
        "Touching".into(),
        Anchor::Bottom,
        [0.0, 0.0],
        [10.0, -1.0],
        Color::new(0.0, 0.0, 0.0, 1.0),
        16.0,
        0.0,
    );

    // =========================================================================

    let mut style = Style::default();
    style.line_width = 2.0;
    style.background_color = Color::new(1.0, 0.2, 0.2, 1.0);

    let bb = BoundingBox::new(-0.5, 11.5, -1.2, 1.2);
    let view = Viewport::from_bounding_box(&bb, SideLength::Width(600));

    let draw_fn = |cr: &cairo::Context| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        l1_not_touching.draw(&style, cr)?;
        l2_not_touching.draw(&style, cr)?;
        text_not_touching.draw(cr)?;

        l1_touching.draw(&style, cr)?;
        l2_touching.draw(&style, cr)?;
        text_touching.draw(cr)?;

        l1_dividing.draw(&style, cr)?;
        l2_dividing.draw(&style, cr)?;
        text_dividing.draw(cr)?;

        l_tangent.draw(&style, cr)?;
        a_tangent.draw(&style, cr)?;
        text_tangent.draw(cr)?;
        return Ok(());
    };

    let fp = file_path("touching");
    view.write_to_file(fp, draw_fn)
        .expect("image could not be created");
}
