use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::f64::consts::PI;

fn main() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let mut line_1: Segment = LineSegment::new([0.0, 0.0], [3.0, 0.0], e, m)
        .expect("segment length is not zero")
        .into();
    let mut line_2: Segment = LineSegment::new([2.0, -0.5], [2.0, 0.5], e, m)
        .expect("segment length is not zero")
        .into();
    let mut arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
        [1.0, 0.0],
        0.5,
        -0.25 * PI,
        1.5 * PI,
        e,
        m,
    )
    .expect("offet angle is not zero")
    .into();

    line_1.line_reflection([0.0, 0.0], [1.0, 0.0]);
    line_2.line_reflection([0.0, 0.0], [1.0, 0.0]);
    arc.line_reflection([0.0, 0.0], [1.0, 0.0]);

    let mut style = Style::default();
    style.line_color = Color::new(0.0, 0.0, 0.0, 1.0);
    style.line_width = 2.0;

    let mut pt_style = Style::default();
    pt_style.line_color = Color::new(0.0, 0.0, 1.0, 1.0);
    pt_style.line_width = 2.0;

    let mut bb = BoundingBox::from(&line_1)
        .union(&BoundingBox::from(&line_2))
        .union(&BoundingBox::from(&arc));
    bb.scale(1.1);

    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let draw_fn = |cr: &cairo::Context| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        line_1.draw(&style, cr)?;
        let text = Text::new(
            "line_1".into(),
            Anchor::BottomRight,
            [2.0, 2.0],
            [0.0, 0.0],
            Color::new(0.0, 0.0, 0.0, 1.0),
            16.0,
            0.0,
        );
        text.draw(cr)?;

        line_2.draw(&style, cr)?;
        let text = Text::new(
            "line_2".into(),
            Anchor::BottomLeft,
            [-2.0, 2.0],
            line_2.stop(),
            Color::new(0.0, 0.0, 0.0, 1.0),
            16.0,
            0.0,
        );
        text.draw(cr)?;

        arc.draw(&style, cr)?;
        let text = Text::new(
            "arc".into(),
            Anchor::BottomRight,
            [2.0, 2.0],
            arc.stop(),
            Color::new(0.0, 0.0, 0.0, 1.0),
            16.0,
            0.0,
        );
        text.draw(cr)?;

        // Points
        let radius = 0.02;
        let c = Contour::circle([2.3, 0.0], radius);
        c.draw(&pt_style, cr)?;
        let c = Contour::circle([2.0, -0.1], radius);
        c.draw(&pt_style, cr)?;
        let c = Contour::circle([1.0, -0.5], radius);
        c.draw(&pt_style, cr)?;

        let i_style = IntersectionStyle::default();

        let i: Intersection<(), ()> = [1.5, 0.0].into();
        i.draw(&i_style, cr)?;

        let i: Intersection<(), ()> = [2.0, 0.0].into();
        i.draw(&i_style, cr)?;

        let i: Intersection<(), ()> = [0.5, 0.0].into();
        i.draw(&i_style, cr)?;
        return Ok(());
    };

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(&format!("docs/intersection_segments.svg"));

    view.write_to_file(path, draw_fn)
        .expect("image could not be created");
}
