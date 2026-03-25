use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::path::PathBuf;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/img/{filename}.svg"));
}

fn main() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let start = [0.0, 2.0];
    let stop = [2.0, 0.0];
    let radius: f64 = 2.0;

    // Counter-clockwise small arc (positive)
    let mut red_arc =
        ArcSegment::from_start_stop_radius(start, stop, radius, true, false, e, m).unwrap();
    red_arc.line_reflection([0.0, 0.0], [1.0, 0.0]);
    let mut red_style = Style::default();
    red_style.line_color = Color::new(1.0, 0.0, 0.0, 1.0);
    red_style.line_width = 2.0;

    // Counter-clockwise large arc (positive)
    let mut blue_arc =
        ArcSegment::from_start_stop_radius(start, stop, radius, true, true, e, m).unwrap();
    blue_arc.line_reflection([0.0, 0.0], [1.0, 0.0]);
    let mut blue_style = Style::default();
    blue_style.line_color = Color::new(0.0, 0.0, 1.0, 1.0);
    blue_style.line_width = 2.0;

    let mut green_arc =
        ArcSegment::from_start_stop_radius(start, stop, radius, false, false, e, m).unwrap();
    green_arc.line_reflection([0.0, 0.0], [1.0, 0.0]);
    let mut green_style = Style::default();
    green_style.line_color = Color::new(0.0, 1.0, 0.0, 1.0);
    green_style.line_width = 2.0;

    // Clockwise small arc (negative)
    let mut yellow_arc =
        ArcSegment::from_start_stop_radius(start, stop, radius, false, true, e, m).unwrap();
    yellow_arc.line_reflection([0.0, 0.0], [1.0, 0.0]);
    let mut yellow_style = Style::default();
    yellow_style.line_color = Color::new(1.0, 1.0, 0.0, 1.0);
    yellow_style.line_width = 2.0;

    let bb = BoundingBox::new(-2.1, 4.1, -4.1, 2.1);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("example_four_arcs");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        red_arc.draw(&red_style, cr)?;
        blue_arc.draw(&blue_style, cr)?;
        green_arc.draw(&green_style, cr)?;
        yellow_arc.draw(&yellow_style, cr)?;
        return Ok(());
    })
    .expect("image creation failed");
}
