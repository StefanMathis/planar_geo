use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::path::PathBuf;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/img/{filename}.svg"));
}

fn main() {
    // Grid
    let bb = BoundingBox::new(0.0, 5.0, 0.0, 5.0);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(300));

    let mut grid_style = Style::default();
    grid_style.line_color = Color {
        r: 0.3,
        g: 0.3,
        b: 0.3,
        a: 1.0,
    };

    let fp = file_path("example_grid");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        return grid(&bb, 1.0, [0.0, 0.0], &grid_style, cr);
    })
    .expect("image creation failed");

    let bb = BoundingBox::new(-0.25, 1.25, -1.25, 0.25);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    let mut cs_style = Style::default();
    cs_style.line_color = Color {
        r: 1.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    cs_style.background_color = Color {
        r: 1.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    cs_style.line_width = 2.0;
    cs_style.line_cap = cairo::LineCap::Square;

    let fp = file_path("example_cs");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;
        return coordinate_system(
            [0.0, 0.0],
            1.0,
            -1.0,
            ArrowHeadSize::Height(0.25),
            &cs_style,
            cr,
        );
    })
    .expect("image creation failed");
}
