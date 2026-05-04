use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;
use std::path::PathBuf;

fn file_path(filename: &str) -> PathBuf {
    return PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("docs/img/{filename}.svg"));
}

fn main() {
    let points = &[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, -0.5],
        [0.5, -0.5],
        [0.5, -1.0],
        [0.0, -1.0],
    ];
    let radii = &[0.5, 0.0, 0.25, 2.0];
    let underlying = Polysegment::from_points(points);
    let iter = Segment::fillet_chain(points, radii);
    let fillet_chain = Polysegment::from_iter(iter);

    let mut underlying_style = Style::default();
    underlying_style.line_color = Color {
        r: 0.3,
        g: 0.3,
        b: 0.3,
        a: 1.0,
    };
    underlying_style.line_width = 2.0;
    underlying_style.line_style = LineStyle::default_dashed();
    let mut fillet_chain_style = Style::default();
    fillet_chain_style.line_width = 2.0;

    let mut bb = BoundingBox::from_points(points.iter().cloned()).expect("has points");
    bb.scale(1.01);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(600));

    let fp = file_path("fillet_chain");
    view.write_to_file(&fp, |cr| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.paint()?;

        underlying.draw(&underlying_style, cr)?;
        fillet_chain.draw(&fillet_chain_style, cr)?;

        return Ok(());
    })
    .expect("image creation failed");
}
