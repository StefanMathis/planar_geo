use planar_geo::prelude::*;

fn line_and_style() -> (Segment, Style) {
    let line: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
        .unwrap()
        .into();
    return (line, Style::default());
}

#[test]
fn to_drawable() {
    let _ = Drawable::from(line_and_style());
}

#[test]
fn to_drawable_ref() {
    let (line, style) = line_and_style();
    let _ = DrawableRef::from((&line, style));
}

#[test]
fn to_drawable_cow() {
    let _ = DrawableCow::from(line_and_style());
    let (line, style) = line_and_style();
    let _ = DrawableCow::from((&line, style));
}
