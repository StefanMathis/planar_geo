use cairo_viewport::{SideLength, Viewport};
use planar_geo::prelude::*;

#[test]
fn test_text_placement() {
    fn create_fn() -> Box<dyn FnOnce(&cairo::Context) -> Result<(), cairo::Error>> {
        let e = DEFAULT_EPSILON;
        let m = DEFAULT_MAX_ULPS;

        let mut style = Style::default();
        style.line_color = Color::new(1.0, 0.5, 0.5, 1.0);
        style.line_width = 2.0;
        let mut txt = Text {
            text: "Center".into(),
            anchor: Anchor::Center,
            fixed_anchor_offset: [0.0, 0.0],
            scaled_anchor_offset: [0.0, 0.0],
            color: Color::new(0.0, 0.0, 0.0, 1.0),
            font_size: 16.0,
            angle: 0.0,
        };

        let contour = Contour::from_points(&[[0.1, 0.1], [1.9, 0.1], [1.9, 0.9], [0.1, 0.9]], e, m);

        // A cross on the right of the contour which shows placement of standalone texts
        let cc = [2.5, 0.5];
        let hori = LineSegment::new([cc[0] - 0.2, cc[1]], [cc[0] + 0.2, cc[1]], e, m)
            .expect("points not identical");
        let vert = LineSegment::new([cc[0], cc[1] - 0.2], [cc[0], cc[1] + 0.2], e, m)
            .expect("points not identical");

        let drawing_fn = move |cr: &cairo::Context| {
            // Set the background to white
            cr.set_source_rgb(1.0, 1.0, 1.0);
            cr.paint()?;

            // Draw the contour first
            contour.draw(&style, cr)?;

            // Draw the cross
            hori.draw(&style, cr)?;
            vert.draw(&style, cr)?;

            // Then draw the contour again with an invisible line color to demonstrate the text placement
            style.line_color = Color::new(0.0, 0.0, 0.0, 0.0);

            let anchors_str = &[
                "Center",
                "TopLeft",
                "Top",
                "TopRight",
                "Right",
                "BottomRight",
                "Bottom",
                "BottomLeft",
                "Left",
            ];
            let anchors = &[
                Anchor::Center,
                Anchor::TopLeft,
                Anchor::Top,
                Anchor::TopRight,
                Anchor::Right,
                Anchor::BottomRight,
                Anchor::Bottom,
                Anchor::BottomLeft,
                Anchor::Left,
            ];
            for (anchor_str, anchor) in anchors_str.iter().zip(anchors.iter()) {
                txt.anchor = *anchor;
                txt.text = anchor_str.to_string();
                style.text = Some(Box::new(txt.clone()));
                contour.draw(&style, cr)?;
            }

            // Show text placement without an associated geometric type
            txt.scaled_anchor_offset = cc;
            let anchors_str = &["TopLeft", "TopRight", "BottomRight", "BottomLeft"];
            let anchors = &[
                Anchor::TopLeft,
                Anchor::TopRight,
                Anchor::BottomRight,
                Anchor::BottomLeft,
            ];
            for (anchor_str, anchor) in anchors_str.iter().zip(anchors.iter()) {
                txt.anchor = *anchor;
                txt.text = anchor_str.to_string();
                txt.draw(cr)?;
            }

            return Ok(());
        };
        return Box::new(drawing_fn);
    }

    let view =
        Viewport::from_bounding_box(&BoundingBox::new(0.0, 3.0, 0.0, 1.0), SideLength::Long(600));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/text_placement.png"),
            create_fn(),
            0.99
        )
        .is_ok()
    );

    view.write_to_file("docs/text_placement.svg", create_fn())
        .unwrap();
}

#[test]
fn anchor_offset_scaling() {
    fn create_fn() -> Box<dyn FnOnce(&cairo::Context) -> Result<(), cairo::Error>> {
        let e = DEFAULT_EPSILON;
        let m = DEFAULT_MAX_ULPS;

        let mut style = Style::default();
        style.line_color = Color::new(1.0, 0.5, 0.5, 1.0);
        style.line_width = 2.0;

        let mut contour =
            Contour::from_points(&[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]], e, m);

        let drawing_fn = move |cr: &cairo::Context| {
            // Set the background to white
            cr.set_source_rgb(1.0, 1.0, 1.0);
            cr.paint()?;

            let txt = Text {
                text: "Scaled offset".into(),
                anchor: Anchor::TopLeft,
                fixed_anchor_offset: [0.0, 0.0],
                scaled_anchor_offset: [0.3, 0.3],
                color: Color::new(0.0, 0.0, 0.0, 1.0),
                font_size: 16.0,
                angle: 0.0,
            };
            style.text = Some(Box::new(txt));
            contour.draw(&style, cr)?;

            let txt = Text {
                text: "Fixed offset".into(),
                anchor: Anchor::TopLeft,
                fixed_anchor_offset: [10.0, 10.0],
                scaled_anchor_offset: [0.0, 0.0],
                color: Color::new(0.0, 0.0, 0.0, 1.0),
                font_size: 16.0,
                angle: 0.0,
            };
            style.line_color = Color::new(0.5, 0.8, 0.5, 1.0);
            style.text = Some(Box::new(txt));
            contour.translate([1.0, 0.0]);
            contour.draw(&style, cr)?;

            return Ok(());
        };
        return Box::new(drawing_fn);
    }

    let mut view =
        Viewport::from_bounding_box(&BoundingBox::new(0.0, 3.0, 0.0, 1.0), SideLength::Long(600));

    view.write_to_file("docs/anchor_offset_scale_1.svg", create_fn())
        .unwrap();

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/anchor_offset_scale_1.png"),
            create_fn(),
            0.99
        )
        .is_ok()
    );

    // Zoom into the image
    view.scale *= 2.0;

    view.write_to_file("docs/anchor_offset_scale_2.svg", create_fn())
        .unwrap();

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/anchor_offset_scale_2.png"),
            create_fn(),
            0.99
        )
        .is_ok()
    );
}
