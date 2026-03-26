use bounding_box::BoundingBox;
use cairo_viewport::*;
use planar_geo::prelude::*;

#[test]
fn test_grid() {
    let mut grid_style = Style::default();
    grid_style.line_color = Color {
        r: 0.3,
        g: 0.3,
        b: 0.3,
        a: 1.0,
    };

    {
        let bb = BoundingBox::new(0.0, 10.0, 0.0, 10.0);
        let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 1.0, [0.0, 0.0], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid_small_spacing.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 0.5, [0.0, 0.0], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid_offset.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 1.0, [0.5, 0.5], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid_offset.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 1.0, [1.5, 1.5], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );
    }
    {
        let bb = BoundingBox::new(-0.5, 4.5, -0.5, 4.5);
        let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid_small.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 1.0, [0.0, 0.0], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );
    }
    {
        let bb = BoundingBox::new(-0.5, 9.5, -0.5, 9.5);
        let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

        assert!(
            view.compare_or_create(
                std::path::Path::new("tests/img/grid_offset.png"),
                |cr| {
                    // Set the background to white
                    cr.set_source_rgb(1.0, 1.0, 1.0);
                    cr.paint()?;

                    // Draw the grid
                    grid(&bb, 1.0, [0.0, 0.0], &grid_style, cr)?;

                    return Ok(());
                },
                0.99
            )
            .is_ok()
        );
    }
}

#[test]
fn test_coordinate_system() {
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

    let bb = BoundingBox::new(-2.0, 2.0, -2.0, 2.0);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/cs_central.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;

                // Draw the coordinate system
                coordinate_system(
                    [0.0, 0.0],
                    1.0,
                    1.0,
                    ArrowHeadSize::Height(0.25),
                    &cs_style,
                    cr,
                )?;

                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/cs_offset.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;

                // Draw the coordinate system
                coordinate_system(
                    [0.5, 0.5],
                    1.0,
                    1.0,
                    ArrowHeadSize::Height(0.25),
                    &cs_style,
                    cr,
                )?;

                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/cs_inverted_x.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;

                // Draw the coordinate system
                coordinate_system(
                    [0.0, 0.0],
                    -1.0,
                    1.0,
                    ArrowHeadSize::Height(0.25),
                    &cs_style,
                    cr,
                )?;

                return Ok(());
            },
            0.99
        )
        .is_ok()
    );

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/cs_inverted_y.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;

                // Draw the coordinate system
                coordinate_system(
                    [0.0, 0.0],
                    1.0,
                    -1.0,
                    ArrowHeadSize::Height(0.25),
                    &cs_style,
                    cr,
                )?;

                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}

#[test]
fn test_grid_and_coordinate_system() {
    let mut grid_style = Style::default();
    grid_style.line_color = Color {
        r: 0.3,
        g: 0.3,
        b: 0.3,
        a: 1.0,
    };
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

    let bb = BoundingBox::new(-5.0, 5.0, -5.0, 5.0);
    let view = Viewport::from_bounding_box(&bb, SideLength::Long(500));

    assert!(
        view.compare_or_create(
            std::path::Path::new("tests/img/grid_and_cs.png"),
            |cr| {
                // Set the background to white
                cr.set_source_rgb(1.0, 1.0, 1.0);
                cr.paint()?;

                // Draw the grid
                grid(&bb, 1.0, [0.0, 0.0], &grid_style, cr)?;

                // Draw the coordinate system
                coordinate_system(
                    [0.0, 0.0],
                    2.0,
                    -2.0,
                    ArrowHeadSize::Height(0.5),
                    &cs_style,
                    cr,
                )?;

                return Ok(());
            },
            0.99
        )
        .is_ok()
    );
}
