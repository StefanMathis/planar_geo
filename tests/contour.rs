use approx::{assert_abs_diff_eq, assert_ulps_eq};
use planar_geo::prelude::*;
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[test]
fn test_from_segments() {
    let mut chain = SegmentChain::new();
    chain.push_back(
        ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, PI, FRAC_PI_2, 0.0, 0)
            .unwrap()
            .into(),
    );
    chain.push_back(
        LineSegment::new([0.0, -2.0], [-2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into(),
    );
    chain.push_back(
        LineSegment::new([-2.0, -2.0], [-2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into(),
    );
    let contour = Contour::new(chain);
    assert_eq!(contour.num_segments(), 3);
}

#[test]
fn test_from_points() {
    // Triangle
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let polygon: Contour = SegmentChain::from_points(vertices).into();
    let vertices: Vec<[f64; 2]> = polygon.points().collect();

    assert_eq!(vertices.len(), 3);
    assert_eq!(vertices[0], [0.0, 0.0]);
    assert_eq!(vertices[1], [1.0, 0.0]);
    assert_eq!(vertices[2], [1.0, 1.0]);
}

#[test]
fn test_rectangle() {
    {
        // This is a line since the y-values are identical
        let contour = Contour::rectangle([0.0, 0.0], [1.0, 0.0]);
        let vertices: Vec<[f64; 2]> = contour.points().collect();
        assert_eq!(vertices.len(), 2);
        assert_eq!(vertices[0], [0.0, 0.0]);
        assert_eq!(vertices[1], [1.0, 0.0]);
    }
    {
        // Construct an identical rectangle from different points
        let contour = Contour::rectangle([0.0, 0.0], [2.0, 2.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );

        let contour = Contour::rectangle([2.0, 2.0], [0.0, 0.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );

        let contour = Contour::rectangle([2.0, 0.0], [0.0, 2.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );

        let contour = Contour::rectangle([0.0, 2.0], [2.0, 0.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );

        let contour = Contour::rectangle([0.0, 2.0], [2.0, 0.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );

        let contour = Contour::rectangle([2.0, 0.0], [0.0, 2.0]);
        assert_eq!(contour.area(), 4.0);
        assert_eq!(
            BoundingBox::from(&contour),
            BoundingBox::new(0.0, 2.0, 0.0, 2.0)
        );
    }
}

#[test]
fn test_cut_on_line() {
    let vertices = [
        [1.0, 0.0],
        [1.0, 1.0],
        [2.0, 1.0],
        [2.0, 2.0],
        [-2.0, 2.0],
        [-2.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
    ];
    let slot_equivalent = Contour::new(SegmentChain::from_points(vertices.as_slice()));

    let cut = SegmentChain::from_points([[-10.0, 1.0], [10.0, 1.0]].as_slice());

    let separated_lines = slot_equivalent.intersection_cut(&cut, 0.0, 0);
    assert_eq!(separated_lines.len(), 4);

    assert_eq!(separated_lines[0].num_segments(), 3);
    for point in separated_lines[0].points() {
        assert!(point[1] <= 1.0);
    }

    assert_eq!(separated_lines[1].num_segments(), 1);
    for point in separated_lines[1].points() {
        assert!(point[1] == 1.0);
    }

    assert_eq!(separated_lines[2].num_segments(), 3);
    for point in separated_lines[2].points() {
        assert!(point[1] >= 1.0);
    }

    assert_eq!(separated_lines[3].num_segments(), 1);
    for point in separated_lines[3].points() {
        assert!(point[1] == 1.0);
    }
}

#[test]
fn test_points() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let polygon = Contour::new(SegmentChain::from_points(vertices));
    let vertices: Vec<[f64; 2]> = polygon.points().collect();
    assert_eq!(vertices.len(), 3);
    assert_eq!(vertices[0], [0.0, 0.0]);
    assert_eq!(vertices[1], [1.0, 0.0]);
    assert_eq!(vertices[2], [0.0, 1.0]);
}

#[test]
fn test_circle() {
    {
        let radius = 2.0;
        let contour = Contour::circle([0.0, 0.0], radius);

        let area_circle = PI * radius.powi(2);
        assert_ulps_eq!(contour.area(), area_circle);

        let bb = BoundingBox::from(&contour);
        assert_ulps_eq!(bb.xmin(), -radius);
        assert_ulps_eq!(bb.ymin(), -radius);
        assert_ulps_eq!(bb.xmax(), radius);
        assert_ulps_eq!(bb.ymax(), radius);
    }
    {
        let contour = planar_geo::prelude::Contour::circle([0.0, 0.0], 20.0);
        assert!(contour.contains_point([0.0, 0.0], 0.0, 0));
        assert!(!contour.contains_point([21.0, 0.0], 0.0, 0));

        // Bug found with the interactive canvas
        assert!(contour.contains_point(
            [0.055555555555542924, -0.9999923706054688],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
        assert!(contour.contains_point(
            [12.055555555555543, -9.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
        assert!(contour.contains_point(
            [8.055555555555543, -6.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
        assert!(contour.contains_point(
            [14.055555555555543, -6.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
    }
}

#[test]
fn test_area() {
    {
        let arc =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, TAU, 0.0, 0)
                .unwrap();
        let contour = Contour::new(arc.into());
        let area_circle = PI * 2.0_f64.powi(2);
        assert_ulps_eq!(contour.area(), area_circle);
    }

    // Square with concave radius in one corner
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::from_center_radius_start_offset_angle(
                [0.0, 0.0],
                2.0,
                PI,
                FRAC_PI_2,
                0.0,
                0,
            )
            .unwrap()
            .into(),
        );
        chain.push_back(
            LineSegment::new([0.0, -2.0], [-2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into(),
        );
        chain.push_back(
            LineSegment::new([-2.0, -2.0], [-2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(chain);
        assert_ulps_eq!(contour.area(), 4.0 - 0.25 * PI * 2.0_f64.powi(2));
    }

    // Square with concave radius in another corner
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::from_center_radius_start_offset_angle(
                [0.0, 0.0],
                2.0,
                0.0,
                -FRAC_PI_2,
                0.0,
                0,
            )
            .unwrap()
            .into(),
        );
        chain.push_back(
            LineSegment::new([0.0, -2.0], [2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into(),
        );
        chain.push_back(
            LineSegment::new([2.0, -2.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(chain);
        assert_ulps_eq!(
            contour.area(),
            4.0 - 0.25 * PI * 2.0_f64.powi(2),
            epsilon = DEFAULT_EPSILON
        );
    }

    // Square with rounded edges
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet([-1.0, 2.0], [-1.0, 0.0], [1.0, 0.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.push_back(
            ArcSegment::fillet([-1.0, 0.0], [1.0, 0.0], [1.0, 2.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 2.0], [-1.0, 2.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.push_back(
            ArcSegment::fillet([1.0, 2.0], [-1.0, 2.0], [-1.0, 0.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(chain);

        assert_ulps_eq!(
            contour.area(),
            4.0 - 2.0 * 0.1f64.powi(2) * (1.0 - 0.25 * PI)
                - 2.0 * 0.2f64.powi(2) * (1.0 - 0.25 * PI),
            epsilon = DEFAULT_EPSILON
        );
    }

    // "Pie" section with a straight connection between the end points of the arc
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 1.0, 0.0, 0)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(chain);
        assert_ulps_eq!(contour.area(), 0.25 * PI - 0.5, epsilon = DEFAULT_EPSILON);
    }

    // "Pie" section with a square cut out at the origin
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet([0.0, -1.0], [-1.0, -1.0], [-1.0, 0.0], 1.0, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.extend_back([-1.0, 0.0]);
        chain.extend_back([-0.5, 0.0]);
        chain.extend_back([-0.5, -0.5]);
        chain.extend_back([0.0, -0.5]);
        chain.extend_back([0.0, -1.0]);
        let contour = Contour::new(chain);
        assert_ulps_eq!(contour.area(), 0.25 * PI - 0.25, epsilon = DEFAULT_EPSILON);
    }
}

#[test]
fn test_contains_point() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(SegmentChain::from_points(vertices));
    assert!(contour.contains_point([0.5, 0.5], 0.0, 0));
    assert!(contour.contains_point([0.5, 1.0], 0.0, 0));
    assert!(!contour.contains_point([0.5, 1.5], 0.0, 0));
}

#[test]
fn test_self_intersection() {
    // Open segment_chain
    {
        let c: Contour =
            SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).into();

        // Intersect the line with itself
        let intersections = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.count(), 0);

        // Intersect the line with itself
        let intersections = c.intersections_contour_par(&c, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.count(), 0);
    }

    // Closed segment_chain
    {
        let c: Contour =
            SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).into();

        // Intersect the line with itself
        let intersections = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.count(), 0);

        // Intersect the line with itself
        let intersections = c.intersections_contour_par(&c, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.count(), 0);
    }
}

#[test]
fn test_centroid() {
    {
        // Rectangles
        let contour = Contour::rectangle([0.0, -2.0], [1.0, 0.0]);
        assert_abs_diff_eq!(&contour.centroid(), &[0.5, -1.0], epsilon = DEFAULT_EPSILON,);

        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        assert_abs_diff_eq!(&contour.centroid(), &[0.5, 1.0], epsilon = DEFAULT_EPSILON,);

        let contour = Contour::rectangle([1.0, -3.0], [2.0, 1.0]);
        assert_abs_diff_eq!(&contour.centroid(), &[1.5, -1.0], epsilon = DEFAULT_EPSILON,);

        // Circle
        let center = [0.0, 0.0];
        let contour = Contour::circle(center, 20.0);
        assert_abs_diff_eq!(&contour.centroid(), &center, epsilon = DEFAULT_EPSILON);

        let center = [1.0, -1.0];
        let contour = Contour::circle(center, 10.0);
        assert_abs_diff_eq!(&contour.centroid(), &center, epsilon = DEFAULT_EPSILON);

        // Quarter-circle => Compare with analytical solution
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet(
                [0.0, 2.0],
                [1.0, 2.0],
                [1.0, 1.0],
                1.0,
                DEFAULT_EPSILON,
                DEFAULT_MAX_ULPS,
            )
            .unwrap()
            .into(),
        );
        chain.extend_back([0.0, 1.0]);
        let contour = Contour::new(chain);
        let centroid = 4.0 / (3.0 * PI);
        assert_abs_diff_eq!(
            &contour.centroid(),
            &[centroid, centroid + 1.0],
            epsilon = DEFAULT_EPSILON
        );
    }

    // Regression test which was introduced after encountering a bug with a bread loaf magnet shape.
    {
        let width: f64 = 20e-3;
        let height = 10e-3;
        let outer_radius: f64 = 50e-3;
        let arc_segment_height =
            outer_radius - 0.5 * (4.0 * outer_radius.powi(2) - width.powi(2)).sqrt();

        let mut chain = SegmentChain::with_capacity(5);
        chain.push_back(
            LineSegment::new(
                [0.0, 0.0],
                [width / 2.0, 0.0],
                DEFAULT_EPSILON,
                DEFAULT_MAX_ULPS,
            )
            .unwrap()
            .into(),
        );
        chain.push_back(
            LineSegment::new(
                [width / 2.0, 0.0],
                [width / 2.0, height],
                DEFAULT_EPSILON,
                DEFAULT_MAX_ULPS,
            )
            .unwrap()
            .into(),
        );
        let arc = ArcSegment::from_start_middle_stop(
            [width / 2.0, height],
            [0.0, height + arc_segment_height],
            [-width / 2.0, height],
            0.0,
            0,
        )
        .unwrap()
        .into();
        chain.push_back(arc);
        chain.push_back(
            LineSegment::new(
                [-width / 2.0, height],
                [-width / 2.0, 0.0],
                DEFAULT_EPSILON,
                DEFAULT_MAX_ULPS,
            )
            .unwrap()
            .into(),
        );
        chain.push_back(
            LineSegment::new(
                [-width / 2.0, 0.0],
                [0.0, 0.0],
                DEFAULT_EPSILON,
                DEFAULT_MAX_ULPS,
            )
            .unwrap()
            .into(),
        );

        let contour = Contour::new(chain);
        approx::assert_abs_diff_eq!(&contour.centroid(), &[0.0, 5.341657e-3], epsilon = 1e-9);
    }
}

#[test]
fn test_contains() {
    {
        // large contains small
        let small = Contour::new(SegmentChain::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(SegmentChain::from_points(&[
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(
            large
                .intersections_composite(&small, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .count(),
            0
        );
        assert!(large.contains(&small, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        // large does not contain small, but they do not intersect and the
        // bounding box of small is within that of large
        let small = Contour::new(SegmentChain::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(SegmentChain::from_points(&[
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(
            large
                .intersections_composite(&small, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .count(),
            0
        );
        assert!(!large.contains(&small, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
}
