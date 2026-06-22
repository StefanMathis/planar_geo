use approx::{assert_abs_diff_eq, assert_relative_eq};
use planar_geo::prelude::*;
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI, SQRT_2, TAU};

#[test]
fn test_from_segments() {
    let mut polysegment = Polysegment::new();
    polysegment.push_back(
        ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, PI, FRAC_PI_2)
            .unwrap()
            .into(),
    );
    polysegment.push_back(LineSegment::new([0.0, -2.0], [-2.0, -2.0]).unwrap().into());
    polysegment.push_back(LineSegment::new([-2.0, -2.0], [-2.0, 0.0]).unwrap().into());
    let contour = Contour::new(polysegment);
    assert_eq!(contour.num_segments(), 5); // 2 glue segments
}

#[test]
fn test_from_points() {
    // Triangle
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let polygon: Contour = Polysegment::from_points(vertices).into();
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
fn test_intersection_cut() {
    {
        // Regression test (found in stem_slot crate)
        let arc = ArcSegment::from_center_radius_start_sweep_angle(
            [-0.0035575567453493533, 0.0017499999999999948],
            0.0010000000000000022,
            3.054659664751011,
            1.6577293156336776,
        )
        .expect("valid arc");
        let contour = Contour::from(arc.clone());
        let cut = Polysegment::from_points([[-10.0, 0.00075], [10.0, 0.00075]].as_slice());
        let separated_lines = contour.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(separated_lines.len(), 2);

        let cut_seg = &separated_lines[0][1];
        match cut_seg {
            Segment::LineSegment(_) => panic!("Must be an ArcSegment"),
            Segment::ArcSegment(cut_arc) => {
                approx::assert_abs_diff_eq!(arc.radius(), cut_arc.radius(), epsilon = 1e-3);
                approx::assert_abs_diff_eq!(arc.center(), cut_arc.center(), epsilon = 1e-3);
                approx::assert_abs_diff_eq!(
                    arc.start_angle(),
                    cut_arc.start_angle(),
                    epsilon = 1e-3
                );
                approx::assert_abs_diff_eq!(
                    arc.sweep_angle(),
                    cut_arc.sweep_angle(),
                    epsilon = 1e-3
                );
            }
        }
    }
    {
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
        let contour = Contour::new(Polysegment::from_points(vertices.as_slice()));

        let cut = Polysegment::from_points([[-10.0, 1.0], [10.0, 1.0]].as_slice());

        let separated_lines = contour.intersection_cut(&cut, 0.0, 0.0);
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
}

#[test]
fn test_points() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let polygon = Contour::new(Polysegment::from_points(vertices));
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
        assert_relative_eq!(contour.area(), area_circle);

        let bb = BoundingBox::from(&contour);
        assert_relative_eq!(bb.xmin(), -radius);
        assert_relative_eq!(bb.ymin(), -radius);
        assert_relative_eq!(bb.xmax(), radius);
        assert_relative_eq!(bb.ymax(), radius);
    }
    {
        let contour = planar_geo::prelude::Contour::circle([0.0, 0.0], 20.0);
        assert!(contour.covers_point([0.0, 0.0], 0.0, 0.0));
        assert!(!contour.covers_point([21.0, 0.0], 0.0, 0.0));

        // Bug found with the interactive canvas
        assert!(contour.covers_point(
            [0.055555555555542924, -0.9999923706054688],
            DEFAULT_EPSILON,
            DEFAULT_MAX_RELATIVE
        ));
        assert!(contour.covers_point(
            [12.055555555555543, -9.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_RELATIVE
        ));
        assert!(contour.covers_point(
            [8.055555555555543, -6.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_RELATIVE
        ));
        assert!(contour.covers_point(
            [14.055555555555543, -6.999992370605469],
            DEFAULT_EPSILON,
            DEFAULT_MAX_RELATIVE
        ));
    }
}

#[test]
fn test_area() {
    {
        let arc =
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, TAU).unwrap();
        let contour = Contour::new(arc.into());
        let area_circle = PI * 2.0_f64.powi(2);
        assert_relative_eq!(contour.area(), area_circle);
    }

    // Square with concave radius in one corner
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, PI, FRAC_PI_2)
                .unwrap()
                .into(),
        );
        polysegment.push_back(LineSegment::new([0.0, -2.0], [-2.0, -2.0]).unwrap().into());
        polysegment.push_back(LineSegment::new([-2.0, -2.0], [-2.0, 0.0]).unwrap().into());
        let contour = Contour::new(polysegment);
        assert_relative_eq!(contour.area(), 4.0 - 0.25 * PI * 2.0_f64.powi(2));
    }

    // Square with concave radius in another corner
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, -FRAC_PI_2)
                .unwrap()
                .into(),
        );
        polysegment.push_back(LineSegment::new([0.0, -2.0], [2.0, -2.0]).unwrap().into());
        polysegment.push_back(LineSegment::new([2.0, -2.0], [2.0, 0.0]).unwrap().into());
        let contour = Contour::new(polysegment);
        assert_relative_eq!(
            contour.area(),
            4.0 - 0.25 * PI * 2.0_f64.powi(2),
            epsilon = DEFAULT_EPSILON
        );
    }

    // Square with rounded edges
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([-1.0, 2.0], [-1.0, 0.0], [1.0, 0.0], 0.1)
                .unwrap()
                .into(),
        );
        polysegment.push_back(
            ArcSegment::fillet([-1.0, 0.0], [1.0, 0.0], [1.0, 2.0], 0.10)
                .unwrap()
                .into(),
        );
        polysegment.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 2.0], [-1.0, 2.0], 0.2)
                .unwrap()
                .into(),
        );
        polysegment.push_back(
            ArcSegment::fillet([1.0, 2.0], [-1.0, 2.0], [-1.0, 0.0], 0.2)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(polysegment);

        assert_relative_eq!(
            contour.area(),
            4.0 - 2.0 * 0.1f64.powi(2) * (1.0 - 0.25 * PI)
                - 2.0 * 0.2f64.powi(2) * (1.0 - 0.25 * PI),
            epsilon = DEFAULT_EPSILON
        );
    }

    // "Pie" section with a straight connection between the end points of the arc
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 1.0)
                .unwrap()
                .into(),
        );
        let contour = Contour::new(polysegment);
        assert_relative_eq!(contour.area(), 0.25 * PI - 0.5, epsilon = DEFAULT_EPSILON);
    }

    // "Pie" section with a square cut out at the origin
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([0.0, -1.0], [-1.0, -1.0], [-1.0, 0.0], 1.0)
                .unwrap()
                .into(),
        );
        polysegment.extend_back([-1.0, 0.0]);
        polysegment.extend_back([-0.5, 0.0]);
        polysegment.extend_back([-0.5, -0.5]);
        polysegment.extend_back([0.0, -0.5]);
        polysegment.extend_back([0.0, -1.0]);
        let contour = Contour::new(polysegment);
        assert_relative_eq!(contour.area(), 0.25 * PI - 0.25, epsilon = DEFAULT_EPSILON);
    }
}

#[test]
fn test_contains_point() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.contains_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[
            [0.0, 0.0],
            [0.0, 0.5],
            [0.1, 0.5],
            [0.1, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.contains_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[[0.0, 0.0], [0.0, 0.5], [0.1, 0.5]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(!contour.contains_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[[0.0, 0.0], [0.0, 0.5], [1.0, 0.0]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(!contour.contains_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.0, 0.5],
            [1.0, 1.0],
            [0.5, 1.0],
            [0.0, 1.0],
            [0.0, 0.5],
        ];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.contains_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.contains_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([-1.0, 0.0], [-2.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([-3.0, -2.0]);
        ps.extend_back([-1.0, 0.0]);
        ps.extend_back([0.0, -1.0]);
        ps.extend_back([0.0, 0.0]);
        ps.extend_back([0.5, 0.0]);
        let c = Contour::new(ps);

        assert!(c.contains_point([-1.5, 0.0], e, m));
        assert!(c.contains_point([0.0, 0.1], e, m));
        assert!(!c.contains_point([0.0, 0.0], e, m));
        assert!(c.contains_point([-1.5, 0.5], e, m));
        assert!(!c.contains_point([-1.0, 0.0], e, m));
        assert!(c.contains_point([-2.0, 0.1], e, m));
        assert!(!c.contains_point([-3.0, -1.0], e, m));
        assert!(c.contains_point([-2.0, -0.2], e, m));
        assert!(!c.contains_point([0.0, -1.0], e, m));
        assert!(!c.contains_point([-5.0, -1.0], e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([0.0, 1.0], [0.0, 0.0], -FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], -FRAC_PI_2)
                .unwrap()
                .into(),
        );
        let c = Contour::new(ps);
        assert!(c.contains_point([0.1, 0.0], e, m));
        assert!(c.contains_point([0.1, 1e-25], e, m));
        assert!(c.contains_point([0.1, -1e-25], e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([0.0, 1.0], [0.0, 0.0], -FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [2.0, 0.5], 0.1)
                .unwrap()
                .into(),
        );
        ps.extend_back([0.0, -1.0]);
        let c = Contour::new(ps);
        assert!(c.contains_point([0.1, 0.0], e, m));
        assert!(c.contains_point([0.1, 1e-25], e, m));
        assert!(c.contains_point([0.1, -1e-25], e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([0.0, 1.0], [0.0, 0.0], -FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [2.0, 0.0], FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.extend_back([0.0, -1.0]);
        let c = Contour::new(ps);
        assert!(c.contains_point([0.1, 0.0], e, m));
        assert!(c.contains_point([0.1, 1e-25], e, m));
        assert!(c.contains_point([0.1, -1e-25], e, m));
    }
}

#[test]
fn test_contains_point_core_contour() {
    // Test of a bug found in the stem_core crate

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;

    let mut ps = Polysegment::new();
    ps.push_back(
        LineSegment::new([0.0, 0.025], [0.15, 0.025])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.025], [0.15, 0.01774999999999992])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.01774999999999992], [0.149, 0.01774999999999992])
            .unwrap()
            .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.14899999999999997, 0.014749999999999961],
            0.0029999999999999593,
            1.5707963267948828,
            1.5707963267948963,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.146, 0.014750000000000001],
            [0.146, 0.0027499999999999985],
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.148, 0.002749999999999999],
            0.001999999999999999,
            3.141592653589793,
            1.5707963267948966,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new([0.148, 0.0007499999999999998], [0.15, 0.00075])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.00075], [0.15, 0.0])
            .unwrap()
            .into(),
    );
    ps.push_back(LineSegment::new([0.15, 0.0], [0.0, 0.0]).unwrap().into());
    ps.push_back(LineSegment::new([0.0, 0.0], [0.0, 0.00075]).unwrap().into());
    ps.push_back(
        LineSegment::new([0.0, 0.00075], [0.002, 0.0007499999999999998])
            .unwrap()
            .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.0019999999999999996, 0.002749999999999999],
            0.001999999999999999,
            4.71238898038469,
            1.5707963267948966,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.003999999999999998, 0.0027499999999999985],
            [0.004, 0.014750000000000001],
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.0010000000000000408, 0.014749999999999961],
            0.0029999999999999593,
            1.3877787807814645e-14,
            1.5707963267948963,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.0009999999999999996, 0.01774999999999992],
            [0.0, 0.01774999999999992],
        )
        .unwrap()
        .into(),
    );
    let c = Contour::from(ps);

    assert!(c.contains_point([0.0085, 0.0027499999999999994], e, m));
    assert!(c.contains_point([0.010499999999999999, 0.0007499999999999998], e, m));
    assert!(c.contains_point([0.013499999999999998, 0.01774999999999992], e, m));
    assert!(c.contains_point([0.0165, 0.014750000000000001], e, m));
    assert!(c.contains_point([0.014499999999999999, 0.0007499999999999998], e, m));
    assert!(c.contains_point([0.012499999999999999, 0.00075], e, m));
    assert!(c.contains_point([0.016499999999999997, 0.0027499999999999985], e, m));
    assert!(!c.contains_point([-2.0, 0.0027499999999999994], e, m));
    assert!(!c.contains_point([0.0, 0.00275], e, m));
}

#[test]
fn test_covers_point() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.covers_point([0.5, 0.5], 0.0, 0.0));
        assert!(contour.covers_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[
            [0.0, 0.0],
            [0.0, 0.5],
            [0.1, 0.5],
            [0.1, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.covers_point([0.5, 0.5], 0.0, 0.0));
        assert!(contour.covers_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[[0.0, 0.0], [0.0, 0.5], [0.1, 0.5]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(!contour.covers_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[[0.0, 0.0], [0.0, 0.5], [1.0, 0.0]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(!contour.covers_point([0.5, 0.5], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let vertices = &[
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.0, 0.5],
            [1.0, 1.0],
            [0.5, 1.0],
            [0.0, 1.0],
            [0.0, 0.5],
        ];
        let contour = Contour::new(Polysegment::from_points(vertices));
        assert!(contour.covers_point([0.5, 0.5], 0.0, 0.0));
        assert!(contour.covers_point([0.5, 1.0], 0.0, 0.0));
        assert!(!contour.covers_point([0.5, 1.5], 0.0, 0.0));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([-1.0, 0.0], [-2.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([-3.0, -2.0]);
        ps.extend_back([-1.0, 0.0]);
        ps.extend_back([0.0, -1.0]);
        ps.extend_back([0.0, 0.0]);
        ps.extend_back([0.5, 0.0]);
        let c = Contour::new(ps);

        assert!(c.covers_point([0.0, 0.1], e, m));
        assert!(c.covers_point([0.0, 0.0], e, m));
        assert!(c.covers_point([-1.5, 0.0], e, m));
        assert!(c.covers_point([-1.0, 0.0], e, m));
        assert!(c.covers_point([-2.0, 0.1], e, m));
        assert!(c.covers_point([-3.0, -1.0], e, m));
        assert!(c.covers_point([-2.0, -0.2], e, m));
        assert!(c.covers_point([0.0, -1.0], e, m));
        assert!(!c.covers_point([-5.0, -1.0], e, m));
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(c.covers_point([0.0, 0.0], e, m));
        assert!(!c.covers_point([0.9, 0.9], e, m));
        assert!(!c.covers_point([0.999999999, 0.999999999], e, m));
        assert!(!c.covers_point([1.0, 1.0], e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([-1.0, 2.0]);
        ps.extend_back([1.0, 2.0]);
        let c = Contour::new(ps);
        assert!(c.covers_point([-1.0, 1.0], e, m));
        assert!(c.covers_point([1.0, 1.0], e, m));
        assert!(c.covers_point([0.5, 1.0], e, m));
    }
}

#[test]
fn test_self_intersection() {
    {
        // Regression test from stem_magnet crate
        let arc1 = ArcSegment::from_center_radius_start_sweep_angle(
            [0.0, -0.06],
            0.06,
            1.3089969389957472,
            0.5235987755982988,
        )
        .unwrap();
        let arc2 = ArcSegment::from_center_radius_start_sweep_angle(
            [-3.469446951953614e-18, -0.027712440737676466],
            0.03,
            2.1148844334769814,
            -1.0881762133641695,
        )
        .unwrap();
        let c: Contour = Polysegment::from_iter([arc1.into(), arc2.into()].into_iter()).into();
        let intersections = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.count(), 0);
    }
    // Open polysegment
    {
        let c: Contour =
            Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).into();

        // Intersect the line with itself
        let intersections = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.count(), 0);

        // Intersect the line with itself
        let intersections = c.intersections_contour_par(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.count(), 0);

        // Intersect the line with itself
        let intersections = c.intersections(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.len(), 0);

        let intersections = c.intersections_par(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.len(), 0);
    }

    // Closed polysegment
    {
        let c: Contour =
            Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).into();

        // Intersect the line with itself
        let intersections = c.intersections_contour(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
        assert_eq!(intersections.count(), 0);

        // Intersect the line with itself
        let intersections = c.intersections_contour_par(&c, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE);
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
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([0.0, 2.0], [1.0, 2.0], [1.0, 1.0], 1.0)
                .unwrap()
                .into(),
        );
        polysegment.extend_back([0.0, 1.0]);
        let contour = Contour::new(polysegment);
        let centroid = 4.0 / (3.0 * PI);
        assert_abs_diff_eq!(
            &contour.centroid(),
            &[centroid, centroid + 1.0],
            epsilon = DEFAULT_EPSILON
        );
    }

    // Regression test which was introduced after encountering a bug with a bread
    // loaf magnet shape.
    {
        let width: f64 = 20e-3;
        let height = 10e-3;
        let outer_radius: f64 = 50e-3;
        let arc_segment_height =
            outer_radius - 0.5 * (4.0 * outer_radius.powi(2) - width.powi(2)).sqrt();

        let mut polysegment = Polysegment::with_capacity(5);
        polysegment.push_back(
            LineSegment::new([0.0, 0.0], [width / 2.0, 0.0])
                .unwrap()
                .into(),
        );
        polysegment.push_back(
            LineSegment::new([width / 2.0, 0.0], [width / 2.0, height])
                .unwrap()
                .into(),
        );
        let arc = ArcSegment::from_start_middle_stop(
            [width / 2.0, height],
            [0.0, height + arc_segment_height],
            [-width / 2.0, height],
        )
        .unwrap()
        .into();
        polysegment.push_back(arc);
        polysegment.push_back(
            LineSegment::new([-width / 2.0, height], [-width / 2.0, 0.0])
                .unwrap()
                .into(),
        );
        polysegment.push_back(
            LineSegment::new([-width / 2.0, 0.0], [0.0, 0.0])
                .unwrap()
                .into(),
        );

        let contour = Contour::new(polysegment);
        approx::assert_abs_diff_eq!(&contour.centroid(), &[0.0, 5.341657e-3], epsilon = 1e-9);
    }
}

#[test]
fn test_intersects() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    let c1: Contour =
        Polysegment::from(ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI).unwrap())
            .into();
    let c2: Contour =
        Polysegment::from(ArcSegment::from_start_center_angle([1.0, 2.0], [1.0, 1.0], PI).unwrap())
            .into();
    let c3: Contour =
        Polysegment::from(ArcSegment::from_start_center_angle([2.0, 2.0], [2.0, 1.0], PI).unwrap())
            .into();
    let c4: Contour = Polysegment::from(
        ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], -PI).unwrap(),
    )
    .into();

    assert_eq!(c1.intersection_cut(c2.polysegment(), e, m).len(), 4);
    assert_eq!(c1.intersection_cut(c3.polysegment(), e, m).len(), 1);
    assert_eq!(c1.intersection_cut(c4.polysegment(), e, m).len(), 2);
}

#[test]
fn test_covers_line_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ]));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 0.0], [0.0, 1.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 1.0], [0.0, 0.0]).unwrap(), e, m));
    }
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        assert!(c.covers_segment(&LineSegment::new([0.1, 0.1], [0.9, 0.9]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([-0.1, 0.5], [0.9, 0.5]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([-0.1, 0.5], [0.9, 1.5]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.1, 0.5], [0.9, 1.5]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 0.0], [2.0, 0.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([-1.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([0.0, 1.0], [1.0, 1.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 0.5], [1.0, 1.5]).unwrap(), e, m));
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(c.covers_segment(&LineSegment::new([-1.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([-0.5, 0.1], [0.5, 0.1]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([-0.5, -0.1], [0.5, -0.1]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 0.1], [0.0, 1.1]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([0.0, 0.1], [0.0, 1.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap(), e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([-1.0, 2.0]);
        ps.extend_back([1.0, 2.0]);
        let c = Contour::new(ps);
        assert!(c.covers_segment(&LineSegment::new([-1.0, 1.0], [1.0, 1.0]).unwrap(), e, m));
        assert!(c.covers_segment(&LineSegment::new([-0.1, 1.0], [0.1, 1.0]).unwrap(), e, m));
        assert!(!c.covers_segment(&LineSegment::new([-0.1, 0.9], [0.1, 0.9]).unwrap(), e, m));
    }
}

#[test]
fn test_contains_line_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        assert!(c.contains_segment(&LineSegment::new([0.1, 0.1], [0.9, 0.9]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-0.1, 0.5], [0.9, 0.5]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-0.1, 0.5], [0.9, 1.5]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.1, 0.5], [0.9, 1.5]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.0], [2.0, 0.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-1.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 1.0], [1.0, 1.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.5], [1.0, 1.5]).unwrap(), e, m));
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(!c.contains_segment(&LineSegment::new([-1.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(c.contains_segment(&LineSegment::new([-0.5, 0.1], [0.5, 0.1]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-0.5, -0.1], [0.5, -0.1]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.1], [0.0, 1.1]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.1], [0.0, 1.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap(), e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([-1.0, 2.0]);
        ps.extend_back([1.0, 2.0]);
        let c = Contour::new(ps);
        assert!(!c.contains_segment(&LineSegment::new([-1.0, 1.0], [1.0, 1.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-0.1, 1.0], [0.1, 1.0]).unwrap(), e, m));
        assert!(!c.contains_segment(&LineSegment::new([-0.1, 0.9], [0.1, 0.9]).unwrap(), e, m));
    }
}

#[test]
fn test_covers_arc_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        assert!(c.covers_segment(&ArcSegment::circle([0.5, 0.5], 0.5).unwrap(), e, m));
        assert!(c.covers_segment(&ArcSegment::circle([0.5, 0.5], 0.4).unwrap(), e, m));
        assert!(!c.covers_segment(&ArcSegment::circle([0.5, 0.5], 0.6).unwrap(), e, m));
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.0, 0.0],
                    1.0,
                    FRAC_PI_2,
                    FRAC_PI_2,
                )
                .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, PI, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.5, 0.5],
                    0.5,
                    1.5 * FRAC_PI_2,
                    PI,
                )
                .unwrap(),
                e,
                m
            )
        );
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.5, 0.5], 0.5, FRAC_PI_2, PI,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(c.covers_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.5, 0.5], 0.5, 0.0, PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.5, 0.5], 0.5, PI, PI).unwrap(),
            e,
            m
        ));
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.covers_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2).unwrap(),
            e,
            m
        ));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, FRAC_PI_2, FRAC_PI_2)
                .unwrap()
                .into(),
        );
        let c = Contour::new(ps);
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.covers_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2).unwrap(),
            e,
            m
        ));
        assert!(!c.covers_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2 + 0.1).unwrap(),
            e,
            m
        ));
        assert!(!c.covers_segment(
            &ArcSegment::from_start_center_angle([1.0, 2.0], [1.0, 1.0], FRAC_PI_2 + 0.1).unwrap(),
            e,
            m
        ));
    }
    {
        let c: Contour = Polysegment::from(ArcSegment::circle([0.0, 0.0], 2.0).unwrap()).into();
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_start_center_angle([2.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(!c.covers_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.1, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(c.covers_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.9, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, PI - 0.1, PI,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, FRAC_PI_2, PI,)
                    .unwrap(),
                e,
                m
            )
        );
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, FRAC_PI_2, PI)
                .unwrap(),
        )
        .into();
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.0, 0.0],
                    2.0,
                    1.25 * FRAC_PI_2,
                    FRAC_PI_2,
                )
                .unwrap(),
                e,
                m
            )
        );
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, -FRAC_PI_2, PI)
                .unwrap(),
        )
        .into();
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.0, 0.0],
                    2.0,
                    -0.75 * FRAC_PI_2,
                    FRAC_PI_2,
                )
                .unwrap(),
                e,
                m
            )
        );
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(
            c.covers_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, PI - 0.2,)
                    .unwrap(),
                e,
                m
            )
        );
    }
}

#[test]
fn test_contains_arc_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        assert!(!c.contains_segment(&ArcSegment::circle([0.5, 0.5], 0.5).unwrap(), e, m));
        assert!(c.contains_segment(&ArcSegment::circle([0.5, 0.5], 0.4).unwrap(), e, m));
        assert!(!c.contains_segment(&ArcSegment::circle([0.5, 0.5], 0.6).unwrap(), e, m));
        assert!(
            !c.contains_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.contains_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.contains_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.0, 0.0],
                    1.0,
                    FRAC_PI_2,
                    FRAC_PI_2,
                )
                .unwrap(),
                e,
                m
            )
        );
        assert!(
            !c.contains_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, PI, FRAC_PI_2,)
                    .unwrap(),
                e,
                m
            )
        );
    }
    {
        let c: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2).unwrap(),
            e,
            m
        ));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, FRAC_PI_2, FRAC_PI_2)
                .unwrap()
                .into(),
        );
        let c = Contour::new(ps);
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([0.0, 1.0], [1.0, 1.0], FRAC_PI_2 + 0.1).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([1.0, 2.0], [1.0, 1.0], FRAC_PI_2 + 0.1).unwrap(),
            e,
            m
        ));
    }
    {
        let c: Contour = Polysegment::from(ArcSegment::circle([0.0, 0.0], 2.0).unwrap()).into();
        assert!(c.contains_segment(
            &ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_start_center_angle([2.0, 0.0], [0.0, 0.0], PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(!c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.1, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.9, 0.5, PI).unwrap(),
            e,
            m
        ));
        assert!(
            !c.contains_segment(
                &ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, PI - 0.1, PI,)
                    .unwrap(),
                e,
                m
            )
        );
    }
}

#[test]
fn test_covers_contour() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;

    // A countour shares a side with another one
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let c2 = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [1.5, 2.0],
            [1.5, 1.0],
        ]));
        assert!(c1.covers_contour(&c2, e, m))
    }
    // A countour covers itself
    {
        let c = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        assert!(c.covers_contour(&c, e, m))
    }
    {
        // large contains small
        let small = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(large.intersections_composite(&small, e, m).count(), 0);
        assert!(large.covers_contour(&small, e, m));
    }
    {
        // large does not contain small, but they do not intersect and the
        // bounding box of small is within that of large
        let small = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(Polysegment::from_points(&[
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(large.intersections_composite(&small, e, m).count(), 0);
        assert!(!large.covers_contour(&small, e, m));
    }
    {
        let c1: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        let c2: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, PI - 0.2)
                .unwrap(),
        )
        .into();
        assert!(c1.covers_contour(&c2, e, m));
    }
    {
        let c1: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        let c2: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle(
                [0.0, 0.0],
                1.0,
                PI - 0.1,
                -(PI - 0.2),
            )
            .unwrap(),
        )
        .into();
        assert!(c1.covers_contour(&c2, e, m));
    }
    {
        let c1: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([-1.0, 0.0], [0.0, 0.0], -PI).unwrap(),
        )
        .into();
        let c2: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, PI - 0.2)
                .unwrap(),
        )
        .into();
        assert!(c1.covers_contour(&c2, e, m));
    }
    {
        let c1: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([-1.0, 0.0], [0.0, 0.0], -PI).unwrap(),
        )
        .into();
        let c2: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle(
                [0.0, 0.0],
                1.0,
                PI - 0.1,
                -(PI - 0.2),
            )
            .unwrap(),
        )
        .into();
        assert!(c1.covers_contour(&c2, e, m));
    }
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [-1.0, -1.0],
            [1.0, -1.0],
            [1.0, 1.0],
            [-1.0, 1.0],
        ]));
        let c2 = Contour::from(ArcSegment::circle([0.0, 0.0], 1.0).unwrap());
        assert!(c1.covers_contour(&c2, e, m));
        assert!(!c2.covers_contour(&c1, e, m));
    }
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [-1.0, -1.0],
            [1.0, -1.0],
            [1.0, 1.0],
            [-1.0, 1.0],
        ]));
        let c2 = Contour::from(ArcSegment::circle([0.0, 0.0], SQRT_2).unwrap());
        assert!(!c1.covers_contour(&c2, e, m));
        assert!(c2.covers_contour(&c1, e, m));
    }
}

#[test]
fn test_contains_contour() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;

    {
        let c1 = Contour::rectangle([0.0, 0.0], [0.15, 0.25]);
        let c2 = Contour::rectangle(
            [0.008499999999999999, 0.0007499999999999998],
            [0.0165, 0.01774999999999992],
        );
        assert!(c1.contains_contour(&c2, e, m));
    }
    {
        let pts = &[
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
            [0.5, 0.0],
        ];
        let radii = &[0.1, 0.1, 0.1, 0.1];
        let c1: Contour = Polysegment::from_fillet_chain(pts, radii).into();

        let pts = &[
            [0.5, 0.2],
            [0.2, 0.0],
            [0.2, 0.8],
            [0.8, 0.8],
            [0.8, 0.2],
            [0.5, 0.2],
        ];
        let radii = &[0.1, 0.1, 0.1, 0.1];
        let c2 = Polysegment::from_fillet_chain(pts, radii).into();
        assert!(c1.contains_contour(&c2, e, m));
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(LineSegment::new([0.0, 0.0], [0.0, 0.3]).unwrap().into());
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.5], 0.2, 1.5 * PI, PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([0.0, 1.0]);
        ps.extend_back([1.0, 1.0]);
        ps.extend_back([1.0, 0.0]);
        let c1 = Contour::from(ps);
        let c2 = Contour::rectangle([0.5, 0.2], [0.8, 0.8]);
        assert!(c1.contains_contour(&c2, e, m));
    }

    // A countour shares a side with another one
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let c2 = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [1.5, 2.0],
            [1.5, 1.0],
        ]));
        assert!(!c1.contains_contour(&c2, e, m));
    }
    // A countour does not contain itself
    {
        let c = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        assert!(!c.contains_contour(&c, e, m));
    }
    {
        // large contains small
        let small = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(large.intersections_composite(&small, e, m).count(), 0);
        assert!(large.contains_contour(&small, e, m));
    }
    {
        // large does not contain small, but they do not intersect and the
        // bounding box of small is within that of large
        let small = Contour::new(Polysegment::from_points(&[
            [1.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [2.0, 1.0],
        ]));
        let large = Contour::new(Polysegment::from_points(&[
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]));

        assert!(BoundingBox::from(&large).contains(&BoundingBox::from(&small)));
        assert_eq!(large.intersections_composite(&small, e, m).count(), 0);
        assert!(!large.contains_contour(&small, e, m));
    }
    {
        let c1: Contour = Polysegment::from(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI).unwrap(),
        )
        .into();
        let c2: Contour = Polysegment::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.1, PI - 0.2)
                .unwrap(),
        )
        .into();
        assert!(!c1.contains_contour(&c2, e, m));
    }
}

#[test]
fn test_contains_segments() {
    // Test of a bug found in the stem_core crate

    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;

    let mut ps = Polysegment::new();
    ps.push_back(
        LineSegment::new([0.0, 0.025], [0.15, 0.025])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.025], [0.15, 0.01774999999999992])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.01774999999999992], [0.149, 0.01774999999999992])
            .unwrap()
            .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.14899999999999997, 0.014749999999999961],
            0.0029999999999999593,
            1.5707963267948828,
            1.5707963267948963,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.146, 0.014750000000000001],
            [0.146, 0.0027499999999999985],
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.148, 0.002749999999999999],
            0.001999999999999999,
            3.141592653589793,
            1.5707963267948966,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new([0.148, 0.0007499999999999998], [0.15, 0.00075])
            .unwrap()
            .into(),
    );
    ps.push_back(
        LineSegment::new([0.15, 0.00075], [0.15, 0.0])
            .unwrap()
            .into(),
    );
    ps.push_back(LineSegment::new([0.15, 0.0], [0.0, 0.0]).unwrap().into());
    ps.push_back(LineSegment::new([0.0, 0.0], [0.0, 0.00075]).unwrap().into());
    ps.push_back(
        LineSegment::new([0.0, 0.00075], [0.002, 0.0007499999999999998])
            .unwrap()
            .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.0019999999999999996, 0.002749999999999999],
            0.001999999999999999,
            4.71238898038469,
            1.5707963267948966,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.003999999999999998, 0.0027499999999999985],
            [0.004, 0.014750000000000001],
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        ArcSegment::from_center_radius_start_sweep_angle(
            [0.0010000000000000408, 0.014749999999999961],
            0.0029999999999999593,
            1.3877787807814645e-14,
            1.5707963267948963,
        )
        .unwrap()
        .into(),
    );
    ps.push_back(
        LineSegment::new(
            [0.0009999999999999996, 0.01774999999999992],
            [0.0, 0.01774999999999992],
        )
        .unwrap()
        .into(),
    );
    let c = Contour::from(ps);

    assert!(
        c.contains_segment(
            &LineSegment::new(
                [0.012499999999999999, 0.00075],
                [0.014499999999999999, 0.0007499999999999998],
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle(
                [0.014499999999999999, 0.002749999999999999],
                0.001999999999999999,
                4.71238898038469,
                1.5707963267948966,
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &LineSegment::new(
                [0.016499999999999997, 0.0027499999999999985],
                [0.0165, 0.014750000000000001],
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle(
                [0.01350000000000004, 0.014749999999999961],
                0.0029999999999999593,
                1.3877787807814645e-14,
                1.5707963267948963,
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &LineSegment::new(
                [0.013499999999999998, 0.01774999999999992],
                [0.0115, 0.01774999999999992],
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle(
                [0.011499999999999958, 0.014749999999999961],
                0.0029999999999999593,
                1.5707963267948828,
                1.5707963267948963,
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &LineSegment::new(
                [0.008499999999999999, 0.014750000000000001],
                [0.0085, 0.0027499999999999985],
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &ArcSegment::from_center_radius_start_sweep_angle(
                [0.010499999999999999, 0.002749999999999999],
                0.001999999999999999,
                3.141592653589793,
                1.5707963267948966,
            )
            .unwrap(),
            e,
            m
        )
    );
    assert!(
        c.contains_segment(
            &LineSegment::new(
                [0.010499999999999999, 0.0007499999999999998],
                [0.012499999999999999, 0.00075],
            )
            .unwrap(),
            e,
            m
        )
    );
}

#[test]
fn test_overlaps_line_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        assert!(c.overlaps_segment(&LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap(), e, m));
        assert!(c.overlaps_segment(&LineSegment::new([-0.5, 0.5], [0.5, 0.5]).unwrap(), e, m));
        assert!(c.overlaps_segment(&LineSegment::new([-1.0, 2.0], [1.0, 0.0]).unwrap(), e, m));
        assert!(c.overlaps_segment(&LineSegment::new([-1.0, 2.0], [0.5, 0.5]).unwrap(), e, m));
        assert!(!c.overlaps_segment(&LineSegment::new([0.5, -0.5], [1.5, 0.5]).unwrap(), e, m));
        assert!(!c.overlaps_segment(&LineSegment::new([0.0, 0.0], [0.0, 1.0]).unwrap(), e, m));
        assert!(!c.overlaps_segment(&LineSegment::new([0.0, -1.0], [0.0, 0.0]).unwrap(), e, m));
        assert!(!c.overlaps_segment(&LineSegment::new([0.0, -2.0], [0.0, -1.0]).unwrap(), e, m));
    }
}

#[test]
fn test_overlaps_arc_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        for (center, start_angle) in [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
            .into_iter()
            .zip([0.0, FRAC_PI_2, 2.0 * FRAC_PI_2, 3.0 * FRAC_PI_2].into_iter())
        {
            assert!(
                c.overlaps_segment(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        center,
                        1.0,
                        start_angle,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                    e,
                    m
                )
            );
            assert!(
                c.overlaps_segment(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        center,
                        1.1,
                        start_angle,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                    e,
                    m
                )
            );
            assert!(
                c.overlaps_segment(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        center,
                        1.1,
                        start_angle + 0.1,
                        FRAC_PI_2 - 0.2,
                    )
                    .unwrap(),
                    e,
                    m
                )
            );
            assert!(
                c.overlaps_segment(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        center,
                        1.1,
                        start_angle - 0.1,
                        FRAC_PI_2 + 0.1,
                    )
                    .unwrap(),
                    e,
                    m
                )
            );
            assert!(
                !c.overlaps_segment(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        center,
                        SQRT_2,
                        start_angle,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                    e,
                    m
                )
            );
        }
    }
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        assert!(
            !c.overlaps_segment(
                &ArcSegment::from_center_radius_start_sweep_angle(
                    [0.1, 0.1],
                    SQRT_2,
                    0.0,
                    FRAC_PI_2,
                )
                .unwrap(),
                e,
                m
            )
        );
    }
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let c2 = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));
        assert!(c1.overlaps_contour(&c2, e, m));
        assert!(c2.overlaps_contour(&c1, e, m));
    }
}

#[test]
fn test_overlaps_contour() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let c2 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [0.5, 1.0],
            [0.5, 0.0],
        ]));
        let c3 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [-1.0, 1.0],
            [-1.0, 0.0],
        ]));
        let c4 = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));

        assert!(c1.overlaps_contour(&c1, e, m));

        assert!(c1.overlaps_contour(&c2, e, m));
        assert!(c2.overlaps_contour(&c1, e, m));

        assert!(!c1.overlaps_contour(&c3, e, m));
        assert!(!c3.overlaps_contour(&c1, e, m));

        assert!(c1.overlaps_contour(&c4, e, m));
        assert!(c4.overlaps_contour(&c1, e, m));
    }
    {
        let c1 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]));
        let c2 = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
        ]));
        assert!(!c1.overlaps_contour(&c2, e, m));
        assert!(!c2.overlaps_contour(&c1, e, m));
    }
    {
        let c1 = Contour::from(ArcSegment::circle([0.0, 0.0], 1.0).unwrap());
        assert!(c1.overlaps_contour(&c1, e, m));

        let c2 = Contour::from(ArcSegment::circle([2.0, 0.0], 1.0).unwrap());
        assert!(!c1.overlaps_contour(&c2, e, m));
        assert!(!c2.overlaps_contour(&c1, e, m));

        let c3 = Contour::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2)
                .unwrap(),
        );
        assert!(c1.overlaps_contour(&c3, e, m));
        assert!(c3.overlaps_contour(&c1, e, m));

        let c4 = Contour::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.2, 0.0], 1.0, 0.0, FRAC_PI_2)
                .unwrap(),
        );
        assert!(c1.overlaps_contour(&c4, e, m));
        assert!(c4.overlaps_contour(&c1, e, m));

        let c5 = Contour::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.2, 0.0], 1.0, FRAC_PI_2, FRAC_PI_2)
                .unwrap(),
        );
        assert!(c1.overlaps_contour(&c5, e, m));
        assert!(c5.overlaps_contour(&c1, e, m));
    }
    {
        let c1 = Contour::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.5, PI).unwrap(),
        );
        let c2 = Contour::from(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 1.0, 0.5 + PI, PI)
                .unwrap(),
        );
        assert!(!c1.overlaps_contour(&c2, e, m));
        assert!(!c2.overlaps_contour(&c1, e, m));
    }
}
