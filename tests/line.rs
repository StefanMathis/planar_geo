use std::f64::consts::{FRAC_PI_4, PI};

use approx;
use planar_geo::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS, prelude::*};

#[test]
fn test_contains_point() {
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        assert!(line.contains_point([1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(!line.contains_point([1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(line.contains_point([1.0, 0.0], 1.0, 0));
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        assert!(line.contains_point([1.0, -3.0f64.sqrt()], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let line = Line::from_point_angle([0.0, -2.0], PI / 3.0);
        assert!(line.contains_point(
            [1.0, 3.0f64.sqrt() - 2.0],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
    }
    {
        // Vertical line
        let line =
            LineSegment::new([0.5, -0.5], [0.5, 0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert!(line.contains_point([0.5, 0.0], 0.0, 0));
    }
    {
        // Vertical line
        let line =
            LineSegment::new([0.5, -0.5], [0.5, 0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert!(!line.contains_point([0.5, 1.0], 0.0, 0));
    }
    {
        let line =
            LineSegment::new([1.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert!(line.contains_point([0.5, 0.5], 0.0, 0));
    }
    {
        // This point is clearly not directly on the line
        let line =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let point = [0.5, 1.0];
        assert!(!line.contains_point(point, 0.0, 0));

        // But it is within this tolerance
        assert!(line.contains_point(point, 1.0, 1));

        // As well as inside this tolerance
        assert!(line.contains_point(point, 1.00001, 0));

        // But not outside this tolerance
        assert!(!line.contains_point(point, 0.99999, 0));

        // Repeat the test with a point located on a straight extension of the original line
        let point = [2.0, 0.0];
        assert!(!line.contains_point(point, 0.0, 0));
        assert!(line.contains_point(point, 1.0, 0));
        assert!(line.contains_point(point, 1.00001, 0));
        assert!(!line.contains_point(point, 0.99999, 0));

        // This point is outside of the circular radius around the endpoint [1.0, 0.0]
        assert!(!line.contains_point([1.1, -0.1], 0.11, 0));
        assert!(line.contains_point([0.5, -0.1], 0.11, 0));
        assert!(!line.contains_point([1.1, -0.1], 0.1, 0));
        assert!(!line.contains_point([1.1, -0.1], 0.09, 0));
    }
    {
        let line =
            LineSegment::new([0.0, 0.0], [1000.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

        // This point is outside of the circular radius around the endpoint [1000.0, 0.0]
        assert!(!line.contains_point([1100.0, -100.0], 110.0, 0));
        assert!(!line.contains_point([1100.0, -100.0], 100.0, 0));
    }
}

#[test]
fn test_identical() {
    {
        let line_segment =
            LineSegment::new([0.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([-2.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert!(!first.identical(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS))
    }
    {
        let line_segment =
            LineSegment::new([0.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([1.0, 1.0], [3.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(first.identical(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS))
    }
}

#[test]
fn test_intersection() {
    {
        let line_segment =
            LineSegment::new([0.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([0.0, 1.0], [2.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::Zero
        );
    }
    {
        let line_segment =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::One([0.0, 0.0])
        );
    }
    {
        let line_segment =
            LineSegment::new([3.0, 0.0], [3.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([0.0, 2.0], [1.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::One([3.0, 2.0])
        );
    }
    {
        let line_segment =
            LineSegment::new([0.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let first = Line::from(&line_segment);

        let line_segment =
            LineSegment::new([-2.0, 0.0], [2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::One([2.0, 2.0])
        );
    }
    {
        let first = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        let second = Line::from_point_angle([0.0, -2.0], PI / 3.0);

        approx::assert_abs_diff_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::One([1.0 / 3.0f64.sqrt(), -1.0])
        );
    }
    {
        let line_1 = Line::from_point_angle([0.11, 0.11], 0.0);
        let line_2 = Line::from_point_angle([0.11, 0.11], std::f64::consts::FRAC_PI_2);
        assert_eq!(
            line_1
                .intersections_line(&line_2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .len(),
            1
        );
        assert_eq!(
            line_2
                .intersections_line(&line_1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .len(),
            1
        );
    }
}

#[test]
fn test_self_intersection() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );
}

#[test]
fn test_intersection_with_contained_line_segment() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&ls, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );
    approx::assert_abs_diff_eq!(
        ls.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );
}

#[test]
fn transformation() {
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.contains_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        line.translate([2.0, -2.0]);
        assert!(line.contains_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.contains_point([0.0, 4.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        line.rotate([2.0, 2.0], 0.5 * PI);
        assert!(line.contains_point([0.0, 4.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let mut line = Line::from_point_angle([1.0, 1.0], -0.25 * PI);
        assert!(!line.contains_point([2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        line.scale(2.0);
        assert!(line.contains_point([2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.contains_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        line.line_reflection([-1.0, 0.0], [1.0, 0.0]);
        assert!(line.contains_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.contains_point([0.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        line.point_reflection([1.0, 0.0]);
        assert!(line.contains_point([0.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
}
