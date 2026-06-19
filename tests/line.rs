use std::f64::consts::{FRAC_PI_4, PI};

use approx;
use planar_geo::{DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE, prelude::*};

#[test]
fn test_convert_to_geo() {
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        let geo = Geometry::from(line);
        let _ = GeometryRef::from(&geo);
        let geo_cow = GeometryCow::from(geo);
        let _ = GeometryRef::from(&geo_cow);
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        let _ = GeometryRef::from(&line);
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        let geo_cow = GeometryCow::from(line);
        let _ = GeometryRef::from(&geo_cow);
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        let geo_cow = GeometryCow::from(&line);
        let _ = GeometryRef::from(&geo_cow);
    }
}

#[test]
fn test_covers_point() {
    {
        let line = Line::from_point_angle([0.0, 1.0], 0.0);
        assert!(line.covers_point([1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(line.covers_point([2.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(line.covers_point([100000.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(!line.covers_point([100000.0, 0.9], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        assert!(line.covers_point([1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(!line.covers_point([1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(line.covers_point([1.0, 0.0], 1.0, 0.0));
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        assert!(line.covers_point([1.0, -3.0f64.sqrt()], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let line = Line::from_point_angle([0.0, -2.0], PI / 3.0);
        assert!(line.covers_point(
            [1.0, 3.0f64.sqrt() - 2.0],
            DEFAULT_EPSILON,
            DEFAULT_MAX_RELATIVE
        ));
    }
    {
        // Vertical line
        let line = LineSegment::new([0.5, -0.5], [0.5, 0.5]).unwrap();
        assert!(line.covers_point([0.5, 0.0], 0.0, 0.0));
    }
    {
        // Vertical line
        let line = LineSegment::new([0.5, -0.5], [0.5, 0.5]).unwrap();
        assert!(!line.covers_point([0.5, 1.0], 0.0, 0.0));
    }
    {
        let line = LineSegment::new([1.0, 0.0], [0.0, 1.0]).unwrap();
        assert!(line.covers_point([0.5, 0.5], 0.0, 0.0));
    }
    {
        // This point is clearly not directly on the line
        let line = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
        let point = [0.5, 1.0];
        assert!(!line.covers_point(point, 0.0, 0.0));

        // But it is within this tolerance
        assert!(line.covers_point(point, 1.0, 0.0));

        // As well as inside this tolerance
        assert!(line.covers_point(point, 1.00001, 0.0));

        // But not outside this tolerance
        assert!(!line.covers_point(point, 0.99999, 0.0));

        // Repeat the test with a point located on a straight extension of the original
        // line
        let point = [2.0, 0.0];
        assert!(!line.covers_point(point, 0.0, 0.0));
        assert!(line.covers_point(point, 1.0, 0.0));
        assert!(line.covers_point(point, 1.00001, 0.0));
        assert!(!line.covers_point(point, 0.99999, 0.0));

        // This point is outside of the circular radius around the endpoint [1.0, 0.0]
        assert!(!line.covers_point([1.1, -0.1], 0.11, 0.0));
        assert!(line.covers_point([0.5, -0.1], 0.11, 0.0));
        assert!(!line.covers_point([1.1, -0.1], 0.1, 0.0));
        assert!(!line.covers_point([1.1, -0.1], 0.09, 0.0));
    }
    {
        let line = LineSegment::new([0.0, 0.0], [1000.0, 0.0]).unwrap();

        // This point is outside of the circular radius around the endpoint [1000.0,
        // 0.0]
        assert!(!line.covers_point([1100.0, -100.0], 110.0, 0.0));
        assert!(!line.covers_point([1100.0, -100.0], 100.0, 0.0));
    }
}

#[test]
fn test_identical() {
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([-2.0, 0.0], [2.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(!first.identical(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE))
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([1.0, 1.0], [3.0, 3.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert!(first.identical(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE))
    }
}

#[test]
fn test_intersection_line_line() {
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 1.0], [2.0, 3.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::Zero
        );
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [0.0, 1.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::One([0.0, 0.0])
        );
    }
    {
        let line_segment = LineSegment::new([3.0, 0.0], [3.0, 1.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 2.0], [1.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::One([3.0, 2.0])
        );
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([-2.0, 0.0], [2.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::One([2.0, 2.0])
        );
    }
    {
        let first = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        let second = Line::from_point_angle([0.0, -2.0], PI / 3.0);

        approx::assert_abs_diff_eq!(
            first.intersections_primitive(&second, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::One([1.0 / 3.0f64.sqrt(), -1.0])
        );
    }
    {
        let line_1 = Line::from_point_angle([0.11, 0.11], 0.0);
        let line_2 = Line::from_point_angle([0.11, 0.11], std::f64::consts::FRAC_PI_2);
        assert_eq!(
            line_1
                .intersections_line(&line_2, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
                .len(),
            1
        );
        assert_eq!(
            line_2
                .intersections_line(&line_1, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE)
                .len(),
            1
        );
    }
}

#[test]
fn test_self_intersection() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
        PrimitiveIntersections::Zero
    );
}

#[test]
fn test_intersection_with_covered_line_segment() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
        PrimitiveIntersections::Zero
    );
    approx::assert_abs_diff_eq!(
        ls.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
        PrimitiveIntersections::Zero
    );
}

#[test]
fn test_intersection_arc_segment() {
    {
        // Regression test from stem_slot crate
        let line = Line {
            a: -0.0,
            b: 1.0,
            c: -0.0020953915850751483,
        };
        let arc = ArcSegment::new(
            [0.0025088728825159476, 0.0029999999999999936],
            0.0009999999999999935,
            4.71238898038469,
            1.4835298641951864,
        )
        .unwrap();
        approx::assert_abs_diff_eq!(
            arc.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE),
            PrimitiveIntersections::One([0.002935116493197851, 0.0020953915850751483]),
            epsilon = 1e-6
        );
    }
}

#[test]
fn test_transformation() {
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        line.translate([2.0, -2.0]);
        assert!(line.covers_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers_point([0.0, 4.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        line.rotate([2.0, 2.0], 0.5 * PI);
        assert!(line.covers_point([0.0, 4.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let mut line = Line::from_point_angle([1.0, 1.0], -0.25 * PI);
        assert!(!line.covers_point([2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        line.scale(2.0);
        assert!(line.covers_point([2.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        line.line_reflection([-1.0, 0.0], [1.0, 0.0]);
        assert!(line.covers_point([2.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers_point([0.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
        line.point_reflection([1.0, 0.0]);
        assert!(line.covers_point([0.0, -2.0], DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    }
}
