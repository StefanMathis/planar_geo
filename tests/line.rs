use std::f64::consts::{FRAC_PI_4, PI};

use approx;
use planar_geo::prelude::*;

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
        assert!(line.covers(&[1.0, 1.0]));
        assert!(line.covers(&[2.0, 1.0]));
        assert!(line.covers(&[100000.0, 1.0]));
        assert!(!line.covers(&[100000.0, 0.9]));
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -FRAC_PI_4);
        assert!(line.covers(&[1.0, -1.0]));
        assert!(!line.covers(&[1.0, 0.0]));
        assert!(line.with_tolerance(1.0, 0.0).covers(&[1.0, 0.0]));
    }
    {
        let line = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        assert!(line.covers(&[1.0, -3.0f64.sqrt()]));
    }
    {
        let line = Line::from_point_angle([0.0, -2.0], PI / 3.0);
        assert!(line.covers(&[1.0, 3.0f64.sqrt() - 2.0],));
    }
    {
        // Vertical line
        let line = LineSegment::new([0.5, -0.5], [0.5, 0.5]).unwrap();
        assert!(line.with_tolerance(0.0, 0.0).covers(&[0.5, 0.0]));
    }
    {
        // Vertical line
        let line = LineSegment::new([0.5, -0.5], [0.5, 0.5]).unwrap();
        assert!(!line.with_tolerance(0.0, 0.0).covers(&[0.5, 1.0]));
    }
    {
        let line = LineSegment::new([1.0, 0.0], [0.0, 1.0]).unwrap();
        assert!(line.with_tolerance(0.0, 0.0).covers(&[0.5, 0.5]));
    }
    {
        // This point is clearly not directly on the line
        let line = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
        let point = [0.5, 1.0];
        assert!(!line.with_tolerance(0.0, 0.0).covers(&point));

        // But it is within this tolerance
        assert!(line.with_tolerance(1.0, 0.0).covers(&point));

        // As well as inside this tolerance
        assert!(line.with_tolerance(1.00001, 0.0).covers(&point,));

        // But not outside this tolerance
        assert!(!line.with_tolerance(0.99999, 0.0).covers(&point));

        // Repeat the test with a point located on a straight extension of the original
        // line
        let point = [2.0, 0.0];
        assert!(!line.with_tolerance(0.0, 0.0).covers(&point));
        assert!(line.with_tolerance(1.0, 0.0).covers(&point,));
        assert!(line.with_tolerance(1.00001, 0.0).covers(&point,));
        assert!(!line.with_tolerance(0.99999, 0.0).covers(&point));

        // This point is outside of the circular radius around the endpoint [1.0, 0.0]
        assert!(!line.with_tolerance(0.11, 0.0).covers(&[1.1, -0.1]));
        assert!(line.with_tolerance(0.11, 0.0).covers(&[0.5, -0.1]));
        assert!(!line.with_tolerance(0.1, 0.0).covers(&[1.1, -0.1]));
        assert!(!line.with_tolerance(0.09, 0.0).covers(&[1.1, -0.1]));
    }
    {
        let line = LineSegment::new([0.0, 0.0], [1000.0, 0.0]).unwrap();

        // This point is outside of the circular radius around the endpoint [1000.0,
        // 0.0]
        assert!(!line.with_tolerance(110.0, 0.0).covers(&[1100.0, -100.0],));
        assert!(!line.with_tolerance(100.0, 0.0).covers(&[1100.0, -100.0],));
    }
}

#[test]
fn test_covers_line_segment() {
    let ls1 = LineSegment::new([0.0, 0.0], [3.0, 0.0]).expect("valid inputs");
    let line = Line::from(&ls1);
    let ls2 = LineSegment::new([1.0, 0.0], [2.0, 0.0]).expect("valid inputs");
    assert!(line.covers(&ls1));
    assert!(line.covers(&ls2));
}

#[test]
fn test_identical() {
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([-2.0, 0.0], [2.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(!first.identical(&second))
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([1.0, 1.0], [3.0, 3.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second));
        assert!(first.identical(&second))
    }
}

#[test]
fn test_intersection_line_line() {
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 1.0], [2.0, 3.0]).unwrap();
        let second = Line::from(&line_segment);

        assert!(first.parallel(&second));
        assert_eq!(
            first.intersections_primitive(&second),
            PrimitiveIntersections::Zero
        );
        assert_eq!(first.intersections(&second), Vec::new());
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [0.0, 1.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second),
            PrimitiveIntersections::One([0.0, 0.0])
        );
        assert_eq!(
            first.intersections(&second),
            vec![Intersection {
                point: [0.0, 0.0],
                left: SegmentKey::from_segment_idx(0),
                right: SegmentKey::from_segment_idx(0),
            }]
        );
    }
    {
        let line_segment = LineSegment::new([3.0, 0.0], [3.0, 1.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([0.0, 2.0], [1.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second),
            PrimitiveIntersections::One([3.0, 2.0])
        );
    }
    {
        let line_segment = LineSegment::new([0.0, 0.0], [2.0, 2.0]).unwrap();
        let first = Line::from(&line_segment);

        let line_segment = LineSegment::new([-2.0, 0.0], [2.0, 2.0]).unwrap();
        let second = Line::from(&line_segment);

        assert_eq!(
            first.intersections_primitive(&second),
            PrimitiveIntersections::One([2.0, 2.0])
        );
    }
    {
        let first = Line::from_point_angle([0.0, 0.0], -PI / 3.0);
        let second = Line::from_point_angle([0.0, -2.0], PI / 3.0);

        approx::assert_abs_diff_eq!(
            first.intersections_primitive(&second),
            PrimitiveIntersections::One([1.0 / 3.0f64.sqrt(), -1.0])
        );
    }
    {
        let line_1 = Line::from_point_angle([0.11, 0.11], 0.0);
        let line_2 = Line::from_point_angle([0.11, 0.11], std::f64::consts::FRAC_PI_2);
        assert_eq!(line_1.intersections_primitive(&line_2).len(), 1);
        assert_eq!(line_2.intersections_primitive(&line_1).len(), 1);
    }
}

#[test]
fn test_self_intersection() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&line),
        PrimitiveIntersections::Zero
    );
}

#[test]
fn test_intersection_with_covered_line_segment() {
    let line = Line::from_point_angle([0.0, 0.0], 0.0);
    let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
    approx::assert_abs_diff_eq!(
        line.intersections_primitive(&ls),
        PrimitiveIntersections::Zero
    );
    approx::assert_abs_diff_eq!(
        ls.intersections_primitive(&line),
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
            arc.intersections_primitive(&line),
            PrimitiveIntersections::One([0.002935116493197851, 0.0020953915850751483]),
            epsilon = 1e-6
        );
    }
}

#[test]
fn test_transformation() {
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers(&[2.0, -2.0]));
        line.translate([2.0, -2.0]);
        assert!(line.covers(&[2.0, -2.0]));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers(&[0.0, 4.0]));
        line.rotate([2.0, 2.0], 0.5 * PI);
        assert!(line.covers(&[0.0, 4.0]));
    }
    {
        let mut line = Line::from_point_angle([1.0, 1.0], -0.25 * PI);
        assert!(!line.covers(&[2.0, 2.0]));
        line.scale(2.0);
        assert!(line.covers(&[2.0, 2.0]));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers(&[2.0, -2.0]));
        line.line_reflection([-1.0, 0.0], [1.0, 0.0]);
        assert!(line.covers(&[2.0, -2.0]));
        assert!(line.covers(&[2.0, -2.0]));
    }
    {
        let mut line = Line::from_point_angle([0.0, 0.0], 0.25 * PI);
        assert!(!line.covers(&[0.0, -2.0]));
        line.point_reflection([1.0, 0.0]);
        assert!(line.covers(&[0.0, -2.0]));
        assert!(line.covers(&[0.0, -2.0]));
    }
}
