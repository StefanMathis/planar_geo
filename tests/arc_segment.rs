use approx::assert_abs_diff_eq;
use std::f64::consts::{FRAC_PI_2, PI, SQRT_2, TAU};

use planar_geo::prelude::*;

#[test]
fn test_contains_point() {
    {
        // Check for rounding error: Point [0.5, 0.0] is almost identical to the start point of the arc
        let arc = ArcSegment::from_start_middle_stop(
            [0.49999999999999967, 0.0],
            [0.8535533905932736, 0.14644660940672583],
            [1.0000000000000002, 0.49999999999999944],
            0.0,
            0,
        )
        .unwrap();
        assert!(arc.contains_point([0.5, 0.0], 0.0, 0));

        let center = [0.0, 1.0];
        let radius = 1.0;
        let start_angle = -4.0 * std::f64::consts::PI; // This is mathematically identical to 0.0
        let offset_angle = -2.5 * std::f64::consts::PI; // This is mathematically identical to -pi/2
        let arc = ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            0.0,
            0,
        )
        .unwrap();
        assert!(arc.contains_angle(-FRAC_PI_2 + 1e-15));

        let p1 = [1.0, 1.0];
        let p2 = [0.0, 0.0];
        let p3 = [1.5, 0.5];
        assert!(arc.contains_point(p1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(arc.contains_point(p2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(!arc.contains_point(p3, DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let pt1 = [-0.9233365375660403, 9.999999999996123e-6];
        let pt2 = [-0.9226294661406899, 2.500031251961854e-6];
        let pt3 = [-0.9219223593595585, 1.582067810090848e-15];
        let arc =
            ArcSegment::from_start_middle_stop(pt1.into(), pt2.into(), pt3.into(), 0.0, 0).unwrap();

        // Due to floating point rounding errors, some points are not "exactly" included
        assert!(!arc.contains_point(pt1.into(), 0.0, 0));
        assert!(!arc.contains_point(pt2.into(), 0.0, 0));
        assert!(!arc.contains_point(pt3.into(), 0.0, 0));
        assert!(arc.contains_point(pt1.into(), DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(arc.contains_point(pt2.into(), DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
        assert!(arc.contains_point(pt3.into(), DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
}

#[test]
fn test_arc_arc_intersection() {
    {
        let center = [0.0, 0.0];
        let radius = 1.0;
        let start_angle = 0.0;
        let offset_angle = 0.5 * std::f64::consts::PI;
        let arc1 = ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            0.0,
            0,
        )
        .unwrap();

        let center = [0.0, 1.0];
        let radius = 1.0;
        let start_angle = 0.0;
        let offset_angle = -0.5 * std::f64::consts::PI;
        let arc2 = ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            0.0,
            0,
        )
        .unwrap();

        let intersections = arc1.intersections_primitive(&arc2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        assert_abs_diff_eq!(
            PrimitiveIntersections::One([3.0_f64.sqrt() / 2.0, 0.5]),
            intersections,
        );
    }
    {
        let center = [0.0, 0.0];
        let radius = 1.0;
        let start_angle = 0.0;
        let offset_angle = TAU;
        let arc1 = ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            0.0,
            0,
        )
        .unwrap();

        let center = [0.0, 1.0];
        let radius = 1.0;
        let start_angle = 0.0;
        let offset_angle = TAU;
        let arc2 = ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            0.0,
            0,
        )
        .unwrap();

        let intersections = arc1.intersections_primitive(&arc2, 0.0, 0);
        assert_eq!(intersections.len(), 2);
        assert_abs_diff_eq!(
            PrimitiveIntersections::Two([
                [-3.0_f64.sqrt() / 2.0, 0.5],
                [3.0_f64.sqrt() / 2.0, 0.5]
            ]),
            intersections,
        );
    }
    {
        let arc1 = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.0,
            PI,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line = Line::from_point_angle([0.0, 0.0], 0.5 * PI);
        let arc2 = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.5 * PI,
            PI,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();

        // Regular intersection
        assert_abs_diff_eq!(
            line.intersections_arc_segment(&arc1, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::One([0.0, 2.0])
        );

        // Degenerate case
        assert_abs_diff_eq!(
            arc1.intersections_arc_segment(&arc2, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
            PrimitiveIntersections::Two([[0.0, 2.0], [-2.0, 0.0]]),
            epsilon = DEFAULT_EPSILON
        );
    }
}

#[test]
fn test_line_segment_arc_intersection() {
    {
        // Arc is touching a line segment
        let arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            1.0,
            0.0,
            std::f64::consts::FRAC_PI_2,
            0.0,
            0,
        )
        .unwrap()
        .into();
        let line: Segment =
            LineSegment::new([1.0, -10.0], [1.0, 10.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let intersections = arc.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.len(), 1);
        assert_abs_diff_eq!(PrimitiveIntersections::One([1.0, 0.0]), intersections,);
    }
    {
        let arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
            [0.5000000000000001, 0.5000000000000002],
            0.5000000000000002,
            4.71238898038469,
            1.5707963267948966,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap()
        .into();
        let line: Segment =
            LineSegment::new([0.5, 1.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let intersections = arc.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.len(), 1);
        assert_abs_diff_eq!(
            PrimitiveIntersections::One([0.9000000000000001, 0.19999999999999984]),
            intersections,
        );
    }
    {
        let arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            20.0,
            0.0,
            6.283185307179586,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap()
        .into();
        let line: Segment = LineSegment::new(
            [-21.0, -0.9999923706054688],
            [0.055555555555542924, -0.9999923706054688],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap()
        .into();
        let intersections = arc.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
        assert_eq!(intersections.len(), 1);
    }
}

#[test]
fn self_intersection() {
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        2.0,
        0.0,
        PI,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();

    // Self-intersection
    assert_eq!(
        arc.intersections_primitive(&arc, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );

    // Intersections with equal primitive
    let arc_cloned = arc.clone();
    assert_abs_diff_eq!(
        arc.intersections_primitive(&arc_cloned, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[2.0, 0.0], [-2.0, 0.0]]),
        epsilon = DEFAULT_EPSILON
    );
}

#[test]
fn test_polygonize_iter_count() {
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        1.0,
        0.0,
        0.75 * PI,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();

    {
        let iter = arc.polygonize(SegmentPolygonizer::MaximumSegmentLength(10.0));
        assert_eq!(iter.count(), 2);

        let iter = arc.polygonize(SegmentPolygonizer::MaximumSegmentLength(-10.0));
        assert_eq!(iter.count(), 1);

        let iter = arc.polygonize(SegmentPolygonizer::MaximumSegmentLength(-1.0));
        assert_eq!(iter.count(), 1);
    }
    {
        let iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(0.5 * PI));
        assert_eq!(iter.count(), 3);
        let iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(0.0));
        assert_eq!(iter.count(), 1);
        let iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(-1.0));
        assert_eq!(iter.count(), 1);
    }
    {
        let iter = arc.polygonize(SegmentPolygonizer::NumberSegments(3));
        assert_eq!(iter.count(), 4);

        let iter = arc.polygonize(SegmentPolygonizer::NumberSegments(0));
        assert_eq!(iter.count(), 1);
    }
}

#[test]
fn test_polygonize_segment_length() {
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        [0.0, 0.0],
        1.0,
        0.0,
        0.75 * PI,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();

    let iter = arc.polygonize(SegmentPolygonizer::MaximumSegmentLength(0.5));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    assert_eq!(chain.num_segments(), 5);

    for segments in chain.as_slices().0.windows(2) {
        assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }

    let iter = arc.polygonize(SegmentPolygonizer::NumberSegments(8));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    for segments in chain.as_slices().0.windows(2) {
        assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }

    let iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(0.1));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    for segments in chain.as_slices().0.windows(2) {
        assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }
}

#[test]
fn test_polygonize_number_segments() {
    let center = [0.0, 0.0];
    let radius = 1.0;
    let start_angle = -4.0 * PI; // This is mathematically identical to 0.0
    let offset_angle = 2.75 * PI; // This is mathematically identical to 3/4*PI
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();

    let option = SegmentPolygonizer::NumberSegments(3);
    let pts: Vec<[f64; 2]> = arc.polygonize(option).collect();
    assert_eq!(pts.len(), 4);
    assert_abs_diff_eq!(&pts[0], &[1.0, 0.0], epsilon = 1e-6);
    assert_abs_diff_eq!(&pts[1], &[0.7071067, 0.7071067], epsilon = 1e-6);
    assert_abs_diff_eq!(&pts[2], &[0.0, 1.0], epsilon = 1e-6);
    assert_abs_diff_eq!(&pts[3], &[-0.707106, 0.707106], epsilon = 1e-6);
}

#[test]
fn test_polygonize_max_arc() {
    let center = [1.0, 2.0];
    let radius = 2.0;
    let start_angle = FRAC_PI_2;
    let offset_angle = PI;
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();

    let option = SegmentPolygonizer::MaximumAngle(FRAC_PI_2 + 0.1);
    let pts: Vec<[f64; 2]> = arc.polygonize(option).collect();
    assert_eq!(pts.len(), 3);
    assert_abs_diff_eq!(&pts[0], &[1.0, 4.0], epsilon = 1e-15);
    assert_abs_diff_eq!(&pts[1], &[-1.0, 2.0], epsilon = 1e-15);
    assert_abs_diff_eq!(&pts[2], &[1.0, 0.0], epsilon = 1e-15);
}

#[test]
fn test_polygonize_max_segment_length() {
    let center = [1.0, 2.0];
    let radius = 2.0;
    let start_angle = FRAC_PI_2;
    let offset_angle = PI;
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();
    assert_abs_diff_eq!(arc.length(), TAU);

    let option = SegmentPolygonizer::MaximumSegmentLength(3.0);
    let pts: Vec<[f64; 2]> = arc.polygonize(option).collect();
    assert_eq!(pts.len(), 3);
    assert_abs_diff_eq!(&pts[0], &[1.0, 4.0]);
    assert_abs_diff_eq!(&pts[1], &[-1.0, 2.0], epsilon = 1e-15);
    assert_abs_diff_eq!(&pts[2], &[1.0, 0.0], epsilon = 1e-15);

    let option = SegmentPolygonizer::MaximumSegmentLength(1.8);
    let pts: Vec<[f64; 2]> = arc.polygonize(option).collect();
    assert_eq!(pts.len(), 5);
    assert_abs_diff_eq!(&pts[0], &[1.0, 4.0], epsilon = 1e-15);
    assert_abs_diff_eq!(&pts[1], &[-0.414213, 3.414213], epsilon = 1e-6);
    assert_abs_diff_eq!(&pts[2], &[-1.0, 2.0], epsilon = 1e-15);
    assert_abs_diff_eq!(&pts[3], &[-0.414213, 0.585786], epsilon = 1e-6);
    assert_abs_diff_eq!(&pts[4], &[1.0, 0.0], epsilon = 1e-15);
}

#[test]
fn test_construct_arc() {
    let center = [0.0, 0.0];
    let radius = 1.0;
    let start_angle = -4.0 * PI; // This is mathematically identical to 0.0
    let offset_angle = -2.5 * PI; // This is mathematically identical to -pi/2
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();
    assert_abs_diff_eq!(arc.start_angle(), 0.0);
    assert_abs_diff_eq!(arc.offset_angle(), -FRAC_PI_2);
    assert_abs_diff_eq!(arc.start(), [1.0, 0.0]);
    assert_abs_diff_eq!(arc.stop(), [0.0, -1.0]);

    let center = [0.0, 0.0];
    let radius = 1.0;
    let start_angle = -4.0 * PI;
    let offset_angle = 2.75 * PI;
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();
    assert_abs_diff_eq!(arc.length(), 0.75 * PI, epsilon = 1e-15);
    assert_abs_diff_eq!(arc.start(), [1.0, 0.0], epsilon = 1e-15);
    assert_abs_diff_eq!(
        arc.stop(),
        [-1.0 / 2.0f64.sqrt(), 1.0 / 2.0f64.sqrt()],
        epsilon = 1e-15
    );

    let center = [0.0, 0.0];
    let radius = 1.0;
    let start_angle = 0.5 * PI; // This is mathematically identical to 0.0
    let arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        -0.5 * PI,
        0.0,
        0,
    )
    .unwrap();
    assert_abs_diff_eq!(arc.length(), 0.5 * PI);
    assert_abs_diff_eq!(arc.start(), [0.0, 1.0]);
    assert_abs_diff_eq!(arc.stop(), [1.0, 0.0]);
}

#[test]
fn test_construct_circle() {
    let circle = ArcSegment::circle([0.5, 0.5], 2.0).unwrap();
    assert_abs_diff_eq!(circle.center(), [0.5, 0.5]);
    assert_abs_diff_eq!(circle.radius(), 2.0);
    assert_abs_diff_eq!(circle.start_angle(), 0.0);
    assert_abs_diff_eq!(circle.stop_angle(), TAU);
    assert_abs_diff_eq!(circle.start(), [2.5, 0.5]);
    assert_abs_diff_eq!(circle.stop(), [2.5, 0.5], epsilon = 1e-3);
    assert_abs_diff_eq!(circle.offset_angle(), TAU);
}

#[test]
fn test_scale() {
    let center = [1.0, 0.0];
    let radius = 3.0;
    let start_angle = 0.5 * PI; // This is mathematically identical to 0.0
    let offset_angle = PI; // This is mathematically identical to -pi/2
    let mut arc = ArcSegment::from_center_radius_start_offset_angle(
        center,
        radius,
        start_angle,
        offset_angle,
        0.0,
        0,
    )
    .unwrap();
    assert_abs_diff_eq!(arc.center(), [1.0, 0.0]);
    assert_abs_diff_eq!(arc.radius(), 3.0);
    assert_abs_diff_eq!(arc.start(), [1.0, 3.0]);
    assert_abs_diff_eq!(arc.stop(), [1.0, -3.0], epsilon = 1e-3);
    assert_abs_diff_eq!(arc.start_angle(), 0.5 * PI);
    assert_abs_diff_eq!(arc.stop_angle(), 1.5 * PI);
    assert_abs_diff_eq!(arc.offset_angle(), PI);

    // Scale the arc
    arc.scale(2.0);
    assert_abs_diff_eq!(arc.center(), [2.0, 0.0]);
    assert_abs_diff_eq!(arc.radius(), 6.0);
    assert_abs_diff_eq!(arc.start(), [2.0, 6.0], epsilon = 1e-15);
    assert_abs_diff_eq!(arc.stop(), [2.0, -6.0], epsilon = 1e-10);
    assert_abs_diff_eq!(arc.start_angle(), 0.5 * PI);
    assert_abs_diff_eq!(arc.stop_angle(), 1.5 * PI);
    assert_abs_diff_eq!(arc.offset_angle(), PI);
}

#[test]
fn test_arc_from_start_middle_stop() {
    let start = [0.0, 1.0];
    let middle = [1.0 - 1.0 / 2.0_f64.sqrt(), 1.0 - 1.0 / 2.0_f64.sqrt()];
    let stop = [1.0, 0.0];
    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, 0.0, 0).unwrap();

    assert_abs_diff_eq!(arc.center(), [1.0, 1.0]);
    assert_abs_diff_eq!(arc.radius(), 1.0);
    assert_abs_diff_eq!(arc.start_angle(), std::f64::consts::PI);
    assert_abs_diff_eq!(arc.offset_angle(), 0.5 * std::f64::consts::PI);
}

#[test]
fn test_fillet_construction() {
    {
        let prev: [f64; 2] = [0.0, 1.0];
        let curr: [f64; 2] = [0.0, 0.0];
        let next: [f64; 2] = [1.0, 0.0];
        let arc = ArcSegment::fillet(prev, curr, next, 0.5, 0.0, 0).unwrap();
        assert_abs_diff_eq!(arc.start(), &[0.0, 0.5]);
        assert_abs_diff_eq!(
            arc.segment_point(0.5),
            [0.5 - 0.5 / 2.0_f64.sqrt(), 0.5 - 0.5 / 2.0_f64.sqrt()]
        );
        assert_abs_diff_eq!(arc.stop(), &[0.5, 0.0]);

        let prev: [f64; 2] = [0.0, 0.0];
        let curr: [f64; 2] = [1.0, 0.0];
        let next: [f64; 2] = [1.0, 1.0];
        let arc = ArcSegment::fillet(prev, curr, next, 0.5, 0.0, 0).unwrap();
        assert_abs_diff_eq!(arc.start(), &[0.5, 0.0]);
        assert_abs_diff_eq!(
            arc.segment_point(0.5),
            [0.5 + 0.5 / 2.0_f64.sqrt(), 0.5 - 0.5 / 2.0_f64.sqrt()],
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(arc.stop(), &[1.0, 0.5], epsilon = 1e-6);
    }
    {
        let prev: [f64; 2] = [1.0, 2.0];
        let curr: [f64; 2] = [0.0, 2.0];
        let next: [f64; 2] = [0.0, 1.0];
        assert!(ArcSegment::fillet(prev, curr, next, 0.0, 0.0, 0).is_err());
    }
}

#[test]
fn test_three_point_arc_almost_on_line() {
    let p = [2.098837955910367, 6.971002742420203];
    let q = [1.0847897424142476, 5.985251140491335];
    let r = [0.07074152891812842, 4.999499538562468];
    assert!(three_point_arc_center(p, q, r).is_ok());
    assert!(three_point_arc_center(q, r, p).is_ok());
    assert!(three_point_arc_center(r, p, q).is_ok());
    assert!(three_point_arc_center(p, r, q).is_ok());
    assert!(three_point_arc_center(r, q, p).is_ok());
    assert!(three_point_arc_center(q, p, r).is_ok());
}

#[test]
fn test_line_reflection() {
    {
        let mut arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.0,
            0.5 * PI,
            0.0,
            0,
        )
        .unwrap();
        assert_abs_diff_eq!(arc.start(), [2.0, 0.0]);
        assert_abs_diff_eq!(arc.stop(), [0.0, 2.0]);

        arc.line_reflection([0.0, 0.0], [0.0, 1.0]);
        assert_abs_diff_eq!(arc.start(), [-2.0, 0.0], epsilon = 1e-15);
        assert_abs_diff_eq!(arc.stop(), [0.0, 2.0], epsilon = 1e-15);
    }
    {
        let mut arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.25 * PI,
            -0.5 * PI,
            0.0,
            0,
        )
        .unwrap();
        assert_abs_diff_eq!(arc.start(), [2.0f64.sqrt(), 2.0f64.sqrt()], epsilon = 1e-15,);
        assert_abs_diff_eq!(arc.stop(), [2.0f64.sqrt(), -2.0f64.sqrt()], epsilon = 1e-15,);
        assert_eq!(arc.offset_angle(), -0.5 * PI);

        arc.line_reflection([0.0, 0.0], [0.0, 1.0]);
        assert_abs_diff_eq!(
            arc.start(),
            [-2.0f64.sqrt(), 2.0f64.sqrt()],
            epsilon = 1e-15,
        );
        assert_abs_diff_eq!(
            arc.stop(),
            [-2.0f64.sqrt(), -2.0f64.sqrt()],
            epsilon = 1e-15,
        );
        assert_eq!(arc.offset_angle(), 0.5 * PI);
    }
    {
        let mut arc = ArcSegment::from_center_radius_start_offset_angle(
            [1.0, 1.0],
            2.0,
            0.5 * PI,
            0.5 * PI,
            0.0,
            0,
        )
        .unwrap();
        assert_abs_diff_eq!(arc.start(), [1.0, 3.0], epsilon = 1e-15);
        assert_abs_diff_eq!(arc.stop(), [-1.0, 1.0], epsilon = 1e-15);

        arc.line_reflection([0.0, 0.0], [1.0, 1.0]);
        assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
        assert_abs_diff_eq!(arc.stop(), [1.0, -1.0], epsilon = 1e-15);
    }
}

#[test]
fn test_bounding_box() {
    {
        let arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.25 * PI,
            0.1 * PI,
            0.0,
            0,
        )
        .unwrap();

        let bb = BoundingBox::from(&arc);
        assert_abs_diff_eq!(bb.xmin(), 0.90798, epsilon = 1e-4);
        assert_abs_diff_eq!(bb.xmax(), 2.0f64.sqrt());
        assert_abs_diff_eq!(bb.ymin(), 2.0f64.sqrt());
        assert_abs_diff_eq!(bb.ymax(), 1.78201, epsilon = 1e-4);
    }
    {
        let arc = ArcSegment::circle([1.0, 1.0], 3.0).unwrap();

        let bb = BoundingBox::from(&arc);
        assert_abs_diff_eq!(bb.xmin(), -2.0);
        assert_abs_diff_eq!(bb.xmax(), 4.0);
        assert_abs_diff_eq!(bb.ymin(), -2.0);
        assert_abs_diff_eq!(bb.ymax(), 4.0);
    }
    {
        let arc = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.25 * PI,
            PI,
            0.0,
            0,
        )
        .unwrap();

        let bb = BoundingBox::from(&arc);
        assert_abs_diff_eq!(bb.xmin(), -2.0);
        assert_abs_diff_eq!(bb.xmax(), 2.0f64.sqrt());
        assert_abs_diff_eq!(bb.ymin(), -2.0f64.sqrt());
        assert_abs_diff_eq!(bb.ymax(), 2.0);
    }
}

#[test]
fn test_vertex_iterator() {
    {
        let line: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, TAU, 0.0, 0)
                .unwrap()
                .into();
        let mut iter = line.polygonize(SegmentPolygonizer::MaximumAngle(FRAC_PI_2));
        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }

    {
        let line: Segment = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.5 * PI,
            PI,
            0.0,
            0,
        )
        .unwrap()
        .into();
        let mut iter = line.polygonize(SegmentPolygonizer::MaximumAngle(0.25 * PI));
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(
            iter.next().unwrap(),
            [-2.0f64.sqrt(), 2.0f64.sqrt()],
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(
            iter.next().unwrap(),
            [-2.0f64.sqrt(), -2.0f64.sqrt()],
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }
    {
        let arc: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, TAU, 0.0, 0)
                .unwrap()
                .into();

        let mut iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(FRAC_PI_2));

        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }

    {
        let arc: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, -TAU, 0.0, 0)
                .unwrap()
                .into();

        let mut iter = arc.polygonize(SegmentPolygonizer::MaximumAngle(FRAC_PI_2));

        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }
}

#[test]
fn test_segment_point() {
    {
        let arc: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, PI, 0.0, 0)
                .unwrap()
                .into();
        assert_abs_diff_eq!(arc.segment_point(0.25), [SQRT_2, SQRT_2]);
        assert_abs_diff_eq!(arc.segment_point(0.75), [-SQRT_2, SQRT_2]);
        assert_abs_diff_eq!(
            arc.segment_point(0.1),
            [2.0 * (0.1 * PI).cos(), 2.0 * (0.1 * PI).sin()]
        );
        assert_abs_diff_eq!(arc.segment_point(-0.25), [2.0, 0.0]);
        assert_abs_diff_eq!(arc.segment_point(1.25), [-2.0, 0.0], epsilon = 1e-15);
    }
    {
        let arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.5 * PI,
            -0.5 * PI,
            0.0,
            0,
        )
        .unwrap()
        .into();
        assert_abs_diff_eq!(arc.segment_point(0.0), [0.0, 2.0]);
        assert_abs_diff_eq!(arc.segment_point(1.0), [2.0, 0.0]);
        assert_abs_diff_eq!(arc.segment_point(2.0), [2.0, 0.0]);
        assert_abs_diff_eq!(arc.segment_point(0.5), [SQRT_2, SQRT_2]);
    }
}

#[test]
fn test_contains_angle() {
    let arc = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, TAU, 0.0, 0)
        .unwrap();

    // This arc contains every angle
    assert!(arc.contains_angle(-0.05002047485831401));
}
