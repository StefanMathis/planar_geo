use std::f64::{INFINITY, NEG_INFINITY, consts::PI};

use planar_geo::prelude::*;

#[test]
fn test_polygonize_number_segments() {
    let start = [0.0, 0.0];
    let stop = [1.0, -1.0];
    let line = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let option = SegmentPolygonizer::InnerSegments(4);
    let pts: Vec<[f64; 2]> = line.polygonize(option).collect();
    assert_eq!(pts.len(), 5);
    approx::assert_abs_diff_eq!(&pts[0], &[0.0, 0.0]);
    approx::assert_abs_diff_eq!(&pts[1], &[0.25, -0.25]);
    approx::assert_abs_diff_eq!(&pts[2], &[0.5, -0.5]);
    approx::assert_abs_diff_eq!(&pts[3], &[0.75, -0.75]);
    approx::assert_abs_diff_eq!(&pts[4], &[1.0, -1.0]);

    // Vertical line
    let start = [0.0, 0.0];
    let stop = [0.0, -1.0];
    let line = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let option = SegmentPolygonizer::InnerSegments(4);
    let pts: Vec<[f64; 2]> = line.polygonize(option).collect();
    assert_eq!(pts.len(), 5);
    approx::assert_abs_diff_eq!(&pts[0], &[0.0, 0.0]);
    approx::assert_abs_diff_eq!(&pts[1], &[0.0, -0.25]);
    approx::assert_abs_diff_eq!(&pts[2], &[0.0, -0.5]);
    approx::assert_abs_diff_eq!(&pts[3], &[0.0, -0.75]);
    approx::assert_abs_diff_eq!(&pts[4], &[0.0, -1.0]);
}

#[test]
fn test_polygonize_max_arc() {
    let start = [0.0, 0.0];
    let stop = [1.0, -1.0];
    let line = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let option = SegmentPolygonizer::MaximumAngle(0.5);
    let pts: Vec<[f64; 2]> = line.polygonize(option).collect();
    assert_eq!(pts.len(), 1);
}

#[test]
fn test_polygonize_max_length() {
    let start = [0.0, 0.0];
    let stop = [0.0, -1.0];
    let line = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let option = SegmentPolygonizer::MaximumSegmentLength(0.21);
    let pts: Vec<[f64; 2]> = line.polygonize(option).collect();
    assert_eq!(pts.len(), 6);
    approx::assert_abs_diff_eq!(&pts[0], &[0.0, 0.0]);
    approx::assert_abs_diff_eq!(&pts[1], &[0.0, -0.2]);
    approx::assert_abs_diff_eq!(&pts[2], &[0.0, -0.4]);
    approx::assert_abs_diff_eq!(&pts[3], &[0.0, -0.6]);
    approx::assert_abs_diff_eq!(&pts[4], &[0.0, -0.8]);
    approx::assert_abs_diff_eq!(&pts[5], &[0.0, -1.0]);
}

#[test]
fn test_line_min_max() {
    let start = [0.0, 0.0];
    let stop = [1.0, -1.0];
    let line = LineSegment::new(start, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    assert_eq!(line.xmin(), 0.0);
    assert_eq!(line.xmax(), 1.0);
    assert_eq!(line.ymin(), -1.0);
    assert_eq!(line.ymax(), 0.0);
}

#[test]
fn test_bounding_box() {
    {
        let line =
            LineSegment::new([0.0, 2.0], [1.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let bb = BoundingBox::from(&line);
        assert_eq!(bb.xmin(), 0.0);
        assert_eq!(bb.xmax(), 1.0);
        assert_eq!(bb.ymin(), -1.0);
        assert_eq!(bb.ymax(), 2.0);
    }
    {
        let line =
            LineSegment::new([1.0, -1.0], [0.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let bb = BoundingBox::from(&line);
        assert_eq!(bb.xmin(), 0.0);
        assert_eq!(bb.xmax(), 1.0);
        assert_eq!(bb.ymin(), -1.0);
        assert_eq!(bb.ymax(), 2.0);
    }
    {
        let line =
            LineSegment::new([1.0, -1.0], [0.0, -1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let bb = BoundingBox::from(&line);
        assert_eq!(bb.xmin(), 0.0);
        assert_eq!(bb.xmax(), 1.0);
        assert_eq!(bb.ymin(), -1.0);
        assert_eq!(bb.ymax(), -1.0);
    }
    {
        let line =
            LineSegment::new([1.0, -1.0], [1.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let bb = BoundingBox::from(&line);
        assert_eq!(bb.xmin(), 1.0);
        assert_eq!(bb.xmax(), 1.0);
        assert_eq!(bb.ymin(), -1.0);
        assert_eq!(bb.ymax(), 2.0);
    }
}

#[test]
fn test_line_reflection() {
    {
        let mut line =
            LineSegment::new([1.0, -1.0], [1.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        line.line_reflection([1.0, -0.5], [-1.0, -0.5]);
        assert_eq!(line.start(), [1.0, 0.0]);
        assert_eq!(line.stop(), [1.0, -3.0]);
    }
    {
        let mut line =
            LineSegment::new([1.0, -1.0], [1.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        line.line_reflection([0.0, 0.0], [1.0, 1.0]);
        assert_eq!(line.start(), [-1.0, 1.0]);
        assert_eq!(line.stop(), [2.0, 1.0]);
    }
}

#[test]
fn test_rotate() {
    let mut line =
        LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    line.rotate([0.0, 0.0], PI);
    approx::assert_abs_diff_eq!(line.start(), [0.0, 0.0]);
    approx::assert_abs_diff_eq!(line.stop(), [0.0, -1.0]);
}

#[test]
fn test_polygonize_points() {
    {
        let line: Segment =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let mut iter = line.polygonize(SegmentPolygonizer::InnerSegments(1));
        assert_eq!(iter.next(), Some([0.0, 0.0]));
        assert_eq!(iter.next(), Some([1.0, 0.0]));
        assert_eq!(iter.next(), None);
    }
    {
        let line: Segment =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let mut iter = line.polygonize(SegmentPolygonizer::MaximumSegmentLength(0.3));
        assert_eq!(iter.next(), Some([0.0, 0.0]));
        assert_eq!(iter.next(), Some([0.25, 0.0]));
        assert_eq!(iter.next(), Some([0.5, 0.0]));
        assert_eq!(iter.next(), Some([0.75, 0.0]));
        assert_eq!(iter.next(), Some([1.0, 0.0]));
        assert_eq!(iter.next(), None);
    }
}

#[test]
fn test_polygonize_iter_count() {
    let line: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
        .unwrap()
        .into();
    {
        let iter = line.polygonize(SegmentPolygonizer::MaximumSegmentLength(10.0));
        assert_eq!(iter.count(), 2);

        let iter = line.polygonize(SegmentPolygonizer::MaximumSegmentLength(-10.0));
        assert_eq!(iter.count(), 1);

        let iter = line.polygonize(SegmentPolygonizer::MaximumSegmentLength(-1.0));
        assert_eq!(iter.count(), 1);
    }
    {
        let iter = line.polygonize(SegmentPolygonizer::MaximumAngle(0.5 * PI));
        assert_eq!(iter.count(), 1);
        let iter = line.polygonize(SegmentPolygonizer::MaximumAngle(0.0));
        assert_eq!(iter.count(), 1);
        let iter = line.polygonize(SegmentPolygonizer::MaximumAngle(-1.0));
        assert_eq!(iter.count(), 1);
    }
    {
        let iter = line.polygonize(SegmentPolygonizer::InnerSegments(3));
        assert_eq!(iter.count(), 4);

        let iter = line.polygonize(SegmentPolygonizer::InnerSegments(0));
        assert_eq!(iter.count(), 1);
    }
}

#[test]
fn test_polygonize_segment_length() {
    let line: Segment = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
        .unwrap()
        .into();

    let iter = line.polygonize(SegmentPolygonizer::MaximumSegmentLength(0.3));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    assert_eq!(chain.num_segments(), 4);

    for segments in chain.as_slices().0.windows(2) {
        approx::assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }

    let iter = line.polygonize(SegmentPolygonizer::InnerSegments(8));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    for segments in chain.as_slices().0.windows(2) {
        approx::assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }

    let iter = line.polygonize(SegmentPolygonizer::MaximumAngle(0.1));
    let pts: Vec<_> = iter.collect();
    let chain = SegmentChain::from_points(&pts);
    for segments in chain.as_slices().0.windows(2) {
        approx::assert_abs_diff_eq!(segments[0].length(), segments[1].length(), epsilon = 1e-14);
    }
}

#[test]
fn test_euclidian_distance_to_point() {
    let line = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    assert_eq!(1.0, line.euclidian_distance_to_point([0.0, 1.0]));
    approx::assert_abs_diff_eq!(1.0, line.euclidian_distance_to_point([0.5, 1.0]));
    approx::assert_abs_diff_eq!(1.25f64.sqrt(), line.euclidian_distance_to_point([1.5, 1.0]));

    approx::assert_abs_diff_eq!(0.0, line.euclidian_distance_to_point([0.5, 0.0]));
}

#[test]
fn test_segment_point() {
    {
        let line =
            LineSegment::new([0.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert_eq!(line.segment_point(0.5), [0.5, 0.5]);
        assert_eq!(line.segment_point(-0.5), [0.0, 0.0]);
        assert_eq!(line.segment_point(1.5), [1.0, 1.0]);
        assert_eq!(line.segment_point(0.2), [0.2, 0.2]);
    }
    {
        let line =
            LineSegment::new([-2.0, 0.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert_eq!(line.segment_point(0.5), [0.0, 0.0]);
        assert_eq!(line.segment_point(-0.5), [-2.0, 0.0]);
        assert_eq!(line.segment_point(1.5), [2.0, 0.0]);
        assert_eq!(line.segment_point(0.2), [-1.2, 0.0]);
    }
}

#[test]
fn test_angle_infinite() {
    let line = LineSegment::new(
        [NEG_INFINITY, NEG_INFINITY],
        [INFINITY, INFINITY],
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();
    approx::assert_abs_diff_eq!(line.angle(), 0.25 * PI);

    let line = LineSegment::new(
        [INFINITY, INFINITY],
        [NEG_INFINITY, NEG_INFINITY],
        DEFAULT_EPSILON,
        DEFAULT_MAX_ULPS,
    )
    .unwrap();
    approx::assert_abs_diff_eq!(line.angle(), -0.75 * PI);
}

#[test]
fn test_contains_point() {
    {
        let line: Segment =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        assert!(line.contains_point([0.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS));
    }
    {
        let line: Segment =
            LineSegment::new([0.5, 1.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        assert!(line.contains_point(
            [0.9000000000000001, 0.19999999999999984],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS
        ));
    }
}

#[test]
fn test_line_segments_intersection() {
    {
        let line = Line::from_point_angle([0.11, 0.11], 0.0);
        let line_segment = LineSegment::new([0.1, 0.9], [0.1, 0.1], 0.0, 0).unwrap();
        assert_eq!(
            line.intersections_line_segment(&line_segment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .len(),
            1
        );
        assert_eq!(
            line_segment
                .intersections_line(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .len(),
            1
        );
    }
    {
        // These two segments intersect
        let line1 =
            LineSegment::new([0.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([1.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_line_segment(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([0.5, 0.5]), intersections);

        // This segment doesn't intersect with line1 because it is parallel to it
        let line3 =
            LineSegment::new([0.0, 1.0], [1.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_line_segment(&line3, 0.0, 0);
        assert_eq!(intersections.len(), 0);

        // This segment doesn't intersect with line1 even though it is not parallel
        let line4 =
            LineSegment::new([0.0, 1.0], [1.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersection = line1.intersections_line_segment(&line4, 0.0, 0);
        assert_eq!(intersection.len(), 0);
    }
    {
        // Segment intersects other segment at the end
        let line1 =
            LineSegment::new([0.0, 0.0], [0.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([0.5, -0.5], [0.5, 1.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_line_segment(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([0.5, 0.0]), intersections);
    }
    {
        // One segment which is included in the other segment
        let line1 =
            LineSegment::new([0.0, 1.0], [0.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([0.0, 0.0], [0.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_line_segment(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 2);
        approx::assert_ulps_eq!(
            PrimitiveIntersections::Two([[0.0, 1.0], [0.0, 2.0]]),
            intersections
        );
    }
    {
        // Overlapping segments
        let line1 =
            LineSegment::new([0.0, 0.0], [0.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([0.0, 1.0], [0.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersection = line1.intersections_line_segment(&line2, 0.0, 0);
        assert_eq!(intersection.len(), 2);
        approx::assert_ulps_eq!(
            PrimitiveIntersections::Two([[0.0, 1.0], [0.0, 2.0]]),
            intersection
        );
    }
    {
        // Touching segments
        let line1 =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([0.0, 1.0], [0.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_line_segment(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([0.0, 1.0]), intersections);
    }
    {
        // One line vertical, other horizontal, one intersection occurs
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2: Segment =
            LineSegment::new([0.5, -0.5], [0.5, 0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let intersections = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([0.5, 0.0]), intersections);
    }
    {
        // One line vertical, other horizontal, no intersection occurs
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.5, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2: Segment =
            LineSegment::new([0.5, -1.5], [0.5, -0.5], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 0);
    }
    {
        // Intersection of two identical vertical lines
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2 = line_1.clone();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 2);
    }
    {
        // Intersection of two parallel vertical lines
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2: Segment =
            LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 0);
    }
    {
        // Intersection of two vertical lines stacked on top of each other
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2: Segment =
            LineSegment::new([0.0, 1.0], [0.0, 2.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 1);
    }
    {
        // Intersection of two vertical lines stacked on top of each other with space in between
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [0.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2: Segment =
            LineSegment::new([0.0, 2.0], [0.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 0);
    }
    {
        // Intersection of two identical horizontal lines
        let line_1: Segment =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line_2 = line_1.clone();

        let results = line_1.intersections_primitive(&line_2, 0.0, 0);
        assert_eq!(results.len(), 2);
    }
    {
        // Regression test for a bug detected with the shapes method of IsSlot
        let line_1: Segment = LineSegment::new(
            [0.003999999999999997, 0.0027499999999999933],
            [0.004, 0.014749999999999997],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap()
        .into();
        let line_2: Segment = LineSegment::new(
            [-0.008, 0.009078996041164352],
            [0.008, 0.009078996041164352],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap()
        .into();
        let intersections = line_1.intersections_primitive(&line_2, 0.0, 0);
        approx::assert_ulps_eq!(
            PrimitiveIntersections::One([0.004, 0.009078996041164352]),
            intersections
        );
    }
    {
        // Parallel lines
        let line_1 =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line_2 =
            LineSegment::new([0.0, 1.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert_eq!(line_1.intersections_primitive(&line_2, 0.0, 0).len(), 0);
    }
    {
        // Parallel lines
        let line_1 =
            LineSegment::new([0.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line_2 =
            LineSegment::new([2.0, 3.0], [3.0, 4.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert_eq!(line_1.intersections_primitive(&line_2, 0.0, 0).len(), 0);
    }
    {
        // Parallel lines
        let line_1 =
            LineSegment::new([0.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line_2 =
            LineSegment::new([3.0, 4.0], [2.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        assert_eq!(line_1.intersections_primitive(&line_2, 0.0, 0).len(), 0);
    }
    {
        let line_1 =
            LineSegment::new([0.0, 1.0], [1.0, -9.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line_2 =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [0.09999999999999998, 0.0]);
    }
    {
        let line_1 =
            LineSegment::new([0.0, 1.0], [1.0, -9.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line_2 =
            LineSegment::new([0.0, 1.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [0.0, 1.0]);
    }
    {
        let line1 =
            LineSegment::new([1.0, 0.0], [1.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([-10.0, 1.0], [10.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_primitive(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([1.0, 1.0]), intersections);
    }
    {
        let line1 =
            LineSegment::new([1.0, 1.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let line2 =
            LineSegment::new([-10.0, 1.0], [10.0, 1.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
        let intersections = line1.intersections_primitive(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 1);
        approx::assert_ulps_eq!(PrimitiveIntersections::One([1.0, 1.0]), intersections);
    }
    {
        // Regression test from the Delaunay triangulation
        let line1: Segment =
            LineSegment::new([6.0, 2.0], [7.0, -3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let line2: Segment =
            LineSegment::new([0.0, 0.0], [1.0, 3.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();
        let intersections = line1.intersections_primitive(&line2, 0.0, 0);
        assert_eq!(intersections.len(), 0);
    }
    {
        // Based on JTS test `testCentralEndpointHeuristicFailure`
        // > Following cases were failures when using the CentralEndpointIntersector heuristic.
        // > This is because one segment lies at a significant angle to the other,
        // > with only one endpoint is close to the other segment.
        // > The CE heuristic chose the wrong endpoint to return.
        // > The fix is to use a new heuristic which out of the 4 endpoints
        // > chooses the one which is closest to the other segment.
        // > This works in all known failure cases.
        let line_1 = LineSegment::new(
            [163.81867067, -211.31840378],
            [165.9174252, -214.1665075],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [2.84139601, -57.95412726],
            [469.59990601, -502.63851732],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [163.81867067, -211.31840378]);
    }
    {
        // Based on JTS test `testCentralEndpointHeuristicFailure2`
        // > Test from Tomas Fa - JTS list 6/13/2012
        // >
        // > Fails using original JTS DeVillers determine orientation test.
        // > Succeeds using DD and Shewchuk orientation
        let line_1 = LineSegment::new(
            [-58.00593335955, -1.43739086465],
            [-513.86101637525, -457.29247388035],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [-215.22279674875, -158.65425425385],
            [-218.1208801283, -160.68343590235],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [-215.22279674875, -158.65425425385]);
    }
    {
        // Based on JTS test `testTomasFa_1`
        // > Test from Tomas Fa - JTS list 6/13/2012
        // >
        // > Fails using original JTS DeVillers determine orientation test.
        // > Succeeds using DD and Shewchuk orientation
        let line_1 = LineSegment::new(
            [-42.0, 163.2],
            [21.2, 265.2],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [-26.2, 188.7],
            [37.0, 290.7],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        assert_eq!(line_1.intersections_primitive(&line_2, 0.0, 0).len(), 0);
    }
    {
        // Based on JTS test `testTomasFa_2`
        //
        // > Test from Tomas Fa - JTS list 6/13/2012
        // >
        // > Fails using original JTS DeVillers determine orientation test.
        let line_1 = LineSegment::new(
            [-5.9, 163.1],
            [76.1, 250.7],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [14.6, 185.0],
            [96.6, 272.6],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        assert_eq!(line_1.intersections_primitive(&line_2, 0.0, 0).len(), 0);
    }
    {
        // Based on JTS test `testLeduc_1`
        //
        // > Test involving two non-almost-parallel lines.
        // > Does not seem to cause problems with basic line intersection algorithm.
        let line_1 = LineSegment::new(
            [305690.0434123494, 254176.46578338774],
            [305601.9999843455, 254243.19999846347],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [305689.6153764265, 254177.33102743194],
            [305692.4999844298, 254171.4999983967],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [305690.0434123494, 254176.46578338774]);
    }
    {
        // Based on JTS test `testGEOS_1()`
        //
        // > Test from strk which is bad in GEOS (2009-04-14).
        let line_1 = LineSegment::new(
            [588750.7429703881, 4518950.493668233],
            [588748.2060409798, 4518933.9452804085],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [588745.824857241, 4518940.742239175],
            [588748.2060437313, 4518933.9452791475],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [588748.2060416829, 4518933.945284994]);
    }
    {
        // Based on JTS test `testGEOS_2()`
        //
        // > Test from strk which is bad in GEOS (2009-04-14).
        let line_1 = LineSegment::new(
            [588743.626135934, 4518924.610969561],
            [588732.2822865889, 4518925.4314047815],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [588739.1191384895, 4518927.235700594],
            [588731.7854614238, 4518924.578370095],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [588733.8306132929, 4518925.319423238]);
    }
    {
        // Based on JTS test `testDaveSkeaCase()`
        //
        // > This used to be a failure case (exception), but apparently works now.
        // > Possibly normalization has fixed this?
        let line_1 = LineSegment::new(
            [2089426.5233462777, 1180182.387733969],
            [2085646.6891757075, 1195618.7333999649],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [1889281.8148903656, 1997547.0560044837],
            [2259977.3672236, 483675.17050843034],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [2087536.6062609926, 1187900.560566967]);
    }
    {
        // Based on JTS test `testCmp5CaseWKT()`
        //
        // > Outside envelope using HCoordinate method.
        let line_1 = LineSegment::new(
            [4348433.262114629, 5552595.478385733],
            [4348440.849387404, 5552599.272022122],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let line_2 = LineSegment::new(
            [4348433.26211463, 5552595.47838573],
            [4348440.8493874, 5552599.27202212],
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        )
        .unwrap();
        let actual = line_1.intersections_primitive(&line_2, 0.0, 0)[0];
        assert_eq!(actual, [4348440.8493874, 5552599.27202212]);
    }
}

#[test]
fn self_intersection() {
    let ls = LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    // Self-intersection
    assert_eq!(
        ls.intersections_primitive(&ls, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Zero
    );

    // Intersections with equal primitive
    let ls_cloned = ls.clone();
    assert_eq!(
        ls.intersections_primitive(&ls_cloned, DEFAULT_EPSILON, DEFAULT_MAX_ULPS),
        PrimitiveIntersections::Two([[0.0, 0.0], [1.0, 0.0]])
    );
}

#[test]
fn test_intersection_line_line_segment() {
    let ls = LineSegment::new([0.0, 1.0], [0.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    let line = Line::from_point_angle([0.5, 0.5], 0.0);

    let intersection = ls.intersections_primitive(&line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    approx::assert_abs_diff_eq!(intersection, PrimitiveIntersections::One([0.0, 0.5]));
}
