use std::f64::consts::{FRAC_PI_2, PI};

use planar_geo::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS, prelude::*};

#[test]
fn test_bounding_line_segment() {
    {
        let segment_chain: SegmentChain =
            LineSegment::new([0.0, 0.0], [1.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
                .unwrap()
                .into();

        let bb = BoundingBox::from(&segment_chain);
        assert_eq!(bb.xmin(), 0.0);
        assert_eq!(bb.xmax(), 1.0);
        assert_eq!(bb.ymin(), 0.0);
        assert_eq!(bb.ymax(), 0.0);
    }
}

#[test]
fn test_bounding_box_arc() {
    {
        let segment_chain: SegmentChain = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            0.5,
            0.0,
            FRAC_PI_2,
            0.0,
            0,
        )
        .unwrap()
        .into();

        let bb = BoundingBox::from(&segment_chain);
        approx::assert_abs_diff_eq!(bb.xmin(), 0.0);
        approx::assert_abs_diff_eq!(bb.xmax(), 0.5);
        approx::assert_abs_diff_eq!(bb.ymin(), 0.0);
        approx::assert_abs_diff_eq!(bb.ymax(), 0.5);
    }
    {
        let segment_chain: SegmentChain = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            0.5,
            -FRAC_PI_2,
            -FRAC_PI_2,
            0.0,
            0,
        )
        .unwrap()
        .into();

        let bb = BoundingBox::from(&segment_chain);
        approx::assert_abs_diff_eq!(bb.xmin(), -0.5);
        approx::assert_abs_diff_eq!(bb.xmax(), 0.0);
        approx::assert_abs_diff_eq!(bb.ymin(), -0.5);
        approx::assert_abs_diff_eq!(bb.ymax(), 0.0);
    }
    {
        let segment_chain: SegmentChain = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.25 * PI,
            0.5 * PI,
            0.0,
            0,
        )
        .unwrap()
        .into();

        let bb = BoundingBox::from(&segment_chain);
        approx::assert_abs_diff_eq!(bb.xmin(), -2.0f64.sqrt(),);
        approx::assert_abs_diff_eq!(bb.xmax(), 2.0f64.sqrt(),);
        approx::assert_abs_diff_eq!(bb.ymin(), 2.0f64.sqrt(),);
        approx::assert_abs_diff_eq!(bb.ymax(), 2.0);
    }
    {
        let segment_chain: SegmentChain = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            2.0,
            0.25 * PI,
            0.1 * PI,
            0.0,
            0,
        )
        .unwrap()
        .into();

        let bb = BoundingBox::from(&segment_chain);
        approx::assert_abs_diff_eq!(bb.xmin(), 0.90798, epsilon = 1e-4);
        approx::assert_abs_diff_eq!(bb.xmax(), 2.0f64.sqrt(),);
        approx::assert_abs_diff_eq!(bb.ymin(), 2.0f64.sqrt(),);
        approx::assert_abs_diff_eq!(bb.ymax(), 1.78201, epsilon = 1e-4);
    }
}

#[test]
fn test_bounding_two_box_shapes() {
    let vertices = [[10e-3, 0.0], [10e-3, 5e-3], [-10e-3, 5e-3], [-10e-3, 0.0]];
    let chain = SegmentChain::from_points(vertices.as_slice());
    let shape_lower = Shape::new(vec![Contour::new(chain)]).unwrap();

    let vertices = vec![
        [10e-3, 5e-3],
        [10e-3, 10e-3],
        [-10e-3, 10e-3],
        [-10e-3, 5e-3],
    ];
    let chain = SegmentChain::from_points(vertices.as_slice());
    let shape_upper = Shape::new(vec![Contour::new(chain)]).unwrap();

    let shapes = vec![shape_lower, shape_upper];
    let bb = BoundingBox::from_bounded_entities(shapes.iter()).unwrap();
    approx::assert_abs_diff_eq!(bb.xmin(), -10e-3);
    approx::assert_abs_diff_eq!(bb.xmax(), 10e-3);
    approx::assert_abs_diff_eq!(bb.ymin(), 0.0);
    approx::assert_abs_diff_eq!(bb.ymax(), 10e-3);
}

#[test]
fn test_intersects() {
    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(-0.5, 0.5, -0.5, 1.5);
    assert!(bb1.intersects(&bb2));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(-0.5, 0.5, -0.5, 1.5);
    assert!(bb1.intersects(&bb2));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(0.2, 0.8, std::f64::NEG_INFINITY, 0.5);
    assert!(bb1.intersects(&bb2));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(0.2, 0.8, 0.5, std::f64::INFINITY);
    assert!(bb1.intersects(&bb2));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(0.2, 0.8, std::f64::NEG_INFINITY, std::f64::INFINITY);
    assert!(bb1.intersects(&bb2));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(std::f64::NEG_INFINITY, std::f64::INFINITY, 0.2, 0.8);
    assert!(bb1.intersects(&bb2));
}

#[test]
fn test_contains() {
    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(-0.5, 1.5, -0.5, 1.5);
    assert!(bb2.approx_contains(&bb1, 0.0, 0));
    assert!(!bb1.approx_contains(&bb2, 0.0, 0));

    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2 = BoundingBox::new(0.2, 1.0, 0.2, 0.8);
    assert!(bb1.approx_contains(&bb2, 0.0, 0));
}

#[test]
fn test_touches() {
    // bb1 touches bb2
    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2: BoundingBox = BoundingBox::new(1.0, 2.0, 0.0, 1.0);
    assert!(bb2.touches(&bb1));
    assert!(bb1.touches(&bb2));

    // bb1 intersects bb2
    let bb1 = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let bb2: BoundingBox = BoundingBox::new(0.8, 2.0, 0.0, 1.0);
    assert!(!bb2.touches(&bb1));
    assert!(!bb1.touches(&bb2));
}
