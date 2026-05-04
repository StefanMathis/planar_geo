use planar_geo::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn test_covers() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;
    {
        let s1: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m)
            .unwrap()
            .into();
        let quarter =
            ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], 0.5 * PI, e, m).unwrap();
        assert!(s1.covers_arc_segment(&quarter, e, m));

        let s2: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m)
            .unwrap()
            .into();
        assert!(!s2.covers_arc_segment(&quarter, e, m));
    }
    {
        let s1: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m)
            .unwrap()
            .into();
        let ls_start_to_middle = LineSegment::new([0.0, 0.0], [0.5, 0.5], e, m).unwrap();
        assert!(s1.covers_line_segment(&ls_start_to_middle, e, m));

        let s2: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m)
            .unwrap()
            .into();
        assert!(!s2.covers_line_segment(&ls_start_to_middle, e, m));
    }
}

#[test]
fn test_fillet_chain() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;
    {
        let mut iter = Segment::fillet_chain(
            &[
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 0.5],
                [0.5, 0.5],
                [0.5, 1.0],
                [1.0, 1.0],
            ],
            &[0.5, 0.5, 2.0],
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [0.5, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([1.0, 0.5], [1.0, 1.0], -FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.5, 1.0], [1.0, 1.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        assert_eq!(iter.next(), None);
        assert_eq!(iter.next(), None);
        assert_eq!(iter.next(), None);
        assert_eq!(iter.next(), None);
    }
    {
        // No radii
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 0.5]], &[]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [1.0, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.0], [1.0, 0.5], e, m)
                    .unwrap()
                    .into()
            )
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Small radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[0.5]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [0.5, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.5], [1.0, 1.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Negative radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[-0.5]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [1.0, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.0], [1.0, 1.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Very large radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[20.0]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            )
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Small and large radius
        let mut iter = Segment::fillet_chain(
            &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            &[0.5, 20.0],
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [0.5, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([1.0, 0.5], [0.5, 0.5], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.5, 1.0], [0.0, 1.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Small and large radius
        let mut iter = Segment::fillet_chain(
            &[[1.0, 0.0], [0.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            &[0.5, 20.0],
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.0], [0.5, 0.0], e, m)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], -FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.0, 0.5], [0.5, 0.5], -FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.5, 1.0], [1.0, 1.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Two points -> radii are ignored
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0]], &[0.5, 20.0]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [1.0, 0.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Empty iterator
        let mut iter = Segment::fillet_chain(&[[1.0, 0.0]], &[0.5, 20.0]);
        assert_eq!(iter.next(), None);
    }
    {
        // Doctest example
        let mut iter = Segment::fillet_chain(
            &[
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 0.5],
                [0.5, 0.5],
                [0.5, 1.0],
                [0.0, 1.0],
            ],
            &[0.5, 0.0, 0.25, 2.0],
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [0.5, 0.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.5], [0.75, 0.5], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.75, 0.5], [0.75, 0.75], -FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.75], [0.25, 0.75], FRAC_PI_2, e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.25, 1.0], [0.0, 1.0], e, m)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
}
