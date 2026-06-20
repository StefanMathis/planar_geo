use planar_geo::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn test_covers() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_RELATIVE;
    {
        let s1: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI)
            .unwrap()
            .into();
        let quarter =
            ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], 0.5 * PI).unwrap();
        assert!(s1.covers_arc_segment(&quarter, e, m));

        let s2: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap().into();
        assert!(!s2.covers_arc_segment(&quarter, e, m));
    }
    {
        let s1: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap().into();
        let ls_start_to_middle = LineSegment::new([0.0, 0.0], [0.5, 0.5]).unwrap();
        assert!(s1.covers_line_segment(&ls_start_to_middle, e, m));

        let s2: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI)
            .unwrap()
            .into();
        assert!(!s2.covers_line_segment(&ls_start_to_middle, e, m));
    }
}

#[test]
fn test_fillet_chain() {
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
            &[0.0, 0.5, 2.0],
        );

        let gap_check = iter.clone();
        let gap_check_skip_one = gap_check.clone();

        // Check if there is a gap between the segments
        for (prev, next) in gap_check.zip(gap_check_skip_one.skip(1)) {
            assert_eq!(prev.stop(), next.start());
        }

        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([1.0, 0.0], [0.5, 0.0], FRAC_PI_2)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.5, 0.5], [0.5, 1.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.5, 1.0], [1.0, 1.0]).unwrap().into())
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
            Some(LineSegment::new([0.0, 0.0], [0.5, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([1.0, 0.5], [0.5, 0.5], FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.5, 1.0], [0.0, 1.0]).unwrap().into()),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
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

        let gap_check = iter.clone();
        let gap_check_skip_one = gap_check.clone();

        // Check if there is a gap between the segments
        for (prev, next) in gap_check.zip(gap_check_skip_one.skip(1)) {
            assert_eq!(prev.stop(), next.start());
        }

        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [0.4999999999999999, 0.0])
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0000000000000002, 0.5000000000000001], [1.0, 0.5])
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.5], [0.9999999999999997, 0.5])
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([1.0, 0.5], [1.0, 1.0], -FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-6
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.4999999999999998, 1.0], [0.5, 1.0])
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.5, 1.0], [1.0, 1.0]).unwrap().into())
        );

        // Assert that the iterator is exhausted
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
            Some(LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([1.0, 0.0], [1.0, 0.5]).unwrap().into())
        );
        assert_eq!(iter.next(), None);
    }
    {
        println!("\n");

        // Small radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[0.5]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.0, 0.0], [0.5, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([1.0, 0.5], [1.0, 1.0]).unwrap().into()),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Negative radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[-0.5]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([1.0, 0.0], [1.0, 1.0]).unwrap().into())
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Very large radius
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], &[20.0]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([0.0, 0.0], [-1.8369701987210297e-16, 0.0])
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], FRAC_PI_2)
                    .unwrap()
                    .into()
            )
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0, 0.9999999999999998], [1.0, 1.0])
                    .unwrap()
                    .into()
            )
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
            Some(LineSegment::new([1.0, 0.0], [0.5, 0.0]).unwrap().into())
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], -FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        // Glue segment
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new(
                    [-5.551115123125783e-17, 0.5000000000000001],
                    [-1.1102230246251565e-16, 0.5]
                )
                .unwrap()
                .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.0, 0.5], [0.5, 0.5], -FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.5, 1.0], [1.0, 1.0],).unwrap().into()),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
    {
        // Two points -> radii are ignored
        let mut iter = Segment::fillet_chain(&[[0.0, 0.0], [1.0, 0.0]], &[0.5, 20.0]);
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap().into()),
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
            Some(LineSegment::new([0.0, 0.0], [0.5, 0.0]).unwrap().into()),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.0], [0.5, 0.5], FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new([1.0000000000000002, 0.5000000000000001], [1.0, 0.5])
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([1.0, 0.5], [0.75, 0.5]).unwrap().into()),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.75, 0.5], [0.75, 0.75], -FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                LineSegment::new(
                    [0.5000000000000006, 0.75],
                    [0.5000000000000006, 0.7499999999999999]
                )
                .unwrap()
                .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(
                ArcSegment::from_start_center_angle([0.5, 0.75], [0.25, 0.75], FRAC_PI_2)
                    .unwrap()
                    .into()
            ),
            epsilon = 1e-8,
        );
        approx::assert_abs_diff_eq!(
            iter.next(),
            Some(LineSegment::new([0.25, 1.0], [0.0, 1.0]).unwrap().into()),
            epsilon = 1e-8,
        );
        assert_eq!(iter.next(), None);
    }
}
