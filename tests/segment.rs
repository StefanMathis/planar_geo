use planar_geo::prelude::*;
use std::f64::consts::PI;

#[test]
pub fn test_contains() {
    {
        let e = DEFAULT_EPSILON;
        let m = DEFAULT_MAX_ULPS;

        let s1: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m)
            .unwrap()
            .into();
        let quarter =
            ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], 0.5 * PI, e, m).unwrap();
        assert!(s1.contains_arc_segment(&quarter, e, m));

        let s2: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m)
            .unwrap()
            .into();
        assert!(!s2.contains_arc_segment(&quarter, e, m));
    }
    {
        let e = DEFAULT_EPSILON;
        let m = DEFAULT_MAX_ULPS;

        let s1: Segment = LineSegment::new([0.0, 0.0], [1.0, 1.0], e, m)
            .unwrap()
            .into();
        let ls_start_to_middle = LineSegment::new([0.0, 0.0], [0.5, 0.5], e, m).unwrap();
        assert!(s1.contains_line_segment(&ls_start_to_middle, e, m));

        let s2: Segment = ArcSegment::from_start_center_angle([0.0, 0.0], [0.0, 1.0], PI, e, m)
            .unwrap()
            .into();
        assert!(!s2.contains_line_segment(&ls_start_to_middle, e, m));
    }
}
