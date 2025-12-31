use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use approx::assert_abs_diff_eq;
use planar_geo::{contour::ArrowHeadSize, prelude::Contour};

#[test]
fn test_arrow_from_tail_length_angle() {
    // Horizontal arrow
    {
        let tail = [0.0, 0.0];
        let length = 4.0;
        let angle = 0.0;
        let arrow =
            Contour::arrow_from_tail_length_angle(tail, length, angle, ArrowHeadSize::Height(1.0))
                .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 6);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [3.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[2], [3.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[3], [4.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[4], [3.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[5], [3.0, 0.0], epsilon = 1e-5);
    }

    {
        let tail = [0.0, 0.0];
        let length = 4.0;
        let angle = PI;
        let arrow =
            Contour::arrow_from_tail_length_angle(tail, length, angle, ArrowHeadSize::Height(1.0))
                .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 6);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [-3.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[2], [-3.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[3], [-4.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[4], [-3.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[5], [-3.0, 0.0], epsilon = 1e-5);
    }

    // Vertical arrow
    {
        let tail = [0.0, 0.0];
        let length = 4.0;
        let angle = FRAC_PI_2;
        let arrow =
            Contour::arrow_from_tail_length_angle(tail, length, angle, ArrowHeadSize::Height(1.0))
                .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 6);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [0.0, 3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[2], [-1.0 / 3.0f64.sqrt(), 3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[3], [0.0, 4.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[4], [1.0 / 3.0f64.sqrt(), 3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[5], [0.0, 3.0], epsilon = 1e-5);
    }

    {
        let tail = [0.0, 0.0];
        let length = 4.0;
        let angle = 3.0 * FRAC_PI_2;
        let arrow =
            Contour::arrow_from_tail_length_angle(tail, length, angle, ArrowHeadSize::Height(1.0))
                .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 6);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [0.0, -3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[2], [1.0 / 3.0f64.sqrt(), -3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[3], [0.0, -4.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[4], [-1.0 / 3.0f64.sqrt(), -3.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[5], [0.0, -3.0], epsilon = 1e-5);
    }

    // Diagonal arrow
    {
        let tail = [0.0, 0.0];
        let length = 2.0 * 2.0f64.sqrt();
        let angle = FRAC_PI_4;
        let arrow = Contour::arrow_from_tail_length_angle(
            tail,
            length,
            angle,
            ArrowHeadSize::Height(2.0f64.sqrt()),
        )
        .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 6);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [1.0, 1.0], epsilon = 1e-5);
        assert_abs_diff_eq!(
            verts[2],
            [1.0 - 1.0 / 3.0f64.sqrt(), 1.0 + 1.0 / 3.0f64.sqrt()],
            epsilon = 1e-5
        );
        assert_abs_diff_eq!(verts[3], [2.0, 2.0], epsilon = 1e-5);
        assert_abs_diff_eq!(
            verts[4],
            [1.0 + 1.0 / 3.0f64.sqrt(), 1.0 - 1.0 / 3.0f64.sqrt()],
            epsilon = 1e-5
        );
        assert_abs_diff_eq!(verts[5], [1.0, 1.0], epsilon = 1e-5);
    }

    // "Tailless" arrow
    {
        let tail = [0.0, 0.0];
        let length = 1.0;
        let angle = 0.0;
        let arrow =
            Contour::arrow_from_tail_length_angle(tail, length, angle, ArrowHeadSize::Height(1.0))
                .unwrap();
        let verts: Vec<[f64; 2]> = arrow.points().collect();

        assert_eq!(verts.len(), 4);
        assert_abs_diff_eq!(verts[0], [0.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[1], [0.0, 1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[2], [1.0, 0.0], epsilon = 1e-5);
        assert_abs_diff_eq!(verts[3], [0.0, -1.0 / 3.0f64.sqrt()], epsilon = 1e-5);
    }
}
