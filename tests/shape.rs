use planar_geo::{error::ShapeConstructorError, prelude::*};
use std::f64::consts::FRAC_PI_2;

#[test]
fn test_new() {
    // Shape with fillets
    {
        let e = planar_geo::DEFAULT_EPSILON;
        let m = planar_geo::DEFAULT_MAX_ULPS;

        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([0.0, 100.0], [0.0, 0.0], [100.0, 0.0], 50.0, e, m)
                .unwrap()
                .into(),
        );
        polysegment.extend_back([100.0, 0.0]);
        polysegment.push_back(
            ArcSegment::fillet([100.0, 0.0], [0.0, 100.0], [0.0, 0.0], 10.0, e, m)
                .unwrap()
                .into(),
        );
        assert!(Shape::new(vec![polysegment.into()]).is_ok());
    }

    // Shape without any hole
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c = Contour::new(Polysegment::from_points(vertices));
        assert!(Shape::new(vec![c]).is_ok());
    }

    // Shape with a single hole
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c1 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
        let c2 = Contour::new(Polysegment::from_points(vertices));

        assert!(Shape::new(vec![c1, c2]).is_ok());
    }

    // Shape with two holes
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c1 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
        let c2 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
        let c3 = Contour::new(Polysegment::from_points(vertices));

        assert!(Shape::new(vec![c1, c2, c3]).is_ok());
    }

    // Fails because the input contour vector is empty
    {
        let err = Shape::new(Vec::new()).unwrap_err();
        assert_eq!(err, ShapeConstructorError::EmptyVec);
    }

    // Fails because one of the contours is empty
    {
        let err = Shape::new(vec![Contour::new(Polysegment::new())]).unwrap_err();
        match err {
            ShapeConstructorError::EmptyContour { input: _, idx } => assert_eq!(idx, 0),
            _ => unreachable!(),
        }
    }
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c = Contour::new(Polysegment::from_points(vertices));

        let err = Shape::new(vec![c, Contour::new(Polysegment::new())]).unwrap_err();
        match err {
            ShapeConstructorError::EmptyContour { input: _, idx } => assert_eq!(idx, 1),
            _ => unreachable!(),
        }
    }
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c = Contour::new(Polysegment::from_points(vertices));

        let err = Shape::new(vec![Contour::new(Polysegment::new()), c]).unwrap_err();
        match err {
            ShapeConstructorError::EmptyContour { input: _, idx } => assert_eq!(idx, 0),
            _ => unreachable!(),
        }
    }

    // Fails because the contour intersects itself
    {
        let vertices = &[[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]];
        let c = Contour::new(Polysegment::from_points(vertices));

        let err = Shape::new(vec![c]).unwrap_err();
        match err {
            ShapeConstructorError::Intersection {
                input: _,
                intersection,
            } => assert_eq!(
                intersection,
                Intersection {
                    point: [0.5, 0.5],
                    left: SegmentKey {
                        contour_idx: 0,
                        segment_idx: 0
                    },
                    right: SegmentKey {
                        contour_idx: 0,
                        segment_idx: 2
                    }
                }
            ),
            _ => unreachable!(),
        }
    }

    // Fails because the two contours intersect
    {
        let vertices = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let c1 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[-1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [-1.0, 3.0]];
        let c2 = Contour::new(Polysegment::from_points(vertices));

        let err = Shape::new(vec![c1, c2]).unwrap_err();
        match err {
            ShapeConstructorError::HoleOutsideContour { input: _, idx: _ } => (),
            _ => unreachable!(),
        }
    }

    // Fails because the second hole is within the first hole
    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c1 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
        let c2 = Contour::new(Polysegment::from_points(vertices));

        let vertices = &[[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8]];
        let c3 = Contour::new(Polysegment::from_points(vertices));

        let err = Shape::new(vec![c1, c2, c3]).unwrap_err();
        match err {
            ShapeConstructorError::HoleInsideHole {
                input: _,
                outer_hole_idx: _,
                inner_hole_idx: _,
            } => (),
            _ => unreachable!(),
        }
    }
}

#[test]
fn test_add_hole() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c = Contour::new(Polysegment::from_points(vertices));

    let mut shape = Shape::new(vec![c]).unwrap();
    assert_eq!(shape.holes().len(), 0);

    // Adding works
    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_ok());
    assert_eq!(shape.holes().len(), 1);

    // Adding fails (intersection)
    let vertices = &[[-1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [-1.0, 3.0]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_err());
    assert_eq!(shape.holes().len(), 1);

    // Adding works
    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_ok());
    assert_eq!(shape.holes().len(), 2);

    // Adding fails (new hole is within first hole)
    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));
    assert!(shape.add_hole(hole).is_err());
    assert_eq!(shape.holes().len(), 2);
}

#[test]
fn test_remove_hole() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.4, 0.1], [0.4, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.6, 0.1], [0.9, 0.1], [0.9, 0.9], [0.6, 0.9]];
    let c3 = Contour::new(Polysegment::from_points(vertices));

    let mut shape = Shape::new(vec![c1, c2, c3]).unwrap();

    assert_eq!(shape.holes().len(), 2);
    assert!(shape.remove_hole(2).is_none());

    assert!(shape.remove_hole(1).is_some());
    assert_eq!(shape.holes().len(), 1);

    assert!(shape.remove_hole(0).is_some());
    assert_eq!(shape.holes().len(), 0);

    assert!(shape.remove_hole(0).is_none());
    assert_eq!(shape.holes().len(), 0);
}

#[test]
fn test_translate() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(&vertices));

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(&vertices));

    let mut shape = Shape::new(vec![c1, c2]).unwrap();
    shape.translate([1.0, 0.0]);

    let vertices: Vec<[f64; 2]> = shape.contour().points().collect();
    assert_eq!(vertices[0], [1.0, 0.0]);
    assert_eq!(vertices[1], [2.0, 0.0]);
    assert_eq!(vertices[2], [2.0, 1.0]);
    assert_eq!(vertices[3], [1.0, 1.0]);
}

#[test]
fn test_rotation() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(&vertices));

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(&vertices));

    let mut shape = Shape::new(vec![c1, c2]).unwrap();
    shape.rotate([0.0, 0.0], FRAC_PI_2);

    let vertices: Vec<[f64; 2]> = shape.contour().points().collect();
    approx::assert_abs_diff_eq!(vertices[0], [0.0, 0.0]);
    approx::assert_abs_diff_eq!(vertices[1], [0.0, 1.0]);
    approx::assert_abs_diff_eq!(vertices[2], [-1.0, 1.0]);
    approx::assert_abs_diff_eq!(vertices[3], [-1.0, 0.0]);
}

#[test]
fn test_rectangle_with_hole() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(&vertices));

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(&vertices));

    let shape = Shape::new(vec![c1, c2]).unwrap();

    assert!(shape.contains_point([0.0, 0.0], DEFAULT_EPSILON, 0));
    assert!(shape.contains_point([0.1, 0.1], DEFAULT_EPSILON, 0));
    assert!(!shape.contains_point([0.11, 0.11], DEFAULT_EPSILON, 0));
    assert!(shape.contains_point([0.11, 0.11], 0.2, 0));
    assert!(!shape.contains_point([0.5, 0.5], DEFAULT_EPSILON, 0));
    assert!(!shape.contains_point([0.0, -0.05], DEFAULT_EPSILON, 0));
    assert!(shape.contains_point([0.0, -0.05], 0.05, 0));
    assert!(!shape.contains_point([0.0, -0.05], 0.02, 0));
}

#[test]
fn test_intersection_with_polysegment() {
    let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let contour = Contour::new(Polysegment::from_points(vertices));

    let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let hole = Contour::new(Polysegment::from_points(vertices));

    let shape = Shape::new(vec![contour, hole]).unwrap();

    let polysegment = Polysegment::from_points(&[
        [-1.0, 1.0],
        [-1.0, 0.5],
        [2.0, 0.5],
        [0.5, -1.0],
        [0.5, 2.0],
        [0.0, 2.0],
    ]);

    let intersections_sc = shape.intersections(&polysegment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    let mut intersections_cs = polysegment.intersections(&shape, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    intersections_cs
        .iter_mut()
        .for_each(|i| *i = Intersection::switch(*i));

    fn slice_approx_contains(slice: &[Intersection], check: &Intersection) -> bool {
        for elem in slice {
            if approx::ulps_eq!(elem, check) {
                return true;
            }
        }
        return false;
    }

    for slice in [intersections_sc.as_slice(), intersections_cs.as_slice()].into_iter() {
        assert_eq!(slice.len(), 8);
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.0, 0.5],
                left: SegmentKey {
                    contour_idx: 0,
                    segment_idx: 3
                },
                right: SegmentKey::from_segment_idx(1)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.1, 0.5],
                left: SegmentKey {
                    contour_idx: 1,
                    segment_idx: 3
                },
                right: SegmentKey::from_segment_idx(1)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.9, 0.5],
                left: SegmentKey {
                    contour_idx: 1,
                    segment_idx: 1
                },
                right: SegmentKey::from_segment_idx(1)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [1.0, 0.5],
                left: SegmentKey::from_segment_idx(1),
                right: SegmentKey::from_segment_idx(1)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.5, 0.0],
                left: SegmentKey {
                    contour_idx: 0,
                    segment_idx: 0
                },
                right: SegmentKey::from_segment_idx(3)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.5, 0.1],
                left: SegmentKey {
                    contour_idx: 1,
                    segment_idx: 0
                },
                right: SegmentKey::from_segment_idx(3)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.5, 0.9],
                left: SegmentKey {
                    contour_idx: 1,
                    segment_idx: 2
                },
                right: SegmentKey::from_segment_idx(3)
            }
        ));
        assert!(slice_approx_contains(
            slice,
            &Intersection {
                point: [0.5, 1.0],
                left: SegmentKey {
                    contour_idx: 0,
                    segment_idx: 2
                },
                right: SegmentKey::from_segment_idx(3)
            }
        ));
    }

    {
        let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let contour = Contour::new(Polysegment::from_points(vertices));
        let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
        let hole = Contour::new(Polysegment::from_points(vertices));
        let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

        let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
        let polysegment = Polysegment::from_points(vertices);

        let intersections = polysegment.intersections(&shape, 0.0, 0);
        assert_eq!(intersections.len(), 4);

        approx::assert_abs_diff_eq!(
            intersections.get(0),
            Some(&Intersection {
                point: [1.0, 0.5],
                left: SegmentKey::from_segment_idx(1),
                right: SegmentKey::new(0, 1)
            })
        );
    }
}

#[test]
fn test_centroid() {
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let shape = Shape::new(vec![contour]).expect("valid input");

        assert_eq!(shape.centroid(), [0.5, 1.0]);
    }
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let hole = Contour::rectangle([0.6, 1.8], [0.9, 0.2]);
        let shape = Shape::new(vec![contour, hole]).expect("valid input");

        approx::assert_abs_diff_eq!(shape.centroid(), [0.421, 1.0], epsilon = 1e-3);
    }
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let shape = Shape::new(vec![contour]).expect("valid input");

        assert_eq!(shape.centroid(), [0.5, 1.0]);
    }
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let hole_1 = Contour::rectangle([0.9, 1.8], [0.6, 0.2]);
        let hole_2 = Contour::rectangle([0.1, 1.8], [0.4, 0.2]);
        let shape = Shape::new(vec![contour, hole_1, hole_2]).expect("valid input");

        approx::assert_abs_diff_eq!(shape.centroid(), [0.5, 1.0], epsilon = 1e-3);
    }
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let hole_1 = Contour::rectangle([0.6, 1.8], [0.9, 0.2]);
        let hole_2 = Contour::rectangle([0.1, 1.8], [0.4, 0.2]);
        let shape = Shape::new(vec![contour, hole_1, hole_2]).expect("valid input");

        approx::assert_abs_diff_eq!(shape.centroid(), [0.5, 1.0], epsilon = 1e-3);
    }
}
