use planar_geo::{error::ShapeConstructorError, prelude::*};
use std::f64::consts::{FRAC_PI_2, PI, SQRT_2};

#[test]
fn test_new() {
    // Shape with fillets
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([0.0, 100.0], [0.0, 0.0], [100.0, 0.0], 50.0)
                .unwrap()
                .into(),
        );
        polysegment.extend_back([100.0, 0.0]);
        polysegment.push_back(
            ArcSegment::fillet([100.0, 0.0], [0.0, 100.0], [0.0, 0.0], 10.0)
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
            ShapeConstructorError::HoleOutsideContour { .. } => (),
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
            ShapeConstructorError::HoleInsideHole { .. } => (),
            _ => unreachable!(),
        }
    }
}

#[test]
fn test_add_hole_line_segments() {
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
fn test_add_hole_arc_segments() {
    {
        let pts = &[
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
            [0.5, 0.0],
        ];
        let radii = &[0.1, 0.1, 0.1, 0.1];
        let c: Contour = Polysegment::from_fillet_chain(pts, radii).into();
        let mut shape = Shape::try_from(c).unwrap();
        assert_eq!(shape.holes().len(), 0);

        let pts = &[
            [0.5, 0.2],
            [0.2, 0.0],
            [0.2, 0.8],
            [0.8, 0.8],
            [0.8, 0.2],
            [0.5, 0.2],
        ];
        let radii = &[0.1, 0.1, 0.1, 0.1];
        let hole = Polysegment::from_fillet_chain(pts, radii).into();
        shape.add_hole(hole).unwrap();
        assert_eq!(shape.holes().len(), 1);
    }
    {
        let mut ps = Polysegment::new();
        ps.push_back(LineSegment::new([0.0, 0.0], [0.0, 0.3]).unwrap().into());
        ps.push_back(
            ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.5], 0.2, 1.5 * PI, PI)
                .unwrap()
                .into(),
        );
        ps.extend_back([0.0, 1.0]);
        ps.extend_back([1.0, 1.0]);
        ps.extend_back([1.0, 0.0]);
        let c = Contour::from(ps);
        let mut shape = Shape::try_from(c).unwrap();
        assert_eq!(shape.holes().len(), 0);

        let hole = Contour::rectangle([0.5, 0.2], [0.8, 0.8]);
        shape.add_hole(hole).unwrap();
        assert_eq!(shape.holes().len(), 1);
    }
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
    approxim::assert_abs_diff_eq!(vertices[0], [0.0, 0.0]);
    approxim::assert_abs_diff_eq!(vertices[1], [0.0, 1.0]);
    approxim::assert_abs_diff_eq!(vertices[2], [-1.0, 1.0]);
    approxim::assert_abs_diff_eq!(vertices[3], [-1.0, 0.0]);
}

#[test]
fn test_rectangle_with_hole() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let c1 = Contour::new(Polysegment::from_points(&vertices));

    let vertices = vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
    let c2 = Contour::new(Polysegment::from_points(&vertices));

    let shape = Shape::new(vec![c1, c2]).unwrap();

    assert!(
        shape
            .with_tolerance(DEFAULT_EPSILON, 0.0)
            .covers(&[0.0, 0.0],)
            .is_ok()
    );
    assert!(
        shape
            .with_tolerance(DEFAULT_EPSILON, 0.0)
            .covers(&[0.1, 0.1])
            .is_ok()
    );
    assert!(
        shape
            .with_tolerance(DEFAULT_EPSILON, 0.0)
            .covers(&[0.11, 0.11])
            .is_err()
    );
    assert!(shape.with_tolerance(0.2, 0.0).covers(&[0.11, 0.11]).is_ok());
    assert!(
        shape
            .with_tolerance(DEFAULT_EPSILON, 0.0)
            .covers(&[0.5, 0.5])
            .is_err()
    );
    assert!(
        shape
            .with_tolerance(DEFAULT_EPSILON, 0.0)
            .covers(&[0.0, -0.05])
            .is_err()
    );
    assert!(
        shape
            .with_tolerance(0.05, 0.0)
            .covers(&[0.0, -0.05])
            .is_ok()
    );
    assert!(
        shape
            .with_tolerance(0.02, 0.0)
            .covers(&[0.0, -0.05])
            .is_err()
    );
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

    let intersections_sc = shape.intersections(&polysegment);
    let mut intersections_cs = polysegment.intersections(&shape);
    intersections_cs
        .iter_mut()
        .for_each(|i| *i = Intersection::switch(*i));

    fn slice_approx_contains(slice: &[Intersection], check: &Intersection) -> bool {
        for elem in slice {
            if approxim::relative_eq!(elem, check) {
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

        let intersections = polysegment.with_tolerance(0.0, 0.0).intersections(&shape);
        assert_eq!(intersections.len(), 4);

        approxim::assert_abs_diff_eq!(
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

        approxim::assert_abs_diff_eq!(shape.centroid(), [0.421, 1.0], epsilon = 1e-3);
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

        approxim::assert_abs_diff_eq!(shape.centroid(), [0.5, 1.0], epsilon = 1e-3);
    }
    {
        let contour = Contour::rectangle([0.0, 2.0], [1.0, 0.0]);
        let hole_1 = Contour::rectangle([0.6, 1.8], [0.9, 0.2]);
        let hole_2 = Contour::rectangle([0.1, 1.8], [0.4, 0.2]);
        let shape = Shape::new(vec![contour, hole_1, hole_2]).expect("valid input");

        approxim::assert_abs_diff_eq!(shape.centroid(), [0.5, 1.0], epsilon = 1e-3);
    }
}

#[test]
fn test_covers_shape() {
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let s1 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s1).is_ok());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.0, 1.0],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s2).is_err());
        assert!(s2.covers(&s1).is_ok());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 0.5],
            [0.0, 0.5],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s2).is_err());
        assert!(s2.covers(&s1).is_err());
    }
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let h = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.9, 0.1],
            [0.9, 0.9],
            [0.1, 0.9],
        ]));
        let s1 = Shape::new(vec![c, h]).unwrap();
        assert!(s1.covers(&s1).is_ok());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.1, 0.0],
            [0.1, 1.0],
            [0.0, 1.0],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s2).is_ok());
        assert!(s2.covers(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.01, 0.1],
            [0.09, 0.1],
            [0.09, 0.9],
            [0.01, 0.9],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s2).is_ok());
        assert!(s2.covers(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.11, 0.11],
            [0.89, 0.11],
            [0.89, 0.89],
            [0.11, 0.89],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.covers(&s2).is_err());
        assert!(s2.covers(&s1).is_err());
    }
    {
        let outer = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let hole = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));
        let shape = Shape::new(vec![outer, hole]).unwrap();
        assert!(shape.covers(&shape).is_ok());
    }
}

#[test]
fn test_covers() {
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let h = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.9, 0.1],
            [0.9, 0.9],
            [0.1, 0.9],
        ]));
        let s1 = Shape::new(vec![c, h]).unwrap();
        assert!(s1.covers(&[0.05, 0.5]).is_ok());
        assert!(s1.covers(&[0.5, 0.5]).is_err());
        assert!(s1.covers(&[0.1, 0.5]).is_ok());
    }
}

#[test]
fn test_contains_point() {
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let h = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.9, 0.1],
            [0.9, 0.9],
            [0.1, 0.9],
        ]));
        let s1 = Shape::new(vec![c, h]).unwrap();
        assert!(s1.contains(&[0.05, 0.5]).is_ok());
        assert!(s1.contains(&[0.5, 0.5]).is_err());
        assert!(s1.contains(&[0.1, 0.5]).is_err());
    }
}

#[test]
fn test_contains_shape() {
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let s1 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.0, 1.0],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s2).is_err());
        assert!(s2.contains(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 0.5],
            [0.0, 0.5],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s2).is_err());
        assert!(s2.contains(&s1).is_err());
    }
    {
        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]));
        let h = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.9, 0.1],
            [0.9, 0.9],
            [0.1, 0.9],
        ]));
        let s1 = Shape::new(vec![c, h]).unwrap();
        assert!(s1.contains(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.1, 0.0],
            [0.1, 1.0],
            [0.0, 1.0],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s2).is_err());
        assert!(s2.contains(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.01, 0.1],
            [0.09, 0.1],
            [0.09, 0.9],
            [0.01, 0.9],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s2).is_ok());
        assert!(s2.contains(&s1).is_err());

        let c = Contour::new(Polysegment::from_points(&[
            [0.11, 0.11],
            [0.89, 0.11],
            [0.89, 0.89],
            [0.11, 0.89],
        ]));
        let s2 = Shape::try_from(c).unwrap();
        assert!(s1.contains(&s2).is_err());
        assert!(s2.contains(&s1).is_err());
    }
}

#[test]
fn test_overlaps_segment() {
    {
        let outer = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let hole = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));
        let shape = Shape::new(vec![outer, hole]).unwrap();

        assert!(
            shape
                .contains_any(&LineSegment::new([0.0, 0.0], [1.0, 1.0]).unwrap())
                .is_ok()
        );
        assert!(
            shape
                .contains_any(&LineSegment::new([-1.0, -1.0], [2.0, 2.0]).unwrap())
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(&LineSegment::new([0.2, 0.2], [0.8, 0.8]).unwrap())
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(&LineSegment::new([0.1, 0.1], [0.9, 0.9]).unwrap())
                .is_ok()
        );
        assert!(
            shape
                .contains_any(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        [0.1, 0.1],
                        0.9,
                        0.0,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                )
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        [0.1, 0.1],
                        SQRT_2,
                        0.0,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                )
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        [0.1, 0.1],
                        0.8,
                        0.0,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                )
                .is_ok()
        );
        assert!(
            shape
                .contains_any(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        [0.1, 0.1],
                        0.8,
                        -0.1,
                        FRAC_PI_2,
                    )
                    .unwrap(),
                )
                .is_ok()
        );
        assert!(
            shape
                .contains_any(
                    &ArcSegment::from_center_radius_start_sweep_angle(
                        [0.1, 0.1],
                        0.8,
                        -0.1,
                        FRAC_PI_2 + 0.2,
                    )
                    .unwrap(),
                )
                .is_ok()
        );
    }
}

#[test]
fn test_overlaps_contour() {
    {
        let outer = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let hole = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));
        let shape = Shape::new(vec![outer, hole]).unwrap();

        assert!(
            shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [0.0, 0.0],
                    [0.0, 1.0],
                    [1.0, 1.0],
                    [1.0, 0.0],
                ])),)
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [0.0, 0.0],
                    [0.0, -1.0],
                    [1.0, -1.0],
                    [1.0, 0.0],
                ])),)
                .is_ok()
        );
        assert!(
            shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [0.0, 0.0],
                    [0.0, 0.1],
                    [1.0, 0.1],
                    [1.0, 0.0],
                ])),)
                .is_ok()
        );
        assert!(
            shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [0.0, 0.0],
                    [0.0, 0.8],
                    [1.0, 0.8],
                    [1.0, 0.0],
                ])),)
                .is_ok()
        );
        assert!(
            !shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [0.1, 0.1],
                    [0.1, 0.9],
                    [0.9, 0.9],
                    [0.9, 0.1],
                ])),)
                .is_ok()
        );
        assert!(
            shape
                .contains_any(&Contour::new(Polysegment::from_points(&[
                    [-0.1, -0.1],
                    [-0.1, 1.1],
                    [1.1, 1.1],
                    [1.1, -0.1],
                ])),)
                .is_ok()
        );
    }
}

#[test]
fn test_overlaps_shape() {
    {
        let outer = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let hole = Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ]));
        let s1 = Shape::new(vec![outer, hole]).unwrap();
        assert!(s1.contains_any(&s1).is_ok());

        let outer = Contour::new(Polysegment::from_points(&[
            [-0.1, -0.1],
            [-0.1, 1.1],
            [1.1, 1.1],
            [1.1, -0.1],
        ]));
        let hole = Contour::new(Polysegment::from_points(&[
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]));
        let s2 = Shape::new(vec![outer, hole]).unwrap();
        assert!(s1.contains_any(&s2).is_err());
        assert!(s2.contains_any(&s1).is_err());

        let s3 = Shape::from_outer(Contour::new(Polysegment::from_points(&[
            [0.1, 0.1],
            [0.1, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ])))
        .unwrap();
        assert!(s1.contains_any(&s3).is_err());
        assert!(s3.contains_any(&s1).is_err());

        let s4 = Shape::from_outer(Contour::new(Polysegment::from_points(&[
            [0.0, 0.1],
            [0.0, 0.9],
            [0.9, 0.9],
            [0.9, 0.1],
        ])))
        .unwrap();
        assert!(s1.contains_any(&s4).is_ok());
        assert!(s4.contains_any(&s1).is_ok());
    }
}

#[test]
fn test_overlaps_various_elements() {
    let contour = Contour::new(Polysegment::from_points(&[
        [-2.0, -2.0],
        [2.0, -2.0],
        [2.0, 2.0],
        [-2.0, 2.0],
    ]));
    let hole = ArcSegment::circle([0.0, 0.0], 1.0).unwrap();
    let outer = Shape::new(vec![contour, hole.into()]).unwrap();

    // Trial geometry elements
    let mut arc =
        ArcSegment::from_center_radius_start_stop_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2).unwrap();
    let mut line_1 = LineSegment::new([0.0, 1.0], [0.0, 0.0]).unwrap();
    let mut line_2 = LineSegment::new([0.0, 0.0], [1.0, 0.0]).unwrap();
    let mut polysegment = Polysegment::new();
    polysegment.push_back(arc.clone().into());
    polysegment.push_back(line_1.clone().into());
    polysegment.push_back(line_2.clone().into());
    let mut contour = Contour::new(polysegment.clone());
    let mut shape = Shape::try_from(contour.clone()).unwrap();

    let mut count = 0.0;
    let max_counter = 200.0;
    let delta_angle = 1.0 / max_counter * std::f64::consts::TAU;

    while count < max_counter {
        if let Ok(e) = outer.contains_any(&arc) {
            panic!(
                "arc overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }
        if let Ok(e) = outer.contains_any(&line_1) {
            panic!(
                "line_1 overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }
        if let Ok(e) = outer.contains_any(&line_2) {
            panic!(
                "line_2 overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }
        if let Ok(e) = outer.contains_any(&polysegment) {
            panic!(
                "polysegment overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }
        if let Ok(e) = outer.contains_any(&contour) {
            panic!(
                "contour overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }
        if let Ok(e) = outer.contains_any(&shape) {
            panic!(
                "shape overlaps for angle {}, reason: {}",
                delta_angle * count,
                e
            );
        }

        // Rotate the elements
        arc.rotate([0.0, 0.0], delta_angle);
        line_1.rotate([0.0, 0.0], delta_angle);
        line_2.rotate([0.0, 0.0], delta_angle);
        polysegment.rotate([0.0, 0.0], delta_angle);
        contour.rotate([0.0, 0.0], delta_angle);
        shape.rotate([0.0, 0.0], delta_angle);

        count += 1.0;
    }
}
