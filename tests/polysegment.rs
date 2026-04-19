use approx::{assert_abs_diff_eq, assert_ulps_eq};
use planar_geo::{polysegment::area_signed, prelude::*};
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[test]
fn test_from_bounding_box() {
    let bounding_box = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let polysegment: Polysegment = bounding_box.into();
    assert_eq!(polysegment.len(), 4);
}

#[test]
fn test_from_iter() {
    let mut polysegment = Polysegment::from_iter(
        [
            ArcSegment::fillet([-1.5, 2.0], [-1.0, 0.0], [1.0, 0.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([-1.0, 0.0], [1.0, 0.0], [1.5, 2.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([1.0, 0.0], [1.5, 2.0], [-1.5, 2.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([1.5, 2.0], [-1.5, 2.0], [-1.0, 0.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
        ]
        .into_iter(),
    );
    assert_eq!(polysegment.len(), 7);

    polysegment.close();
    assert_eq!(polysegment.len(), 8);
}

#[test]
fn test_to_bounding_box() {
    {
        let bounding_box = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
        let polysegment: Polysegment = bounding_box.into();
        let bb = BoundingBox::from(&polysegment);
        assert_eq!(bounding_box, bb);
    }
    {
        // Empty polysegment
        let polysegment = Polysegment::new();
        let bb = BoundingBox::from(&polysegment);
        assert_eq!(bb.xmin(), 0.0);
        assert_eq!(bb.xmax(), 0.0);
        assert_eq!(bb.ymin(), 0.0);
        assert_eq!(bb.ymax(), 0.0);
    }
}

#[test]
fn test_area_calculation() {
    {
        // Triangle
        let vertices: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        approx::assert_abs_diff_eq!(area_signed(&mut vertices.into_iter()), 0.5);
    }

    {
        // Rectangle
        let vertices: Vec<[f64; 2]> = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]];
        approx::assert_abs_diff_eq!(area_signed(&mut vertices.into_iter()), 2.0);
    }

    {
        // Closed rectangle (last vertex equals first vertex)
        let vertices: Vec<[f64; 2]> =
            vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
        approx::assert_abs_diff_eq!(area_signed(&mut vertices.into_iter()), 2.0);
    }

    {
        let pt1 = [0.0, 0.0];
        let pt2 = [1.0, 0.0];
        let pt3 = [1.0, 1.0];
        assert_eq!(area_signed([pt1, pt2, pt3].into_iter()), 0.5);
        assert_eq!(area_signed([pt2, pt1, pt3].into_iter()), -0.5);
    }
    {
        let pt1 = [0.0, 0.0];
        let pt2 = [1.0, 0.0];
        let pt3 = [1.0, 1.0];
        let pt4 = [0.0, 1.0];
        assert_eq!(area_signed([pt1, pt2, pt3, pt4].into_iter()), 1.0);
        assert_eq!(area_signed([pt1, pt2, pt3, pt4].into_iter().rev()), -1.0);
    }
}

#[test]
fn test_from_points() {
    // Triangle
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let path = Polysegment::from_points(&vertices);
    let vertices: Vec<[f64; 2]> = path.points().collect();

    assert_eq!(vertices.len(), 3);
    assert_eq!(vertices[0], [0.0, 0.0]);
    assert_eq!(vertices[2], [1.0, 1.0]);

    // Three consecutive vertices are identical -> unite them
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let path = Polysegment::from_points(&vertices);
    assert_eq!(path.len(), 2);

    // Failed to create a polysegment because the vertices are identical
    let vertices = vec![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
    let empty_polysegment = Polysegment::from_points(&vertices);
    assert_eq!(empty_polysegment.num_segments(), 0);
}

#[test]
fn test_push_and_pop() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0]];
    let mut polysegment = Polysegment::from_points(&vertices);

    // Currently, we have one segment (a single line segment)
    assert_eq!(polysegment.len(), 1);

    // Add another line segment which is directly connected to the polysegment
    let segment: Segment =
        LineSegment::new([1.0, 0.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();

    polysegment.push_back(segment.clone());
    assert_eq!(polysegment.len(), 2);
    let removed_segment = polysegment.pop_back().unwrap();
    assert_eq!(segment, removed_segment);

    // Add another line segment which is not directly connected to the polysegment
    let segment: Segment =
        LineSegment::new([3.0, 0.0], [4.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();
    polysegment.push_back(segment.clone());
    assert_eq!(polysegment.len(), 3);
    let removed_segment = polysegment.pop_back().unwrap();
    assert_eq!(segment, removed_segment);

    // Check the connection line segment which was put between the original
    // polysegment and the previously added (and removed) line segment
    let removed_segment = polysegment.pop_back().unwrap();
    let segment: Segment =
        LineSegment::new([1.0, 0.0], [3.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();
    assert_eq!(segment, removed_segment);
}

#[test]
fn test_push_extend() {
    let mut polysegment = Polysegment::new();
    polysegment.push_back(
        ArcSegment::fillet([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], 0.5, 0.0, 0)
            .unwrap()
            .into(),
    );
    polysegment.extend_front([0.0, 0.0]);
    polysegment.extend_back([1.0, 1.0]);
    polysegment.extend_back([0.0, 1.0]);
    polysegment.extend_back([0.0, 0.0]);

    // Check the segments
    assert_eq!(polysegment.len(), 5);
    assert_ulps_eq!(polysegment[0].start(), [0.0, 0.0]);
    assert_ulps_eq!(polysegment[0].stop(), [0.5, 0.0]);
    assert_ulps_eq!(polysegment[1].start(), [0.5, 0.0]);
    assert_ulps_eq!(polysegment[1].stop(), [1.0, 0.5]);
    assert_ulps_eq!(polysegment[2].start(), [1.0, 0.5]);
    assert_ulps_eq!(polysegment[2].stop(), [1.0, 1.0]);
    assert_ulps_eq!(polysegment[3].start(), [1.0, 1.0]);
    assert_ulps_eq!(polysegment[3].stop(), [0.0, 1.0]);
    assert_ulps_eq!(polysegment[4].start(), [0.0, 1.0]);
    assert_ulps_eq!(polysegment[4].stop(), [0.0, 0.0]);

    // Now check the vertices
    let verts: Vec<[f64; 2]> = polysegment.points().collect();
    assert_eq!(verts.len(), 6);
    assert_ulps_eq!(verts[0], [0.0, 0.0]);
    assert_ulps_eq!(verts[1], [0.5, 0.0]);
    assert_ulps_eq!(verts[2], [1.0, 0.5]);
    assert_ulps_eq!(verts[3], [1.0, 1.0]);
    assert_ulps_eq!(verts[4], [0.0, 1.0]);
    assert_ulps_eq!(verts[5], [0.0, 0.0]);
}

#[test]
fn test_single_cut() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0]];
    let line = Polysegment::from_points(&vertices);

    let cut_vertices = vec![[0.5, -0.5], [0.5, 0.5]];
    let cut = Polysegment::from_points(&cut_vertices);

    let separated_lines = line.intersection_cut(&cut, 0.0, 0);

    // This cut results in three separate entitites
    assert_eq!(separated_lines.len(), 2);
    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    assert_eq!(verts[0], [0.0, 0.0]);
    assert_eq!(verts[1], [0.5, 0.0]);

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    assert_eq!(verts[0], [0.5, 0.0]);
    assert_eq!(verts[1], [1.0, 0.0]);
}

#[test]
fn test_cutted_slot_contour() {
    let line = Polysegment::from_points(&[
        [0.001, 0.0],
        [0.001, 0.001],
        [0.004158825456719743, 0.001],
        [0.007653668647301796, 0.007053245970574237],
        [-0.007653668647301796, 0.007053245970574237],
        [-0.004158825456719743, 0.001],
        [-0.001, 0.001],
        [-0.001, 0.0],
        [0.001, 0.0],
    ]);
    let cutter = Polysegment::from_points(&[
        [-0.015307337294603592, 0.001],
        [0.015307337294603592, 0.001],
    ]);

    let separated_lines = line.intersection_cut(&cutter, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    assert_eq!(separated_lines.len(), 5);

    let vertices = &[
        vec![[0.001, 0.0], [0.001, 0.001]],
        vec![[0.001, 0.001], [0.004158825456719743, 0.001]],
        vec![
            [0.004158825456719743, 0.001],
            [0.007653668647301796, 0.007053245970574237],
            [-0.007653668647301796, 0.007053245970574237],
            [-0.004158825456719743, 0.001],
        ],
        vec![[-0.004158825456719743, 0.001], [-0.001, 0.001]],
        vec![[-0.001, 0.001], [-0.001, 0.0]],
    ];
    for (vertices, line) in vertices.iter().zip(separated_lines.iter()) {
        for (v1, v2) in vertices.iter().zip(line.points()) {
            assert_eq!(*v1, v2);
        }
    }
}

#[test]
fn test_cut_vertical_horizontal() {
    let vertical_line = Polysegment::from_points(&[[0.0, 0.0], [0.0, 1.0]]);
    let horizontal_line = Polysegment::from_points(&[[-0.5, 0.5], [0.5, 0.5]]);

    // Cut the vertical line with the horizontal line
    let separated_lines =
        vertical_line.intersection_cut(&horizontal_line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    assert_eq!(verts[0], [0.0, 0.0]);
    assert_eq!(verts[1], [0.0, 0.5]);

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    assert_eq!(verts[0], [0.0, 0.5]);
    assert_eq!(verts[1], [0.0, 1.0]);

    // Cut the horizontal line with the vertical line
    let separated_lines =
        horizontal_line.intersection_cut(&vertical_line, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    assert_eq!(verts[0], [-0.5, 0.5]);
    assert_eq!(verts[1], [0.0, 0.5]);

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    assert_eq!(verts[0], [0.0, 0.5]);
    assert_eq!(verts[1], [0.5, 0.5]);
}

#[test]
fn test_second_segment_intersection_before_first() {
    let vertices = [[0.0, 0.0], [2.0, 0.0]];
    let cutted = Polysegment::from_points(vertices.as_slice());

    let cut_vertices = vec![[1.5, 1.0], [1.5, -1.0], [0.5, -1.0], [0.5, 1.0]];
    let cut = Polysegment::from_points(&cut_vertices);

    let separated_lines = cutted.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    assert_eq!(separated_lines.len(), 3);

    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    assert_ulps_eq!(verts[0], [0.0, 0.0]);
    assert_ulps_eq!(verts[1], [0.5, 0.0]);

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    assert_ulps_eq!(verts[0], [0.5, 0.0]);
    assert_ulps_eq!(verts[1], [1.5, 0.0]);

    let verts: Vec<[f64; 2]> = separated_lines[2].points().collect();
    assert_ulps_eq!(verts[0], [1.5, 0.0]);
    assert_ulps_eq!(verts[1], [2.0, 0.0]);

    // Circle
    // =====================================================================================
    let cutted: Polysegment =
        ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2, 0.0, 0)
            .unwrap()
            .into();

    let cut_vertices = vec![[0.5, 1.0], [0.0, 0.0], [1.0, 0.5]];
    let mut cut = Polysegment::from_points(&cut_vertices);

    let separated_lines = cutted.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    assert_eq!(separated_lines.len(), 3);

    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    approx::assert_abs_diff_eq!(
        verts.last().unwrap(),
        &[0.89442719, 0.4472135],
        epsilon = 1e-5
    );

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    approx::assert_abs_diff_eq!(
        verts.last().unwrap(),
        &[0.4472135, 0.89442719],
        epsilon = 1e-5
    );

    // Invert the cut
    cut.reverse();

    let separated_lines = cutted.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);
    assert_eq!(separated_lines.len(), 3);

    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    approx::assert_abs_diff_eq!(
        verts.last().unwrap(),
        &[0.89442719, 0.4472135],
        epsilon = 1e-5
    );

    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    approx::assert_abs_diff_eq!(
        verts.last().unwrap(),
        &[0.4472135, 0.89442719],
        epsilon = 1e-5
    );
}

#[test]
fn test_intersection_cut() {
    let mut polysegment = Polysegment::new();
    polysegment.push_back(
        ArcSegment::fillet([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], 0.5, 0.0, 0)
            .unwrap()
            .into(),
    );
    polysegment.extend_front([0.0, 0.0]);
    polysegment.extend_back([1.0, 1.0]);
    polysegment.extend_back([0.0, 1.0]);

    let cut_vertices = vec![[0.5, 1.0], [1.0, 0.0], [1.0, -0.5], [0.5, -0.5], [0.5, 0.5]];
    let cut = Polysegment::from_points(&cut_vertices);

    let separated_lines = polysegment.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

    // This cut results in four separate entitites
    assert_eq!(separated_lines.len(), 4);

    // Check the vertices of the first entry
    let verts: Vec<[f64; 2]> = separated_lines[0].points().collect();
    assert_eq!(verts[0], [0.0, 0.0]);
    assert_ulps_eq!(&verts[1], &[0.5, 0.0]);

    // Check the vertices of the second entry
    let verts: Vec<[f64; 2]> = separated_lines[1].points().collect();
    assert_abs_diff_eq!(verts[0], [0.5, 0.0]);
    assert_abs_diff_eq!(&verts[1], &[0.9, 0.2], epsilon = 1e-15);

    // Check the vertices of the third entry
    let verts: Vec<[f64; 2]> = separated_lines[2].points().collect();
    assert_abs_diff_eq!(&verts[0], &[0.9, 0.2], epsilon = 1e-15);
    assert_abs_diff_eq!(&verts[1], &[1.0, 0.5], epsilon = 1e-15);
    assert_eq!(verts[2], [1.0, 1.0]);
    assert_eq!(verts[3], [0.5, 1.0]);

    // Check the vertices of the first entry
    let verts: Vec<[f64; 2]> = separated_lines[3].points().collect();
    assert_eq!(verts[0], [0.5, 1.0]);
    assert_eq!(verts[1], [0.0, 1.0]);
}

#[test]
fn test_self_intersection() {
    // Open polysegment
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 0.5, 0.0, 0)
                .unwrap()
                .into(),
        );
        polysegment.extend_front([0.0, 0.0]);

        // Intersect the line with itself
        let intersections =
            polysegment.intersections_polysegment(&polysegment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

        // No self-intersection
        assert_eq!(intersections.count(), 0);
    }

    // Closed polysegment
    {
        let mut polysegment = Polysegment::new();
        polysegment.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 0.5, 0.0, 0)
                .unwrap()
                .into(),
        );
        polysegment.extend_front([0.0, 0.0]);
        polysegment.extend_back([0.0, 0.0]);

        // Intersect the line with itself
        let intersections =
            polysegment.intersections_polysegment(&polysegment, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

        // Polysegment touches itself at the end
        assert_eq!(intersections.count(), 1);

        // Intersect the line with itself
        let intersections = polysegment.intersections_polysegment_par(
            &polysegment,
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        );

        // No self-intersection
        assert_eq!(intersections.count(), 1);
    }
}

#[test]
fn test_intersections_primitive() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;
    {
        let mut ps = Polysegment::new();
        ps.push_back(
            ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], PI, e, m)
                .unwrap()
                .into(),
        );
        ps.push_back(
            ArcSegment::from_start_center_angle([-1.0, 0.0], [-2.0, 0.0], PI, e, m)
                .unwrap()
                .into(),
        );
        ps.extend_back([-2.0, -2.0]);
        ps.extend_back([-1.0, 0.0]);
        ps.extend_back([0.0, -1.0]);
        ps.extend_back([0.0, 0.0]);
        ps.extend_back([0.5, 0.0]);

        let ls = LineSegment::new([-1.5, 0.0], [-10.0, 0.0], e, m).unwrap();
        let mut intersections = ps.intersections_primitive(&ls, e, m);
        assert_eq!(
            intersections.next(),
            Some(Intersection {
                point: [-3.0, 0.0],
                left: SegmentKey::from_segment_idx(1),
                right: SegmentKey::from_segment_idx(0)
            })
        );
        assert_eq!(
            intersections.next(),
            Some(Intersection {
                point: [-3.0, 0.0],
                left: SegmentKey::from_segment_idx(2),
                right: SegmentKey::from_segment_idx(0)
            })
        );
        assert_eq!(intersections.next(), None);
    }
}

#[test]
fn test_append() {
    // Appending and closing
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = Polysegment::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = Polysegment::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has 6 vertices
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);

    // Appending and not closing
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = Polysegment::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = Polysegment::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has five elements
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);

    // Appending and closing non-closed paths
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = Polysegment::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = Polysegment::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has four elements
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);
}

#[test]
fn test_intersection_cut_fillets() {
    fn cut_by_height(polysegment: &Polysegment, height: f64) {
        let bounding_box = BoundingBox::new(-2.0, 2.0, height, 4.0);
        let new_polylines = polysegment.intersection_cut(
            &Polysegment::from(&bounding_box),
            DEFAULT_EPSILON,
            DEFAULT_MAX_ULPS,
        );
        assert_eq!(3, new_polylines.len());

        // Connect the third and first line
        let mut line = new_polylines[2].clone();
        line.append(&mut new_polylines[0].clone());

        for line in [line, new_polylines[1].clone()].iter() {
            assert_ulps_eq!(height, line.front().unwrap().start()[1]);
            assert_ulps_eq!(height, line.back().unwrap().stop()[1]);
        }
    }

    let polysegment = Polysegment::from_iter(
        [
            ArcSegment::fillet([-1.5, 2.0], [-1.0, 0.0], [1.0, 0.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([-1.0, 0.0], [1.0, 0.0], [1.5, 2.0], 0.1, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([1.0, 0.0], [1.5, 2.0], [-1.5, 2.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
            ArcSegment::fillet([1.5, 2.0], [-1.5, 2.0], [-1.0, 0.0], 0.2, 0.0, 0)
                .unwrap()
                .into(),
        ]
        .into_iter(),
    );

    cut_by_height(&polysegment, 1e-3);
    cut_by_height(&polysegment, 0.5e-3);
    cut_by_height(&polysegment, 0.2e-3);
    cut_by_height(&polysegment, 0.1e-3);
    cut_by_height(&polysegment, 0.01e-3);
}

#[test]
fn test_rotational_pattern() {
    let org_polysegment = Polysegment::from_points(&[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0]]);
    {
        // Do nothing
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.rotational_pattern([0.0, 1.0], FRAC_PI_2, 0);
        assert_eq!(org_polysegment, mod_polysegment);
    }
    {
        // Single repetition
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.rotational_pattern([0.0, 1.0], FRAC_PI_2, 1);

        assert_eq!(mod_polysegment.num_segments(), 5);
        let pts: Vec<_> = mod_polysegment.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [0.0, 2.0]);
        assert_ulps_eq!(pts[4], [0.0, 3.0]);
        assert_ulps_eq!(pts[5], [-1.0, 3.0]);
    }
    {
        // Two repetitions
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.rotational_pattern([0.0, 1.0], FRAC_PI_2, 2);

        assert_eq!(mod_polysegment.num_segments(), 8);
        let pts: Vec<_> = mod_polysegment.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [0.0, 2.0]);
        assert_ulps_eq!(pts[4], [0.0, 3.0]);
        assert_ulps_eq!(pts[5], [-1.0, 3.0]);
        assert_ulps_eq!(pts[6], [-1.0, 1.0]);
        assert_ulps_eq!(pts[7], [-2.0, 1.0]);
        assert_ulps_eq!(pts[8], [-2.0, 0.0]);
    }
}

#[test]
fn test_translational_pattern() {
    let org_polysegment = Polysegment::from_points(&[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0]]);
    {
        // Do nothing
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.translational_pattern([0.0, 1.0], 0);
        assert_eq!(org_polysegment, mod_polysegment);
    }
    {
        // Single repetition
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.translational_pattern([0.0, 1.0], 1);

        assert_eq!(mod_polysegment.num_segments(), 5);
        let pts: Vec<_> = mod_polysegment.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [1.0, 2.0]);
        assert_ulps_eq!(pts[4], [2.0, 2.0]);
        assert_ulps_eq!(pts[5], [2.0, 3.0]);
    }
    {
        // Two repetitions
        let mut mod_polysegment = org_polysegment.clone();
        mod_polysegment.translational_pattern([0.0, 1.0], 2);

        assert_eq!(mod_polysegment.num_segments(), 8);
        let pts: Vec<_> = mod_polysegment.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [1.0, 2.0]);
        assert_ulps_eq!(pts[4], [2.0, 2.0]);
        assert_ulps_eq!(pts[5], [2.0, 3.0]);
        assert_ulps_eq!(pts[6], [1.0, 3.0]);
        assert_ulps_eq!(pts[7], [2.0, 3.0]);
        assert_ulps_eq!(pts[8], [2.0, 4.0]);
    }
}

#[test]
fn test_polygonize() {
    {
        let arc: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, TAU, 0.0, 0)
                .unwrap()
                .into();
        let path = Polysegment::from(arc);

        let mut iter = path.polygonize(Polygonizer::PerType {
            arc: SegmentPolygonizer::MaximumAngle(FRAC_PI_2),
            straight: SegmentPolygonizer::InnerSegments(1),
        });

        approx::assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }

    {
        let arc: Segment =
            ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 0.0, -TAU, 0.0, 0)
                .unwrap()
                .into();
        let path = Polysegment::from(arc);

        let mut iter = path.polygonize(Polygonizer::PerType {
            arc: SegmentPolygonizer::MaximumAngle(FRAC_PI_2),
            straight: SegmentPolygonizer::InnerSegments(1),
        });

        approx::assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [0.0, -2.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [-2.0, 0.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [0.0, 2.0], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(iter.next().unwrap(), [2.0, 0.0], epsilon = 1e-10);
        assert_eq!(iter.next(), None);
    }
}

#[test]
fn test_straight_line_polysegment() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let polysegment = Polysegment::from_points(&vertices);

    // Points which are covered in the first segment
    assert!(polysegment.covers_point([0.5, 0.0], 0.0, 0));
    assert!(!polysegment.covers_point([0.5, 0.1], 0.0, 0));
    assert!(polysegment.covers_point([0.5, 0.1], 0.11, 0));
    assert!(!polysegment.covers_point([0.5, 0.2], 0.11, 0));

    // Points which are covered in the second segment
    assert!(polysegment.covers_point([1.0, 0.5], 0.0, 0));
    assert!(!polysegment.covers_point([1.1, 0.5], 0.0, 0));
    assert!(polysegment.covers_point([1.1, 0.5], 0.11, 0));
    assert!(!polysegment.covers_point([1.2, 0.5], 0.11, 0));

    // Points which are covered in both segments
    assert!(polysegment.covers_point([1.0, 0.0], 0.0, 0));
    assert!(!polysegment.covers_point([0.98, 0.02], 0.0, 0));
    assert!(polysegment.covers_point([0.98, 0.02], 0.11, 0));
    assert!(!polysegment.covers_point([1.1, -0.1], 0.11, 0));
}

#[test]
fn test_intersection_polysegments() {
    let first = Polysegment::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let second = Polysegment::from_points(&[[0.0, 0.5], [2.0, 0.5]]);

    let vec: Vec<_> = first.intersections_polysegment(&second, 0.0, 0).collect();
    assert_eq!(vec.len(), 1);

    let i = vec[0];
    assert_eq!(
        i,
        Intersection {
            point: [1.0, 0.5],
            left: SegmentKey::from_segment_idx(1),
            right: SegmentKey::from_segment_idx(0)
        }
    );
}

#[test]
fn test_intersection_line_polysegment() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let polysegment = Polysegment::from_points(&vertices);

    let line = Line::from_point_angle([0.5, 0.5], 0.0);

    {
        let mut intersections = polysegment.intersections_primitive(&line, 0.0, 0);

        approx::assert_abs_diff_eq!(
            intersections.next(),
            Some(Intersection {
                point: [1.0, 0.5],
                left: SegmentKey::from_segment_idx(1),
                right: Default::default()
            })
        );
        approx::assert_abs_diff_eq!(
            intersections.next(),
            Some(Intersection {
                point: [0.0, 0.5],
                left: SegmentKey::from_segment_idx(3),
                right: Default::default()
            })
        );
        assert!(intersections.next().is_none());
    }
    {
        let mut intersections = polysegment.intersections(&line, 0.0, 0).into_iter();

        approx::assert_abs_diff_eq!(
            intersections.next(),
            Some(Intersection {
                point: [1.0, 0.5],
                left: SegmentKey::from_segment_idx(1),
                right: Default::default()
            })
        );
        approx::assert_abs_diff_eq!(
            intersections.next(),
            Some(Intersection {
                point: [0.0, 0.5],
                left: SegmentKey::from_segment_idx(3),
                right: Default::default()
            })
        );
        assert!(intersections.next().is_none());
    }
}

#[test]
fn covers_segment() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let mut polyseg = Polysegment::new();
    polyseg.push_back(
        ArcSegment::from_start_center_angle([-1.0, 0.0], [0.0, 0.0], 0.5 * -PI, e, m)
            .unwrap()
            .into(),
    );
    polyseg.push_back(
        ArcSegment::from_start_center_angle([0.0, 1.0], [0.0, 0.0], 0.5 * -PI, e, m)
            .unwrap()
            .into(),
    );
    polyseg.extend_back([2.0, 0.0]);
    polyseg.extend_back([4.0, 0.0]);
    polyseg.extend_back([4.0, 2.0]);
    polyseg.extend_back([4.0, 4.0]);

    {
        // Contains line segment
        assert!(polyseg.covers_segment(
            &LineSegment::new([1.0, 0.0], [2.0, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([2.0, 0.0], [1.0, 0.0], e, m).unwrap(),
            e,
            m
        ));

        assert!(polyseg.covers_segment(
            &LineSegment::new([1.5, 0.0], [2.5, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([2.5, 0.0], [1.5, 0.0], e, m).unwrap(),
            e,
            m
        ));

        assert!(polyseg.covers_segment(
            &LineSegment::new([1.0, 0.0], [4.0, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([4.0, 0.0], [1.0, 0.0], e, m).unwrap(),
            e,
            m
        ));

        assert!(polyseg.covers_segment(
            &LineSegment::new([2.1, 0.0], [3.9, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([3.9, 0.0], [2.1, 0.0], e, m).unwrap(),
            e,
            m
        ));

        assert!(!polyseg.covers_segment(
            &LineSegment::new([1.0, 0.2], [2.0, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(!polyseg.covers_segment(
            &LineSegment::new([1.0, 0.2], [2.0, 0.2], e, m).unwrap(),
            e,
            m
        ));

        assert!(polyseg.covers_segment(
            &LineSegment::new([4.0, 4.0], [4.0, 0.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([4.0, 0.0], [4.0, 4.0], e, m).unwrap(),
            e,
            m
        ));

        assert!(polyseg.covers_segment(
            &LineSegment::new([4.0, 4.0], [4.0, 3.0], e, m).unwrap(),
            e,
            m
        ));
        assert!(polyseg.covers_segment(
            &LineSegment::new([4.0, 3.0], [4.0, 4.0], e, m).unwrap(),
            e,
            m
        ));
    }
    {
        // Contains arc segment
        let quarter_arc =
            ArcSegment::from_start_center_angle([-1.0, 0.0], [0.0, 0.0], 0.5 * -PI, e, m).unwrap();
        assert!(polyseg.covers_segment(&quarter_arc, e, m));

        let quarter_arc_shifted = ArcSegment::from_center_radius_start_offset_angle(
            [0.0, 0.0],
            1.0,
            0.75 * PI,
            -0.5 * PI,
            e,
            m,
        )
        .unwrap();
        assert!(polyseg.covers_segment(&quarter_arc_shifted, e, m));
    }
}

#[test]
fn covers_segment_circle() {
    let e = DEFAULT_EPSILON;
    let m = DEFAULT_MAX_ULPS;

    let circle = ArcSegment::circle([0.0, 0.0], 1.0).unwrap();
    let polyseg = Polysegment::from(circle);

    let arc =
        ArcSegment::from_start_center_angle([-1.0, 0.0], [0.0, 0.0], 0.5 * -PI, e, m).unwrap();
    assert!(polyseg.covers_segment(&arc, e, m));

    let arc = ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], 1.5 * -PI, e, m).unwrap();
    assert!(polyseg.covers_segment(&arc, e, m));

    let arc = ArcSegment::from_start_center_angle([1.0, 0.0], [0.0, 0.0], 1.5 * PI, e, m).unwrap();
    assert!(polyseg.covers_segment(&arc, e, m));

    let arc =
        ArcSegment::from_center_radius_start_stop_angle([0.0, 0.0], 1.0, -1.0, 1.0, e, m).unwrap();
    assert!(polyseg.covers_segment(&arc, e, m));

    let arc =
        ArcSegment::from_center_radius_start_stop_angle([0.0, 0.0], 1.0, 1.0, -1.0, e, m).unwrap();
    assert!(polyseg.covers_segment(&arc, e, m));

    let arc =
        ArcSegment::from_center_radius_start_stop_angle([0.0, 0.0], 2.0, 1.0, -1.0, e, m).unwrap();
    assert!(!polyseg.covers_segment(&arc, e, m));

    let arc =
        ArcSegment::from_center_radius_start_stop_angle([0.5, 0.0], 1.0, 1.0, -1.0, e, m).unwrap();
    assert!(!polyseg.covers_segment(&arc, e, m));
}
