use approx::{assert_abs_diff_eq, assert_ulps_eq};
use planar_geo::{prelude::*, segment_chain::area_signed};
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, TAU};

#[test]
fn test_from_bounding_box() {
    let bounding_box = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
    let chain: SegmentChain = bounding_box.into();
    assert_eq!(chain.len(), 4);
}

#[test]
fn test_from_iter() {
    let mut chain = SegmentChain::from_iter(
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
    assert_eq!(chain.len(), 7);

    chain.close();
    assert_eq!(chain.len(), 8);
}

#[test]
fn test_to_bounding_box() {
    {
        let bounding_box = BoundingBox::new(0.0, 1.0, 0.0, 1.0);
        let chain: SegmentChain = bounding_box.into();
        let bb = BoundingBox::from(&chain);
        assert_eq!(bounding_box, bb);
    }
    {
        // Empty chain
        let chain = SegmentChain::new();
        let bb = BoundingBox::from(&chain);
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
    let path = SegmentChain::from_points(&vertices);
    let vertices: Vec<[f64; 2]> = path.points().collect();

    assert_eq!(vertices.len(), 3);
    assert_eq!(vertices[0], [0.0, 0.0]);
    assert_eq!(vertices[2], [1.0, 1.0]);

    // Three consecutive vertices are identical -> unite them
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let path = SegmentChain::from_points(&vertices);
    assert_eq!(path.len(), 2);

    // Failed to create a segment_chain because the vertices are identical
    let vertices = vec![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
    let empty_segment_chain = SegmentChain::from_points(&vertices);
    assert_eq!(empty_segment_chain.num_segments(), 0);
}

#[test]
fn test_push_and_pop() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0]];
    let mut segment_chain = SegmentChain::from_points(&vertices);

    // Currently, we have one segment (a single line segment)
    assert_eq!(segment_chain.len(), 1);

    // Add another line segment which is directly connected to the segment_chain
    let segment: Segment =
        LineSegment::new([1.0, 0.0], [2.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();

    segment_chain.push_back(segment.clone());
    assert_eq!(segment_chain.len(), 2);
    let removed_segment = segment_chain.pop_back().unwrap();
    assert_eq!(segment, removed_segment);

    // Add another line segment which is not directly connected to the segment_chain
    let segment: Segment =
        LineSegment::new([3.0, 0.0], [4.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();
    segment_chain.push_back(segment.clone());
    assert_eq!(segment_chain.len(), 3);
    let removed_segment = segment_chain.pop_back().unwrap();
    assert_eq!(segment, removed_segment);

    // Check the connection line segment which was put between the original segment_chain and the previously added (and removed) line segment
    let removed_segment = segment_chain.pop_back().unwrap();
    let segment: Segment =
        LineSegment::new([1.0, 0.0], [3.0, 0.0], DEFAULT_EPSILON, DEFAULT_MAX_ULPS)
            .unwrap()
            .into();
    assert_eq!(segment, removed_segment);
}

#[test]
fn test_push_extend() {
    let mut chain = SegmentChain::new();
    chain.push_back(
        ArcSegment::fillet([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], 0.5, 0.0, 0)
            .unwrap()
            .into(),
    );
    chain.extend_front([0.0, 0.0]);
    chain.extend_back([1.0, 1.0]);
    chain.extend_back([0.0, 1.0]);
    chain.extend_back([0.0, 0.0]);

    // Check the segments
    assert_eq!(chain.len(), 5);
    assert_ulps_eq!(chain[0].start(), [0.0, 0.0]);
    assert_ulps_eq!(chain[0].stop(), [0.5, 0.0]);
    assert_ulps_eq!(chain[1].start(), [0.5, 0.0]);
    assert_ulps_eq!(chain[1].stop(), [1.0, 0.5]);
    assert_ulps_eq!(chain[2].start(), [1.0, 0.5]);
    assert_ulps_eq!(chain[2].stop(), [1.0, 1.0]);
    assert_ulps_eq!(chain[3].start(), [1.0, 1.0]);
    assert_ulps_eq!(chain[3].stop(), [0.0, 1.0]);
    assert_ulps_eq!(chain[4].start(), [0.0, 1.0]);
    assert_ulps_eq!(chain[4].stop(), [0.0, 0.0]);

    // Now check the vertices
    let verts: Vec<[f64; 2]> = chain.points().collect();
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
    let line = SegmentChain::from_points(&vertices);

    let cut_vertices = vec![[0.5, -0.5], [0.5, 0.5]];
    let cut = SegmentChain::from_points(&cut_vertices);

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
    let line = SegmentChain::from_points(&[
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
    let cutter = SegmentChain::from_points(&[
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
    let vertical_line = SegmentChain::from_points(&[[0.0, 0.0], [0.0, 1.0]]);
    let horizontal_line = SegmentChain::from_points(&[[-0.5, 0.5], [0.5, 0.5]]);

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
    let cutted = SegmentChain::from_points(vertices.as_slice());

    let cut_vertices = vec![[1.5, 1.0], [1.5, -1.0], [0.5, -1.0], [0.5, 1.0]];
    let cut = SegmentChain::from_points(&cut_vertices);

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
    let cutted: SegmentChain =
        ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 1.0, 0.0, FRAC_PI_2, 0.0, 0)
            .unwrap()
            .into();

    let cut_vertices = vec![[0.5, 1.0], [0.0, 0.0], [1.0, 0.5]];
    let mut cut = SegmentChain::from_points(&cut_vertices);

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
    let mut chain = SegmentChain::new();
    chain.push_back(
        ArcSegment::fillet([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], 0.5, 0.0, 0)
            .unwrap()
            .into(),
    );
    chain.extend_front([0.0, 0.0]);
    chain.extend_back([1.0, 1.0]);
    chain.extend_back([0.0, 1.0]);

    let cut_vertices = vec![[0.5, 1.0], [1.0, 0.0], [1.0, -0.5], [0.5, -0.5], [0.5, 0.5]];
    let cut = SegmentChain::from_points(&cut_vertices);

    let separated_lines = chain.intersection_cut(&cut, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

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
    // Open segment_chain
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 0.5, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.extend_front([0.0, 0.0]);

        // Intersect the line with itself
        let intersections =
            chain.intersections_segment_chain(&chain, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

        // No self-intersection
        assert_eq!(intersections.count(), 0);
    }

    // Closed segment_chain
    {
        let mut chain = SegmentChain::new();
        chain.push_back(
            ArcSegment::fillet([1.0, 0.0], [1.0, 1.0], [0.0, 1.0], 0.5, 0.0, 0)
                .unwrap()
                .into(),
        );
        chain.extend_front([0.0, 0.0]);
        chain.extend_back([0.0, 0.0]);

        // Intersect the line with itself
        let intersections =
            chain.intersections_segment_chain(&chain, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

        // Segment chain touches itself at the end
        assert_eq!(intersections.count(), 1);

        // Intersect the line with itself
        let intersections =
            chain.intersections_segment_chain_par(&chain, DEFAULT_EPSILON, DEFAULT_MAX_ULPS);

        // No self-intersection
        assert_eq!(intersections.count(), 1);
    }
}

#[test]
fn test_append() {
    // Appending and closing
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = SegmentChain::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = SegmentChain::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has 6 vertices
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);

    // Appending and not closing
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = SegmentChain::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = SegmentChain::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has five elements
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);

    // Appending and closing non-closed paths
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let mut path1 = SegmentChain::from_points(&vertices);
    let vertices = vec![[0.0, 1.0], [-0.5, 0.5]];
    let mut path2 = SegmentChain::from_points(&vertices);
    path1.append(&mut path2);
    let verts: Vec<[f64; 2]> = path1.points().collect();

    // The new path now has four elements
    assert_eq!(verts.len(), 5);
    assert_eq!(verts[4], [-0.5, 0.5]);
}

#[test]
fn test_intersection_cut_fillets() {
    fn cut_by_height(segment_chain: &SegmentChain, height: f64) {
        let bounding_box = BoundingBox::new(-2.0, 2.0, height, 4.0);
        let new_polylines = segment_chain.intersection_cut(
            &SegmentChain::from(&bounding_box),
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

    let chain = SegmentChain::from_iter(
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

    cut_by_height(&chain, 1e-3);
    cut_by_height(&chain, 0.5e-3);
    cut_by_height(&chain, 0.2e-3);
    cut_by_height(&chain, 0.1e-3);
    cut_by_height(&chain, 0.01e-3);
}

#[test]
fn test_rotational_pattern() {
    let org_chain = SegmentChain::from_points(&[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0]]);
    {
        // Do nothing
        let mut mod_chain = org_chain.clone();
        mod_chain.rotational_pattern([0.0, 1.0], FRAC_PI_2, 0);
        assert_eq!(org_chain, mod_chain);
    }
    {
        // Single repetition
        let mut mod_chain = org_chain.clone();
        mod_chain.rotational_pattern([0.0, 1.0], FRAC_PI_2, 1);

        assert_eq!(mod_chain.num_segments(), 5);
        let pts: Vec<_> = mod_chain.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [0.0, 2.0]);
        assert_ulps_eq!(pts[4], [0.0, 3.0]);
        assert_ulps_eq!(pts[5], [-1.0, 3.0]);
    }
    {
        // Two repetitions
        let mut mod_chain = org_chain.clone();
        mod_chain.rotational_pattern([0.0, 1.0], FRAC_PI_2, 2);

        assert_eq!(mod_chain.num_segments(), 8);
        let pts: Vec<_> = mod_chain.points().collect();
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
    let org_chain = SegmentChain::from_points(&[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0]]);
    {
        // Do nothing
        let mut mod_chain = org_chain.clone();
        mod_chain.translational_pattern([0.0, 1.0], 0);
        assert_eq!(org_chain, mod_chain);
    }
    {
        // Single repetition
        let mut mod_chain = org_chain.clone();
        mod_chain.translational_pattern([0.0, 1.0], 1);

        assert_eq!(mod_chain.num_segments(), 5);
        let pts: Vec<_> = mod_chain.points().collect();
        assert_ulps_eq!(pts[0], [1.0, 1.0]);
        assert_ulps_eq!(pts[1], [2.0, 1.0]);
        assert_ulps_eq!(pts[2], [2.0, 2.0]);
        assert_ulps_eq!(pts[3], [1.0, 2.0]);
        assert_ulps_eq!(pts[4], [2.0, 2.0]);
        assert_ulps_eq!(pts[5], [2.0, 3.0]);
    }
    {
        // Two repetitions
        let mut mod_chain = org_chain.clone();
        mod_chain.translational_pattern([0.0, 1.0], 2);

        assert_eq!(mod_chain.num_segments(), 8);
        let pts: Vec<_> = mod_chain.points().collect();
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
        let path = SegmentChain::from(arc);

        let mut iter = path.polygonize(Polygonizer::PerType {
            arc: SegmentPolygonizer::MaximumAngle(FRAC_PI_2),
            straight: SegmentPolygonizer::NumberSegments(1),
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
        let path = SegmentChain::from(arc);

        let mut iter = path.polygonize(Polygonizer::PerType {
            arc: SegmentPolygonizer::MaximumAngle(FRAC_PI_2),
            straight: SegmentPolygonizer::NumberSegments(1),
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
fn test_straight_line_chain() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
    let chain = SegmentChain::from_points(&vertices);

    // Points which are contained in the first segment
    assert!(chain.contains_point([0.5, 0.0], 0.0, 0));
    assert!(!chain.contains_point([0.5, 0.1], 0.0, 0));
    assert!(chain.contains_point([0.5, 0.1], 0.11, 0));
    assert!(!chain.contains_point([0.5, 0.2], 0.11, 0));

    // Points which are contained in the second segment
    assert!(chain.contains_point([1.0, 0.5], 0.0, 0));
    assert!(!chain.contains_point([1.1, 0.5], 0.0, 0));
    assert!(chain.contains_point([1.1, 0.5], 0.11, 0));
    assert!(!chain.contains_point([1.2, 0.5], 0.11, 0));

    // Points which are contained in both segments
    assert!(chain.contains_point([1.0, 0.0], 0.0, 0));
    assert!(!chain.contains_point([0.98, 0.02], 0.0, 0));
    assert!(chain.contains_point([0.98, 0.02], 0.11, 0));
    assert!(!chain.contains_point([1.1, -0.1], 0.11, 0));
}

#[test]
fn test_intersection_segment_chains() {
    let first = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let second = SegmentChain::from_points(&[[0.0, 0.5], [2.0, 0.5]]);

    let vec: Vec<_> = first
        .intersections_segment_chain(&second, 0.0, 0)
        .collect();
    assert_eq!(vec.len(), 1);

    let i = vec[0];
    assert_eq!(
        i,
        Intersection {
            point: [1.0, 0.5],
            left: SegmentIdx(1),
            right: SegmentIdx(0)
        }
    );
}

#[test]
fn test_intersection_line_segment_chain() {
    let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let chain = SegmentChain::from_points(&vertices);

    let line = Line::from_point_angle([0.5, 0.5], 0.0);
    let mut intersections = chain.intersections_primitive(&line, 0.0, 0);

    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {
            point: [1.0, 0.5],
            left: SegmentIdx(1),
            right: ()
        })
    );
    approx::assert_abs_diff_eq!(
        intersections.next(),
        Some(Intersection {
            point: [0.0, 0.5],
            left: SegmentIdx(3),
            right: ()
        })
    );
    assert!(intersections.next().is_none());
}
