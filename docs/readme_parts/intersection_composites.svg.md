_This image was created with examples/intersection_segments.rs_

```rust
use planar_geo::prelude::*;
use std::f64::consts::PI;

// Abbreviated to make examples more concise
let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

let line_1: Segment = LineSegment::new([0.0, 0.0], [3.0, 0.0], e, m)
    .expect("segment length is not zero").into();
let line_2: Segment = LineSegment::new([2.0, -0.5], [2.0, 0.5], e, m)
    .expect("segment length is not zero").into();
let arc: Segment = ArcSegment::from_center_radius_start_offset_angle(
    [1.0, 0.0],
    0.5,
    -0.25*PI,
    1.5*PI,
    e,
    m,
).expect("offet angle is not zero").into();

// Are the following points part of the respective segment?
let pt1 = [2.3, 0.0];
let pt2 = [2.0, 0.1];
let pt3 = [1.0, 0.5];

assert!(line_1.contains_point(pt1, e, m));
assert!(!line_1.contains_point(pt2, e, m));
assert!(line_2.contains_point(pt2, e, m));
assert!(arc.contains_point(pt3, e, m));

// Find intersections between the segments. The order doesn't matter, i.e.
// line_1.intersections_primitive(&line_2) produces the same result as
// line_2.intersections_primitive(&line_1).

// line_1 and line_2 intersect once in point (2, 0)
assert_eq!(line_1.intersections_primitive(&line_2, e, m), PrimitiveIntersections::One([2.0, 0.0]));

// line_1 and arc intersect twice (in (0.5, 0.0) and (1.5, 0.0))
assert_eq!(line_1.intersections_primitive(&arc, e, m), PrimitiveIntersections::Two([[1.5, 0.0], [0.5, 0.0]]));

// line_2 and arc don't intersect at all
assert_eq!(line_2.intersections_primitive(&arc, e, m), PrimitiveIntersections::Zero);
```

It is also possible to calculate the intersections between composite types, as
shown in `examples/intersection_composites.rs`: