```rust
use planar_geo::prelude::*;
use std::f64::consts::{PI, FRAC_PI_2};
use approx;

// For brevity
let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

// Construct an arc, a line and another arc using various constructors. These
// constructors can fail for invalid input data, see the expect() strings.
let first_arc =
    ArcSegment::from_center_radius_start_offset_angle([1.5, 0.0], 1.5, PI, -FRAC_PI_2, e, m)
        .expect("radius is positive and offset angle is not zero");
let line = LineSegment::new([1.5, 1.5], [3.5, 1.5], e, m)
    .expect("segment length is not zero");
let second_arc = ArcSegment::from_start_center_angle([3.5, 1.5], [3.5, 0.0], -FRAC_PI_2, e, m)
    .expect("radius is positive and offset angle is not zero");

// Build a segment chain from these three segments
let mut chain = SegmentChain::new();

// Initially, the chain is empty (has no segments)
assert_eq!(chain.num_segments(), 0);

// Now add the three segments
chain.push_back(first_arc.into());
chain.push_back(line.into());
chain.push_back(second_arc.into());
assert_eq!(chain.num_segments(), 3);

// Create a contour out of the chain. If start and end of the chain are not
// identical, a line segment is automatically added.
let contour = Contour::new(chain);

// During the conversion to a contour, a line segment has been added which closes the chain
assert_eq!(contour.num_segments(), 4);

// Create a second contour by creating multiple line segments directly from vertices.
let hole = Contour::new(SegmentChain::from_points(&[[1.5, 0.2], [3.5, 0.2], [3.5, 1.3], [1.5, 1.3]]));

// First element of the vector is interpreted as outer contour, all further
// elements are holes. The resulting shape is visualized below the code snippet.
let shape = Shape::new(vec![contour, hole]).expect("inputs form a valid shape");

// Calculate the surface area: Area of outer contour minus hole area
let quarter_circle_area = 0.25 * PI * (1.5f64).powi(2);
let center_rect_area = 1.5 * 2.0;
let hole_area = 1.1 * 2.0;

// Exact match except for floating point rounding errors (arcs are not
// approximated as polylines)
approx::assert_abs_diff_eq!(shape.area(), 2.0 * quarter_circle_area + center_rect_area - hole_area, epsilon = 1e-15);

// Length of the hole contour
approx::assert_abs_diff_eq!(shape.holes()[0].length(), 2.0 * (1.1 + 2.0), epsilon = 1e-15);

// Length of the first arc -> Is the first segment of the contour
// 0.5 * 1.5 * PI = Quarter circle circumference times radius.
approx::assert_abs_diff_eq!(shape.contour()[0].length(), 0.5 * 1.5 * PI, epsilon = 1e-15); 

// Centroid of the hole
approx::assert_abs_diff_eq!(shape.holes()[0].centroid(), [2.5, 0.75], epsilon = 1e-15);
```

## Transformation

The [`Transformation`] trait allows translating, rotating, scaling and mirroring
of all [`Primitive`] and [`Composite`] types:

```rust
use planar_geo::prelude::*;
use std::f64::consts::PI;

// For brevity
let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

// Translate a point
let mut pt = [2.0, 2.0];
pt.translate([1.0, 2.0]);

assert_eq!(pt, [3.0, 4.0]);

// Rotate a line segment
let mut ls = LineSegment::new([1.0, 0.0], [2.0, 0.0], e, m).expect("points are not equal");
ls.rotate([1.0, 1.0], PI);

approx::assert_abs_diff_eq!(ls.start(), [1.0, 2.0], epsilon = 1e-15);
approx::assert_abs_diff_eq!(ls.stop(), [0.0, 2.0], epsilon = 1e-15);

// Scale a segment chain
let mut chain = SegmentChain::from_points(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
chain.scale(3.0);

let mut pts = chain.points();
assert_eq!(pts.next(), Some([0.0, 0.0]));
assert_eq!(pts.next(), Some([3.0, 0.0]));
assert_eq!(pts.next(), Some([3.0, 3.0]));

// Mirror a contour
let mut contour = Contour::new(chain);
contour.line_reflection([0.0, 0.0], [0.0, 1.0]);

let mut pts = contour.points();
assert_eq!(pts.next(), Some([0.0, 0.0]));
assert_eq!(pts.next(), Some([-3.0, 0.0]));
assert_eq!(pts.next(), Some([-3.0, 3.0]));
```

## Intersections

A major feature of this crate are the various methods available to find
collisions and intersections between different geometric types. 

For example, the following code shows intersections between the segments shown
in this image: