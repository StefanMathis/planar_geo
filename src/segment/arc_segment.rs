/*!
A module containing the [`ArcSegment`] struct, one of the available variants
for the [`Segment`](crate::segment::Segment) enum.

The [`ArcSegment`] is one of the "primitive" geometric types defined in this
crate (implementing the [`Primitive`] type). Wrapped in the
[`Segment`](crate::segment::Segment) enum, it serves as a fundamental building
block of composite geometric types such as the
[`SegmentChain`](crate::segment_chain::SegmentChain).
See the module documentation of [segment](crate::segment) for more.

Most users should interact with this module through the [`ArcSegment`] type
itself; see its documentation for details on construction, invariants, and
usage.
*/

use compare_variables::compare_variables;

use approx::ulps_eq;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    Rotation2, Transformation,
    primitive::{Primitive, PrimitiveIntersections},
};
use bounding_box::BoundingBox;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

/**
A circular, directed arc between two points with a positive radius and a center
point.

*/
#[cfg_attr(
    docsrs,
    doc = "\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/example_arc_segment.svg \"Arc segment\")"
)]
#[cfg_attr(
    not(docsrs),
    doc = "\n\n![>> Example image missing, copy folder docs from crate root to doc root folder (where index.html is) to display the image <<](../../docs/example_arc_segment.svg)"
)]
/**

The endpoints of an arc can be identical (in this case, the arc segment is a
circle), but the radius must be positive.

Of particular interest are the various constructors available for an
[`ArcSegment`]. These always start with `from_` and can construct an instance
of [`ArcSegment`] from various inputs. Some examples:
- From its center, radius, start and offset angle: [`ArcSegment::from_center_radius_start_offset_angle`].
- From its center, radius, start and stop angle: [`ArcSegment::from_center_radius_start_stop_angle`].
- From three points on the arc: [`ArcSegment::from_start_middle_stop`].
- From its start, center and offset angle: [`ArcSegment::from_start_center_angle`].
These constructors can fail if either the radius becomes non-positive or
infinite (which degenerates the arc segment to a line segment) or if the
offset angle becomes zero. The latter is checked using [`approx::ulps_eq`],
which is why all constructors take an `epsilon` and a `max_ulps` argument. See
the documentation of [approx] for more.

See the [module documentation](crate::segment::arc_segment) for more.

# Serialization and deserialization

When the `serde` feature is enabled, line segments can be serialized and
deserialized using the [`serde`] crate.
 */
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct ArcSegment {
    center: [f64; 2],
    radius: f64,
    start_angle: f64,
    offset_angle: f64,
}

impl ArcSegment {
    /**
    Returns the center of `self`.
     */
    pub fn center(&self) -> [f64; 2] {
        return self.center;
    }

    /**
    Returns the radius of `self`.
     */
    pub fn radius(&self) -> f64 {
        return self.radius;
    }

    /**
    Returns the start angle of `self` in rad.
     */
    pub fn start_angle(&self) -> f64 {
        return self.start_angle;
    }

    /**
    Returns the offset angle of `self` in rad. This is the angle covered by the
    arc segment. Add this value to the start angle to get the stop angle.
     */
    pub fn offset_angle(&self) -> f64 {
        return self.offset_angle;
    }

    /**
    Returns the stop angle of `self` in read. This is the sum of
    [`ArcSegment::start_angle`] and [`ArcSegment::offset_angle`].

    # Examples

    ```
    use planar_geo::segment::ArcSegment;

    let arc = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 1.0, -1.0, 0.0, 0).unwrap();
    assert_eq!(arc.stop_angle(), 0.0);
    ```
     */
    pub fn stop_angle(&self) -> f64 {
        return self.start_angle() + self.offset_angle();
    }

    /**
    Returns whether `self` is mathematically positive. Positive means a
    counter-clockwise rotation around [`ArcSegment::center`] from
    [`ArcSegment::start`] to [`ArcSegment::stop`].
     */
    pub fn is_positive(&self) -> bool {
        return self.offset_angle().is_sign_positive();
    }

    /**
    Creates an [`ArcSegment`] representing a full circle, i.e. where
    [`ArcSegment::start_angle`] is 0 and [`ArcSegment::offset_angle`] is 2 times
    the circle number pi. Fails if radius is not positive.

    # Examples

    ```
    use planar_geo::segment::ArcSegment;

    let circle = ArcSegment::circle([1.0, 1.0], 2.0).unwrap();
    approx::assert_abs_diff_eq!(circle.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(circle.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radii are not positive
    assert!(ArcSegment::circle([1.0, 1.0], 0.0).is_err());
    assert!(ArcSegment::circle([1.0, 1.0], -2.0).is_err());
    ```
     */
    pub fn circle(center: [f64; 2], radius: f64) -> crate::error::Result<Self> {
        return Self::from_center_radius_start_offset_angle(center, radius, 0.0, TAU, 0.0, 0);
    }

    /**
    Creates an [`ArcSegment`] from its `center`, `radius`, `start_angle` and
    `offset_angle`. This fails in the following cases:
    - `radius` is not positive.
    - `offset_angle` is approximately zero (checked with [`approx::ulps_eq`] and
    the provided `epsilon` and `max_ulps`).

    # Examples

    ```
    use std::f64::consts::{FRAC_PI_2, TAU};
    use planar_geo::segment::ArcSegment;

    // Successful creation
    let arc = ArcSegment::from_center_radius_start_offset_angle([1.0, 1.0], 2.0, 0.0, -FRAC_PI_2, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, -1.0], epsilon = 1e-15);

    // Example forming a full circle
    let arc = ArcSegment::from_center_radius_start_offset_angle([1.0, 1.0], 2.0, 0.0, TAU, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radius is not positive
    assert!(ArcSegment::from_center_radius_start_offset_angle([1.0, 1.0], -2.0, 0.0, -FRAC_PI_2, 0.0, 0).is_err());

    // Offset angle is zero
    assert!(ArcSegment::from_center_radius_start_offset_angle([1.0, 1.0], 2.0, 0.0, 0.0, 0.0, 0).is_err());
    ```
     */
    pub fn from_center_radius_start_offset_angle(
        center: [f64; 2],
        radius: f64,
        start_angle: f64,
        offset_angle: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        compare_variables!(0.0 < radius)?;
        if approx::ulps_eq!(offset_angle, 0.0, epsilon = epsilon, max_ulps = max_ulps) {
            let mut c = center;
            c.translate([radius * start_angle.cos(), radius * start_angle.sin()]);
            return Err(crate::error::ErrorType::PointsIdentical { start: c, stop: c }.into());
        }

        let start_angle = start_angle.rem_euclid(TAU);
        let offset_angle = if offset_angle.abs() == TAU {
            TAU * offset_angle.signum()
        } else {
            offset_angle % TAU
        };

        return Ok(ArcSegment {
            center,
            radius,
            start_angle,
            offset_angle,
        });
    }

    /**
    Creates an [`ArcSegment`] from its `center`, `radius`, `start_angle` and
    `stop_angle`. This fails in the following cases:
    - `radius` is not positive.
    - `start_angle` is approximately equal to stop_angle (checked with
    [`approx::ulps_eq`] and the provided `epsilon` and `max_ulps`).

    # Examples

    ```
    use std::f64::consts::{FRAC_PI_2, TAU};
    use planar_geo::segment::ArcSegment;

    // Successful creation
    let arc = ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 0.0, FRAC_PI_2, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, 3.0], epsilon = 1e-15);

    // Example forming a full circle
    let arc = ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 0.0, TAU, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radius is not positive
    assert!(ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], -2.0, 0.0, -FRAC_PI_2, 0.0, 0).is_err());

    // Offset angle is zero
    assert!(ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 1.0, 1.0, 0.0, 0).is_err());
    ```
     */
    pub fn from_center_radius_start_stop_angle(
        center: [f64; 2],
        radius: f64,
        start_angle: f64,
        stop_angle: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        let offset_angle = stop_angle - start_angle;
        return Self::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            epsilon,
            max_ulps,
        );
    }

    /**
    Creates an [`ArcSegment`] from its `start`, `stop` and an arbitrary point
    `middle` located somewhere on the arc between `start` and `stop`. This fails
    in the following cases:
    - `middle` is identical to `start` or `stop`.
    - The three points are collinear (which results in the radius becoming
    infinite).

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the offset angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(arc.center(), [2.0, 1.0], epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.radius(), 1.0, epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.start_angle(), 2.0*std::f64::consts::PI, epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.offset_angle(), -0.5*std::f64::consts::PI, epsilon = DEFAULT_EPSILON);

    // Middle is identical to start or stop
    assert!(ArcSegment::from_start_middle_stop(start, start, stop, 0.0, 0).is_err());
    assert!(ArcSegment::from_start_middle_stop(start, stop, stop, 0.0, 0).is_err());

    // All three points are collinear
    assert!(ArcSegment::from_start_middle_stop([0.0, 0.0], [1.0, 0.0], [2.0, 0.0], 0.0, 0).is_err());
    ```
    */
    pub fn from_start_middle_stop(
        start: [f64; 2],
        middle: [f64; 2],
        stop: [f64; 2],
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        let center = three_point_arc_center(start, middle, stop)?;
        let is_positive = three_point_arc_is_positive(start, middle, stop);

        let start_angle = ((start[1] - center[1]).atan2(start[0] - center[0])).rem_euclid(TAU);
        let stop_angle = ((stop[1] - center[1]).atan2(stop[0] - center[0])).rem_euclid(TAU);
        let offset_angle = calculate_offset_angle(start_angle, stop_angle, is_positive);

        let radius = ((start[0] - center[0]).powi(2) + (start[1] - center[1]).powi(2)).sqrt();
        return ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            epsilon,
            max_ulps,
        );
    }

    /**
    Creates an [`ArcSegment`] from its `start`, its `offset_angle` and its
    `center`. This fails in the following cases:
    - `start` is identical to `center`.
    - `offset_angle` is apprixmately zero.

    ```
    use planar_geo::segment::ArcSegment;

    // Clockwise arc => the offset angle is negative
    let start = [0.0, 2.0];
    let center = [0.0, 0.0];
    let angle = 1.5*std::f64::consts::PI;

    let arc = ArcSegment::from_start_center_angle(start, center, angle, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.radius(), 2.0);

    // start and center are identical
    assert!(ArcSegment::from_start_center_angle(start, start, angle, 0.0, 0).is_err());

    // Offset angle is not zero
    assert!(ArcSegment::from_start_center_angle(start, center, 0.0, 0.0, 0).is_err());
    ```
     */
    pub fn from_start_center_angle(
        start: [f64; 2],
        center: [f64; 2],
        offset_angle: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        let start_angle = (start[1] - center[1]).atan2(start[0] - center[0]);
        let radius = ((start[0] - center[0]).powi(2) + (start[1] - center[1]).powi(2)).sqrt();
        return ArcSegment::from_center_radius_start_offset_angle(
            center,
            radius,
            start_angle,
            offset_angle,
            epsilon,
            max_ulps,
        );
    }

    /**
    Creates an [`ArcSegment`] as a ["fillet"](https://en.wikipedia.org/wiki/Fillet_(mechanics))
    (rounding) of the corner formed by connecting `start` to `corner` and
    `corner` to `stop`. This fails in the following cases:
    - `start`, `corner` and `stop` are collinear
    - `radius` is not positive.
    - the arc does not meet the corner line segments because the `radius` is too
    large.

    # Examples

    ```
    use planar_geo::prelude::ArcSegment;

    let arc = ArcSegment::fillet([0.0, 1.0], [1.0, 1.0], [1.0, 0.0], 0.5, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.radius(), 0.5, epsilon = 1e-6);
    approx::assert_abs_diff_eq!(arc.start(), [0.5, 1.0], epsilon = 1e-6);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, 0.5], epsilon = 1e-6);
    ```
     */
    pub fn fillet(
        start: [f64; 2],
        corner: [f64; 2],
        stop: [f64; 2],
        radius: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> crate::error::Result<Self> {
        // The code within this function is based on:
        // https://stackoverflow.com/questions/24771828/algorithm-for-creating-rounded-corners-in-a-polygon

        // Check if two neighboring points are identical. In that case,
        // the points are also collinear
        if ulps_eq!(start, corner, epsilon = epsilon, max_ulps = max_ulps)
            || ulps_eq!(corner, stop, epsilon = epsilon, max_ulps = max_ulps)
        {
            return Err(crate::error::ErrorType::Collinear([start, corner, stop]).into());
        }

        // Check if the radius is zero or negative
        compare_variables!(radius > 0.0)?;

        let start_angle = (corner[1] - start[1]).atan2(corner[0] - start[0]);
        let stop_angle = (corner[1] - stop[1]).atan2(corner[0] - stop[0]);

        let pos_offset_angle: f64;
        let neg_offset_angle: f64;

        // Values for a counter-clockwise arc
        if stop_angle < start_angle {
            pos_offset_angle = stop_angle + TAU - start_angle;
        } else {
            pos_offset_angle = stop_angle - start_angle;
        }

        // Values for a clockwise arc
        if stop_angle > start_angle {
            neg_offset_angle = start_angle + TAU - stop_angle;
        } else {
            neg_offset_angle = start_angle - stop_angle;
        }

        // Angle between the sides (absolute value must be smaller than 180 degree)
        let mut offset_angle: f64;
        let dir: f64;
        if pos_offset_angle.abs() < PI {
            offset_angle = pos_offset_angle;
            dir = 0.5;
        } else {
            offset_angle = neg_offset_angle;
            dir = -0.5;
        }
        offset_angle = offset_angle % TAU;

        // Segment length between fillet start point and cornerent point
        let mut segment_len = radius / (offset_angle / 2.0).tan().abs();

        // Distance of cornerent point to startious and stop point
        let delta_pc = ((corner[0] - start[0]).powi(2) + (corner[1] - start[1]).powi(2)).sqrt();
        let delta_cn = ((corner[0] - stop[0]).powi(2) + (corner[1] - stop[1]).powi(2)).sqrt();

        // Check if the radius is small enough to produce rounded corners
        let min_dist: f64;
        if delta_pc > delta_cn {
            min_dist = delta_cn;
        } else {
            min_dist = delta_pc;
        }
        let radius = if segment_len > min_dist {
            segment_len = min_dist;
            segment_len * (offset_angle / 2.0).tan()
        } else {
            radius
        };

        // Distance between cornerent point and intersections of the normals at fillet start and end point
        let delta_ci = (radius.powi(2) + segment_len.powi(2)).sqrt();

        // Fillet start and end point coordinates
        let xp = corner[0] - (corner[0] - start[0]) * segment_len / delta_pc;
        let yp = corner[1] - (corner[1] - start[1]) * segment_len / delta_pc;
        let start_fillet = [xp, yp]; // Fillet start point

        let xn = corner[0] - (corner[0] - stop[0]) * segment_len / delta_cn;
        let yn = corner[1] - (corner[1] - stop[1]) * segment_len / delta_cn;
        let stop_fillet = [xn, yn]; // Fillet end point

        // Distance between fillet center and cornerent point
        let delta_x = corner[0] * 2.0 - xp - xn;
        let delta_y = corner[1] * 2.0 - yp - yn;
        let delta_cc = (delta_x.powi(2) + delta_y.powi(2)).sqrt();

        // Fillet center coordinates
        let xc = corner[0] - delta_x * delta_ci / delta_cc;
        let yc = corner[1] - delta_y * delta_ci / delta_cc;

        // Calculate the middle point of the fillet with center and angle
        let xm = (start_angle + offset_angle * dir).cos() * radius + xc;
        let ym = (start_angle + offset_angle * dir).sin() * radius + yc;
        let middle_fillet = [xm, ym];

        return ArcSegment::from_start_middle_stop(
            start_fillet,
            middle_fillet,
            stop_fillet,
            epsilon,
            max_ulps,
        );
    }

    /**
    Returns if `angle` is covered by `self`.
    ```
    use planar_geo::prelude::*;

    let center = [0.0, 1.0];
    let radius = 1.0;
    let start_angle = -4.0*std::f64::consts::PI; // This is mathematically identical to 0.0
    let offset_angle = -2.5*std::f64::consts::PI; // This is mathematically identical to -pi/2
    let arc = ArcSegment::from_center_radius_start_offset_angle(center, radius, start_angle, offset_angle, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    assert!(arc.contains_angle(-2.25*std::f64::consts::PI));
    assert!(!arc.contains_angle(-2.75*std::f64::consts::PI));

    // A circle contains all angles
    let start_angle = -4.0*std::f64::consts::PI; // This is mathematically identical to 0.0
    let offset_angle = -2.0*std::f64::consts::PI; // This is mathematically identical to -pi/2
    let arc = ArcSegment::from_center_radius_start_offset_angle(center, radius, start_angle, offset_angle, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();

    assert!(arc.contains_angle(0.0));
    assert!(arc.contains_angle(12.0));
    ```
     */
    pub fn contains_angle(&self, angle: f64) -> bool {
        let angle = angle % TAU;
        let start_angle = self.start_angle();
        let offset_angle = self.offset_angle();
        let stop_angle = start_angle + offset_angle;

        // Arc is counter-clockwise
        if offset_angle > 0.0 {
            return (angle >= start_angle && angle <= stop_angle)
                || (angle + TAU >= start_angle && angle + TAU <= stop_angle)
                || (angle - TAU >= start_angle && angle - TAU <= stop_angle);
        } else {
            return (angle > stop_angle && angle < start_angle)
                || (angle + TAU >= stop_angle && angle + TAU <= start_angle)
                || (angle - TAU >= stop_angle && angle - TAU <= start_angle);
        }
    }

    /**
    Returns the `start` of `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the offset angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, DEFAULT_EPSILON, DEFAULT_MAX_ULPS).unwrap();
    approx::assert_abs_diff_eq!(arc.start(), start);
    ```
     */
    pub fn start(&self) -> [f64; 2] {
        let mut c = self.center.clone();
        c.translate([
            self.radius * self.start_angle().cos(),
            self.radius * self.start_angle().sin(),
        ]);
        return c;
    }

    /**
    Returns the `stop` of `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the offset angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!(arc.stop(), stop, epsilon = 1e-14);
    ```
     */
    pub fn stop(&self) -> [f64; 2] {
        let mut c = self.center.clone();
        c.translate([
            self.radius * self.stop_angle().cos(),
            self.radius * self.stop_angle().sin(),
        ]);
        return c;
    }

    /**
    Returns the number of points in the segment. For a [`ArcSegment`], this is
    always 3.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the offset angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, 0.0, 0).unwrap();
    assert_eq!(arc.number_points(), 3);
    ```
     */
    pub fn number_points(&self) -> usize {
        return 3;
    }

    /**
    Returns the length of `self`. This is the (absolute) product of
    [`ArcSegment::radius`] and [`ArcSegment::offset_angle`] and represents the
    part of the circumference covered by `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the offset angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop, 0.0, 0).unwrap();
    approx::assert_abs_diff_eq!((arc.radius() * arc.offset_angle()).abs(), arc.length());
    ```
     */
    pub fn length(&self) -> f64 {
        return (self.radius() * self.offset_angle()).abs();
    }

    /**
    Returns the centroid / center of mass of `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    let circle = ArcSegment::circle([0.0, -3.0], 2.0).unwrap();

    // In case of a circle, the centroid is the center
    approx::assert_abs_diff_eq!(circle.centroid(), [0.0, -3.0]);
    ```
     */
    pub fn centroid(&self) -> [f64; 2] {
        return crate::CentroidData::from(self).into();
    }

    /**
    Reverses the direction of `self` (i.e., exchanges start and stop point)

    ```
    use planar_geo::prelude::*;

    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let mut arc = ArcSegment::from_start_middle_stop(start, middle, stop, 0.0, 0).unwrap();
    let start = arc.start();
    let stop = arc.stop();
    arc.reverse();
    assert_eq!(arc.start(), stop);
    assert_eq!(arc.stop(), start);
    ```
     */
    pub fn reverse(&mut self) {
        self.start_angle = self.stop_angle();
        self.offset_angle = -self.offset_angle;
    }

    /**
    Returns a point on the segment defined by its normalized position on it.

    For example, `normalized = 0` returns the start point, `normalized = 1`
    returns the end point and `normalized = 0.5` returns the middle point of the
    segment. The input `normalized` is clamped to [0, 1].

    # Examples

    ```
    use planar_geo::prelude::*;

    let arc = ArcSegment::circle([0.0, 0.0], 2.0).unwrap();

    // Middle point of the segment
    approx::assert_abs_diff_eq!(arc.segment_point(0.25), [0.0, 2.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.segment_point(0.5), [-2.0, 0.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.segment_point(0.75), [0.0, -2.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.segment_point(-1.0), [2.0, 0.0], epsilon = 1e-15); // Start point
    approx::assert_abs_diff_eq!(arc.segment_point(1.5), [2.0, 0.0], epsilon = 1e-15); // Stop point
    ```
     */
    pub fn segment_point(&self, normalized: f64) -> [f64; 2] {
        let normalized = normalized.clamp(0.0, 1.0);
        let angle = self.start_angle() + normalized * self.offset_angle;
        return [
            self.center[0] + self.radius * angle.cos(),
            self.center[1] + self.radius * angle.sin(),
        ];
    }

    /**
    Returns the points of a polygon chain which approximates `self`. The number
    of points is defined by the [`SegmentPolygonizer`] (see its docstring). The
    points are regularily distributed over the segment, which means that two
    subsequent points always have the same euclidian distance from each other.

    # Examples

    ```
    use planar_geo::prelude::*;

    let arc = ArcSegment::circle([0.0, 0.0], 2.0).unwrap();

    let mut iter = arc.polygonize(SegmentPolygonizer::NumberSegments(4));

    approx::assert_abs_diff_eq!(iter.next(), Some(arc.segment_point(0.0)), epsilon = 1e-15);
    approx::assert_abs_diff_eq!(iter.next(), Some(arc.segment_point(0.25)), epsilon = 1e-15);
    approx::assert_abs_diff_eq!(iter.next(), Some(arc.segment_point(0.5)), epsilon = 1e-15);
    approx::assert_abs_diff_eq!(iter.next(), Some(arc.segment_point(0.75)), epsilon = 1e-15);
    approx::assert_abs_diff_eq!(iter.next(), Some(arc.segment_point(1.0)), epsilon = 1e-15);
    assert!(iter.next().is_none());
    ```
     */
    pub fn polygonize<'a>(
        &'a self,
        polygonizer: super::SegmentPolygonizer,
    ) -> super::PolygonPointsIterator<'a> {
        let radius = self.radius();

        let num_segs = match polygonizer {
            super::SegmentPolygonizer::NumberSegments(segs) => segs,
            super::SegmentPolygonizer::MaximumSegmentLength(max_len) => {
                if max_len <= 0.0 {
                    0
                } else {
                    let max_angle = 2.0 * (0.5 * max_len / radius).clamp(-1.0, 1.0).asin();
                    (self.offset_angle().abs() / max_angle).ceil() as usize
                }
            }
            super::SegmentPolygonizer::MaximumAngle(max_angle) => {
                if max_angle <= 0.0 {
                    0
                } else {
                    (self.offset_angle().abs() / max_angle).ceil() as usize
                }
            }
        };

        return super::PolygonPointsIterator {
            index: 0,
            num_segs,
            segment: super::SegmentRef::ArcSegment(self),
        };
    }

    /**
    Returns an iterator over the start and stop point of `self`.

    This is a shorthand for `self.polygonize( Polygonizer::NumberSegments(1))`.

    # Examples

    ```
    use planar_geo::segment::ArcSegment;

    let arc = ArcSegment::from_center_radius_start_offset_angle([0.0, 0.0], 2.0, 1.0, -1.0, 0.0, 0).unwrap();

    let mut iter = arc.points();

    assert_eq!(iter.next(), Some(arc.segment_point(0.0)));
    assert_eq!(iter.next(), Some(arc.segment_point(1.0)));
    assert!(iter.next().is_none());
    ```
     */
    pub fn points<'a>(&'a self) -> super::PolygonPointsIterator<'a> {
        return self.polygonize(super::SegmentPolygonizer::NumberSegments(1));
    }

    /// Algorithm from https://cp-algorithms.com/geometry/circle-line-intersection.html
    pub(crate) fn intersections_line_circle(
        &self,
        a: f64,
        b: f64,
        c: f64,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        let x0 = -a * c / (a.powi(2) + b.powi(2));
        let y0 = -b * c / (a.powi(2) + b.powi(2));

        // Argument for the sqrt function
        let arg = self.radius().powi(2) - c.powi(2) / (a.powi(2) + b.powi(2));

        let mut intersections = PrimitiveIntersections::Zero;

        if arg < 0.0 {
            // If the argument is negative, no intersection exists
        } else if arg == 0.0 {
            // Special case: Argument of the sqrt-function is exactly zero => point equals center + (x0, y0)
            let mut pt = self.center();
            pt.translate([x0, y0]);
            if self.contains_point(pt, epsilon, max_ulps) {
                intersections.push(pt);
            }
        } else {
            // Intersects exist
            let m = (arg / (a.powi(2) + b.powi(2))).sqrt();

            // Calculate the two possible points
            let xa = x0 + b * m;
            let ya = y0 - a * m;
            let xb = x0 - b * m;
            let yb = y0 + a * m;

            // Check if the points are actually in the arc. If a point is not in the
            // arc, return None
            let mut pta = self.center();
            pta.translate([xa, ya]);

            let mut ptb = self.center();
            ptb.translate([xb, yb]);

            let angle_a = ya.atan2(xa);
            let angle_b = yb.atan2(xb);

            if self.contains_angle(angle_a) {
                intersections.push(pta);
            }
            if self.contains_angle(angle_b) {
                intersections.push(ptb);
            }
        }
        return intersections;
    }
}

impl crate::primitive::private::Sealed for ArcSegment {}

impl Primitive for ArcSegment {
    fn contains_point(&self, point: [f64; 2], epsilon: f64, max_ulps: u32) -> bool {
        if ulps_eq!(self.start(), point, epsilon = epsilon, max_ulps = max_ulps) {
            return true;
        }
        if ulps_eq!(self.stop(), point, epsilon = epsilon, max_ulps = max_ulps) {
            return true;
        }

        // Quick first check: If point p is outside the segment bounding box, it can't be part of the arc segment.
        if !BoundingBox::from(self).approx_contains_point(point, epsilon, max_ulps) {
            return false;
        }

        let shifted_pt = [point[0] - self.center()[0], point[1] - self.center()[1]];

        if ulps_eq!(
            self.radius().powi(2),
            shifted_pt[0].powi(2) + shifted_pt[1].powi(2),
            epsilon = epsilon,
            max_ulps = max_ulps
        ) {
            let angle = shifted_pt[1].atan2(shifted_pt[0]);
            return self.contains_angle(angle);
        } else {
            return false;
        }
    }

    fn intersections_line(
        &self,
        line: &crate::line::Line,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        line.intersections_arc_segment(self, epsilon, max_ulps)
    }

    fn intersections_line_segment(
        &self,
        line_segment: &crate::segment::LineSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        line_segment.intersections_arc_segment(self, epsilon, max_ulps)
    }

    fn intersections_arc_segment(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        // If the segments are identical, they do not intersect by definition
        if std::ptr::eq(self, arc_segment) {
            return PrimitiveIntersections::Zero;
        }

        // Algorithm from https://cp-algorithms.com/geometry/circle-circle-intersection.html
        let radius_arc1 = self.radius();
        let radius_arc2 = arc_segment.radius();
        let center_arc1 = self.center();
        let center_arc2 = arc_segment.center();

        // Assume the first circle is in the origin and shift the second circle accordingly.
        // The new center coordinates of arc 2 / circle 2 therefore are:
        let x2 = center_arc2[0] - center_arc1[0];
        let y2 = center_arc2[1] - center_arc1[1];

        // Represent the first circle as a line
        let a = -2.0 * x2;
        let b = -2.0 * y2;
        let c = x2.powi(2) + y2.powi(2) + (radius_arc1).powi(2) - (radius_arc2).powi(2);

        // Check if the two circles are actually identical. This is the case, if
        // c == 0. In that case, return the "common" end points (if there are
        // any)
        if c == 0.0 {
            let mut intersections = PrimitiveIntersections::Zero;
            if self.contains_angle(arc_segment.start_angle()) {
                intersections.push(arc_segment.start());
            }
            if self.contains_angle(arc_segment.stop_angle()) {
                intersections.push(arc_segment.stop());
            }
            if intersections.len() == 2 {
                return intersections;
            }
            if arc_segment.contains_angle(self.start_angle()) {
                intersections.push(self.start());
            }
            if intersections.len() == 2 {
                return intersections;
            }
            if arc_segment.contains_angle(self.stop_angle()) {
                intersections.push(self.stop());
            }
            return intersections;
        } else {
            return self.intersections_line_circle(a, b, c, epsilon, max_ulps);
        }
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_ulps: u32,
    ) -> PrimitiveIntersections {
        other.intersections_arc_segment(self, epsilon, max_ulps)
    }
}

impl Transformation for ArcSegment {
    fn translate(&mut self, shift: [f64; 2]) {
        self.center.translate(shift);
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let t = Rotation2::new(angle);
        let pt = [self.center[0] - center[0], self.center[1] - center[1]];
        self.center = t * pt;
        self.center.translate([center[0], center[1]]);

        self.start_angle = (self.start_angle + angle).rem_euclid(TAU);
    }

    fn scale(&mut self, factor: f64) {
        self.center.scale(factor);
        self.radius *= factor;
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        // Treatment of special case: vertical line
        if start[0] == stop[0] {
            self.center = [-self.center[0] + 2.0 * start[0], self.center[1]];

        // Treatment of special case: Horizontal line
        } else if start[1] == stop[1] {
            self.center = [self.center[0], -self.center[1] + 2.0 * start[1]];

        // All other cases
        } else {
            // Solve the line equation
            let m = (stop[1] - start[1]) / (stop[0] - start[0]);
            let c = (stop[0] * start[1] - start[0] * stop[1]) / (stop[0] - start[0]);

            let d = (self.center[0] + (self.center[1] - c) * m) / (1.0 + m.powi(2));
            self.center = [
                2.0 * d - self.center[0],
                2.0 * d * m - self.center[1] + 2.0 * c,
            ];
        };

        // Offset angle is simply inverted
        self.offset_angle = -self.offset_angle;

        // Calculate the angle between start and stop
        let line_angle = (stop[1] - start[1]).atan2(stop[0] - start[0]);
        let diff = line_angle - self.start_angle;
        self.start_angle = (self.start_angle + 2.0 * diff).rem_euclid(TAU);
    }
}

impl From<&ArcSegment> for BoundingBox {
    fn from(value: &ArcSegment) -> BoundingBox {
        let c = value.center();

        let mut xmin_pt = c.clone();
        xmin_pt.translate([-value.radius(), 0.0]);

        let mut xmax_pt = c.clone();
        xmax_pt.translate([value.radius(), 0.0]);

        let mut ymin_pt = c.clone();
        ymin_pt.translate([0.0, -value.radius()]);

        let mut ymax_pt = c.clone();
        ymax_pt.translate([0.0, value.radius()]);

        let start = value.start();
        let stop = value.stop();

        let xmin = if value.contains_angle(PI) {
            xmin_pt[0]
        } else {
            // Smaller of the two values
            if start[0] < stop[0] {
                start[0]
            } else {
                stop[0]
            }
        };

        let xmax = if value.contains_angle(0.0) {
            xmax_pt[0]
        } else {
            // Larger of the two values
            if start[0] > stop[0] {
                start[0]
            } else {
                stop[0]
            }
        };

        let ymin = if value.contains_angle(-FRAC_PI_2) {
            ymin_pt[1]
        } else {
            // Smaller of the two values
            if start[1] < stop[1] {
                start[1]
            } else {
                stop[1]
            }
        };

        let ymax = if value.contains_angle(FRAC_PI_2) {
            ymax_pt[1]
        } else {
            // Larger of the two values
            if start[1] > stop[1] {
                start[1]
            } else {
                stop[1]
            }
        };

        return BoundingBox::new(xmin, xmax, ymin, ymax);
    }
}

impl From<&ArcSegment> for crate::CentroidData {
    fn from(value: &ArcSegment) -> Self {
        let inv_3 = 1.0 / 3.0;

        let center = value.center();
        let angle = value.offset_angle();
        let radius = value.radius();

        // centroid and area for shape 1
        let x = (value.start()[0] + value.stop()[0]) * inv_3;
        let y = (value.start()[1] + value.stop()[1]) * inv_3;
        let area = 0.5 * (value.start()[0] * value.stop()[1] - value.stop()[0] * value.start()[1]);
        let data_shape_1 = crate::CentroidData { area, x, y };

        // centroid and area for shape 2
        // Centroid of a circular segment is calculated by moving the segment so that:
        // a) Center equals origin
        // b) circular segment is symmetric around the y-axis
        // Calculate the centroid and reverse the operation.
        let half_angle = 0.5 * angle;
        let area = half_angle * radius.powi(2);
        let centroid_transformed = [2.0 / 3.0 * radius * half_angle.sin() / half_angle, 0.0];

        let middle_angle = value.start_angle() + 0.5 * value.offset_angle();
        let t = Rotation2::new(middle_angle);
        let mut centroid_movec_back = t * centroid_transformed;
        centroid_movec_back.translate([center[0], center[1]]);

        let data_shape_2 = Self {
            area,
            x: centroid_movec_back[0],
            y: centroid_movec_back[1],
        };

        // centroid and area for shape 3
        let x = (value.start()[0] + value.stop()[0] + center[0]) * inv_3;
        let y = (value.start()[1] + value.stop()[1] + center[1]) * inv_3;
        let area = 0.5
            * ((value.stop()[0] - value.start()[0]) * (center[1] - value.start()[1])
                - (center[0] - value.start()[0]) * (value.stop()[1] - value.start()[1]));

        let data_shape_3 = Self { area, x, y };

        // Apply A = A1 + A2 - A3
        return data_shape_1.union(&data_shape_2).subtract(&data_shape_3);
    }
}

/**
Calculates the center from three points of an arc: `start` point of the arc, an
arbitrary `middle` point somewhere on the path and the arc `stop` point. If the
three points are collinear, the center is undefined and therefore `None` is
returned.

The algorithm is taken from:
<https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points>
```
use planar_geo::prelude::*;

let start = [1.0, 0.0];
let middle = [1.0/2.0_f64.sqrt(), -1.0/2.0_f64.sqrt()];
let stop = [0.0, -1.0];

let res = three_point_arc_center(start, middle, stop);
approx::assert_abs_diff_eq!(res.unwrap(), [0.0, 0.0], epsilon = 1e-15);

let start = [1.0, 0.0];
let middle = [0.5, 0.0];
let stop = [0.0, 0.0];

let res = three_point_arc_center(start, middle, stop);
assert!(res.is_err());
```
 */
pub fn three_point_arc_center(
    start: [f64; 2],
    middle: [f64; 2],
    stop: [f64; 2],
) -> crate::error::Result<[f64; 2]> {
    fn determinant(m: [f64; 9]) -> f64 {
        return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
            + m[2] * (m[3] * m[7] - m[4] * m[6]);
    }

    if geometry_predicates::orient2d(start.into(), middle.into(), stop.into()) == 0.0 {
        return Err(crate::error::ErrorType::Collinear([start, middle, stop]).into());
    }

    let start_sq = start[0].powi(2) + start[1].powi(2);
    let middle_sq = middle[0].powi(2) + middle[1].powi(2);
    let stop_sq = stop[0].powi(2) + stop[1].powi(2);

    let m11_det = determinant([
        start[0], start[1], 1.0, middle[0], middle[1], 1.0, stop[0], stop[1], 1.0,
    ]);
    let m12_det = determinant([
        start_sq, start[1], 1.0, middle_sq, middle[1], 1.0, stop_sq, stop[1], 1.0,
    ]);
    let m13_det = determinant([
        start_sq, start[0], 1.0, middle_sq, middle[0], 1.0, stop_sq, stop[0], 1.0,
    ]);

    let x0 = 0.5 * m12_det / m11_det;
    let y0 = -0.5 * m13_det / m11_det;
    return Ok([x0, y0]);
}

/**
Returns if the arc defined by the three given points is mathematically positive.

```
use planar_geo::prelude::*;

let start = [1.0, 0.0];
let middle = [1.0/2.0_f64.sqrt(), 1.0/2.0_f64.sqrt()];
let stop = [0.0, 1.0];
assert!(three_point_arc_is_positive(start, middle, stop));

let start = [1.0, 0.0];
let middle = [1.0/2.0_f64.sqrt(), -1.0/2.0_f64.sqrt()];
let stop = [0.0, -1.0];
assert!(!three_point_arc_is_positive(start, middle, stop));
```
 */
pub fn three_point_arc_is_positive(start: [f64; 2], middle: [f64; 2], stop: [f64; 2]) -> bool {
    /*
    Cross product of the two vectors a x b:
    a = [stop[0] - start[0], stop[1] - start[1], 0.0]
    b = [middle[0] - start[0], middle[1] - start[1], 0.0]
    Only the last element is compared to zero
    */
    let z = (stop[0] - start[0]) * (middle[1] - start[1])
        - (stop[1] - start[1]) * (middle[0] - start[0]);
    return z < 0.0;
}

/**
Calculates the offset angle from a `start_angle` and a `stop_angle`. The arc
between a start and a stop angle can either be mathematically `positive`
(i.e. counter-clockwise) or negative (not `positive`, i.e. clockwise).
 */
pub fn calculate_offset_angle(start_angle: f64, stop_angle: f64, positive: bool) -> f64 {
    let offset_angle: f64;
    if positive {
        if start_angle < stop_angle {
            offset_angle = stop_angle - start_angle;
        } else {
            offset_angle = TAU - (start_angle - stop_angle);
        }
    } else {
        if start_angle > stop_angle {
            offset_angle = stop_angle - start_angle;
        } else {
            offset_angle = (stop_angle - start_angle) - TAU;
        }
    }
    return offset_angle;
}
