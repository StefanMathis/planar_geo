/*!
A module containing the [`ArcSegment`] struct, one of the available variants
for the [`Segment`](crate::segment::Segment) enum.

The [`ArcSegment`] is one of the "primitive" geometric types defined in this
crate (implementing the [`Primitive`] type). Wrapped in the
[`Segment`](crate::segment::Segment) enum, it serves as a fundamental building
block of composite geometric types such as the
[`Polysegment`](crate::polysegment::Polysegment).
See the module documentation of [segment](crate::segment) for more.

Most users should interact with this module through the [`ArcSegment`] type
itself; see its documentation for details on construction, invariants, and
usage.
*/

use approx::{relative_eq, relative_ne};
use compare_variables::compare_variables;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    CentroidData, Rotation2, Transformation,
    primitive::{Primitive, PrimitiveIntersections},
};
use bounding_box::{BoundingBox, ToBoundingBox};

use std::f64::consts::{FRAC_PI_2, PI, TAU};

/**
A circular, directed arc between two points with a positive radius and a center
point.

*/
#[doc = ""]
#[cfg_attr(
    feature = "doc-images",
    doc = "![Arc segment example][example_arc_segment]"
)]
#[cfg_attr(
    feature = "doc-images",
    embed_doc_image::embed_doc_image("example_arc_segment", "docs/img/example_arc_segment.svg")
)]
#[cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images not enabled**. Compile docs with
    `cargo doc --features 'doc-images'` and Rust version >= 1.54."
)]
/**

The endpoints of an arc can be identical (in this case, the arc segment is a
circle), but the radius must be positive.

Of particular interest are the various constructors available for an
[`ArcSegment`]. These always start with `from_` and can construct an
[`ArcSegment`] from various inputs. Some examples:
- From its center, radius, start angle and sweep angle: [`ArcSegment::from_center_radius_start_sweep_angle`].
- From its center, radius, start angle and stop angle: [`ArcSegment::from_center_radius_start_stop_angle`].
- From three points on the arc: [`ArcSegment::from_start_middle_stop`].
- From its start, center and sweep angle: [`ArcSegment::from_start_center_angle`].
- From its radius, start and stop points: [`ArcSegment::from_start_stop_radius`].

These constructors can fail if either the radius becomes non-positive or
infinite (which degenerates the arc segment to a line segment) or if the
sweep angle becomes zero. The latter is checked using [`approx::relative_eq`],
which is why all constructors take an `epsilon` and a `max_relative` argument. See
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
    sweep_angle: f64,
}

impl ArcSegment {
    /**
    Creates an [`ArcSegment`] from its `center`, `radius`, `start_angle` and
    `sweep_angle`. This fails in the following cases:
    - `radius` is not positive.
    - `sweep_angle` is approximately zero (checked with [`approx::relative_eq`] and
    the provided `epsilon` and `max_relative`).

    This function is an alias for
    [`ArcSegment::from_center_radius_start_sweep_angle`].

    # Examples

    ```
    use std::f64::consts::{FRAC_PI_2, TAU};
    use planar_geo::segment::ArcSegment;

    // Successful creation
    let arc = ArcSegment::new([1.0, 1.0], 2.0, 0.0, -FRAC_PI_2).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, -1.0], epsilon = 1e-15);

    // Example forming a full circle
    let arc = ArcSegment::new([1.0, 1.0], 2.0, 0.0, TAU).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radius is not positive
    assert!(ArcSegment::new([1.0, 1.0], -2.0, 0.0, -FRAC_PI_2).is_err());

    // sweep angle is zero
    assert!(ArcSegment::new([1.0, 1.0], 2.0, 0.0, 0.0).is_err());
    ```
     */
    pub fn new(
        center: [f64; 2],
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    ) -> crate::error::Result<Self> {
        return Self::from_center_radius_start_sweep_angle(
            center,
            radius,
            start_angle,
            sweep_angle,
        );
    }

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
    Returns the sweep angle of `self` in rad. This is the angle covered by the
    arc segment. Add this value to the start angle to get the stop angle.
     */
    pub fn sweep_angle(&self) -> f64 {
        return self.sweep_angle;
    }

    /**
    Returns the stop angle of `self` in read. This is the sum of
    [`ArcSegment::start_angle`] and [`ArcSegment::sweep_angle`].

    # Examples

    ```
    use planar_geo::segment::ArcSegment;

    let arc = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 1.0, -1.0).expect("valid_arc");
    assert_eq!(arc.stop_angle(), 0.0);
    ```
     */
    pub fn stop_angle(&self) -> f64 {
        return self.start_angle() + self.sweep_angle();
    }

    /**
    Returns whether `self` is mathematically positive. Positive means a
    counter-clockwise rotation around [`ArcSegment::center`] from
    [`ArcSegment::start`] to [`ArcSegment::stop`].
     */
    pub fn is_positive(&self) -> bool {
        return self.sweep_angle().is_sign_positive();
    }

    /**
    Creates an [`ArcSegment`] representing a full circle, i.e. where
    [`ArcSegment::start_angle`] is 0 and [`ArcSegment::sweep_angle`] is 2 times
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
        return Self::from_center_radius_start_sweep_angle(center, radius, 0.0, TAU);
    }

    /**
    Creates an [`ArcSegment`] from its `center`, `radius`, `start_angle` and
    `sweep_angle`. This fails in the following cases:
    - `radius` is not positive.
    - `sweep_angle` is approximately zero (checked with [`approx::relative_eq`] and
    the provided `epsilon` and `max_relative`).

    # Examples

    ```
    use std::f64::consts::{FRAC_PI_2, TAU};
    use planar_geo::segment::ArcSegment;

    // Successful creation
    let arc = ArcSegment::from_center_radius_start_sweep_angle([1.0, 1.0], 2.0, 0.0, -FRAC_PI_2).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, -1.0], epsilon = 1e-15);

    // Example forming a full circle
    let arc = ArcSegment::from_center_radius_start_sweep_angle([1.0, 1.0], 2.0, 0.0, TAU).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radius is not positive
    assert!(ArcSegment::from_center_radius_start_sweep_angle([1.0, 1.0], -2.0, 0.0, -FRAC_PI_2).is_err());

    // sweep angle is zero
    assert!(ArcSegment::from_center_radius_start_sweep_angle([1.0, 1.0], 2.0, 0.0, 0.0).is_err());
    ```
     */
    pub fn from_center_radius_start_sweep_angle(
        center: [f64; 2],
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    ) -> crate::error::Result<Self> {
        compare_variables!(0.0 < radius)?;
        if sweep_angle == 0.0 {
            return Err(crate::error::ErrorType::ZeroArcSweep.into());
        }

        let start_angle = start_angle.rem_euclid(TAU);
        let sweep_angle = if sweep_angle.abs() == TAU {
            TAU * sweep_angle.signum()
        } else {
            sweep_angle % TAU
        };

        return Ok(ArcSegment {
            center,
            radius,
            start_angle,
            sweep_angle,
        });
    }

    /**
    Creates an [`ArcSegment`] from its `center`, `radius`, `start_angle` and
    `stop_angle`. This fails in the following cases:
    - `radius` is not positive.
    - `start_angle` is approximately equal to `stop_angle` (checked with
    [`approx::relative_eq`] and the provided `epsilon` and `max_relative`).

    # Examples

    ```
    use std::f64::consts::{FRAC_PI_2, TAU};
    use planar_geo::segment::ArcSegment;

    // Successful creation
    let arc = ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 0.0, FRAC_PI_2).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [1.0, 3.0], epsilon = 1e-15);

    // Example forming a full circle
    let arc = ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 0.0, TAU).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [3.0, 1.0], epsilon = 1e-15);
    approx::assert_abs_diff_eq!(arc.stop(), [3.0, 1.0], epsilon = 1e-15);

    // Radius is not positive
    assert!(ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], -2.0, 0.0, -FRAC_PI_2).is_err());

    // sweep angle is zero
    assert!(ArcSegment::from_center_radius_start_stop_angle([1.0, 1.0], 2.0, 1.0, 1.0).is_err());
    ```
     */
    pub fn from_center_radius_start_stop_angle(
        center: [f64; 2],
        radius: f64,
        start_angle: f64,
        stop_angle: f64,
    ) -> crate::error::Result<Self> {
        let sweep_angle = stop_angle - start_angle;
        return Self::from_center_radius_start_sweep_angle(
            center,
            radius,
            start_angle,
            sweep_angle,
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

    // Clockwise arc => the sweep angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.center(), [2.0, 1.0], epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.radius(), 1.0, epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.start_angle(), 2.0*std::f64::consts::PI, epsilon = DEFAULT_EPSILON);
    approx::assert_abs_diff_eq!(arc.sweep_angle(), -0.5*std::f64::consts::PI, epsilon = DEFAULT_EPSILON);

    // Middle is identical to start or stop
    assert!(ArcSegment::from_start_middle_stop(start, start, stop).is_err());
    assert!(ArcSegment::from_start_middle_stop(start, stop, stop).is_err());

    // All three points are collinear
    assert!(ArcSegment::from_start_middle_stop([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]).is_err());
    ```
    */
    pub fn from_start_middle_stop(
        start: [f64; 2],
        middle: [f64; 2],
        stop: [f64; 2],
    ) -> crate::error::Result<Self> {
        let center = three_point_arc_center(start, middle, stop)?;
        let is_positive = three_point_arc_is_positive(start, middle, stop);

        let start_angle = ((start[1] - center[1]).atan2(start[0] - center[0])).rem_euclid(TAU);
        let stop_angle = ((stop[1] - center[1]).atan2(stop[0] - center[0])).rem_euclid(TAU);
        let sweep_angle = calculate_sweep_angle(start_angle, stop_angle, is_positive);

        let radius = ((start[0] - center[0]).powi(2) + (start[1] - center[1]).powi(2)).sqrt();
        return ArcSegment::from_center_radius_start_sweep_angle(
            center,
            radius,
            start_angle,
            sweep_angle,
        );
    }

    /**
    Creates an [`ArcSegment`] from its `start`, its `sweep_angle` and its
    `center`. This fails in the following cases:
    - `start` is identical to `center`.
    - `sweep_angle` is apprixmately zero.

    ```
    use planar_geo::segment::ArcSegment;

    // Clockwise arc => the sweep angle is negative
    let start = [0.0, 2.0];
    let center = [0.0, 0.0];
    let angle = 1.5*std::f64::consts::PI;

    let arc = ArcSegment::from_start_center_angle(start, center, angle).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.radius(), 2.0);

    // start and center are identical
    assert!(ArcSegment::from_start_center_angle(start, start, angle).is_err());

    // sweep angle is zero
    assert!(ArcSegment::from_start_center_angle(start, center, 0.0).is_err());
    ```
     */
    pub fn from_start_center_angle(
        start: [f64; 2],
        center: [f64; 2],
        sweep_angle: f64,
    ) -> crate::error::Result<Self> {
        let start_angle = (start[1] - center[1]).atan2(start[0] - center[0]);
        let radius = ((start[0] - center[0]).powi(2) + (start[1] - center[1]).powi(2)).sqrt();
        return ArcSegment::from_center_radius_start_sweep_angle(
            center,
            radius,
            start_angle,
            sweep_angle,
        );
    }

    /// Creates an [`ArcSegment`] from its `start` point, `stop` point and
    /// `radius`.
    ///
    /// This fails if `start` is identical to `stop`. If the radius is
    /// smaller than half the euclidian distance between `start` and `stop`, it
    /// the radius is set to that value in order to avoid rejection of valid
    /// arc segments due to rounding errors.
    ///   
    /// Four different arc segments can be constructed from the input
    /// parameters. Therefore, it is also necessary to specify the
    /// arc direction (`positive` means mathematically positive /
    /// counter-clockwise) and whether to pick the `large_arc` or not. The
    /// example below shows all four possibilities:
    ///
    /// ```
    /// use planar_geo::segment::ArcSegment;
    /// use std::f64::consts::PI;
    ///
    /// let start = [0.0, 2.0];
    /// let stop = [2.0, 0.0];
    /// let radius: f64 = 2.0;
    ///
    /// let red_arc = ArcSegment::from_start_stop_radius(start, stop, radius, true, false).expect("valid_arc");
    /// approx::assert_abs_diff_eq!(red_arc.center(), [2.0, 2.0]);
    /// approx::assert_abs_diff_eq!(red_arc.sweep_angle(), 0.5 * PI, epsilon = 1e-3);
    ///
    /// let blue_arc = ArcSegment::from_start_stop_radius(start, stop, radius, true, true).expect("valid_arc");
    /// approx::assert_abs_diff_eq!(blue_arc.center(), [0.0, 0.0]);
    /// approx::assert_abs_diff_eq!(blue_arc.sweep_angle(), 1.5 * PI, epsilon = 1e-3);
    ///
    /// let green_arc = ArcSegment::from_start_stop_radius(start, stop, radius, false, false).expect("valid_arc");
    /// approx::assert_abs_diff_eq!(green_arc.center(), [0.0, 0.0]);
    /// approx::assert_abs_diff_eq!(green_arc.sweep_angle(), -0.5 * PI, epsilon = 1e-3);
    ///
    /// let yellow_arc = ArcSegment::from_start_stop_radius(start, stop, radius, false, true).expect("valid_arc");
    /// approx::assert_abs_diff_eq!(yellow_arc.center(), [2.0, 2.0]);
    /// approx::assert_abs_diff_eq!(yellow_arc.sweep_angle(), -1.5 * PI, epsilon = 1e-3);
    /// ```
    #[doc = ""]
    #[cfg_attr(feature = "doc-images", doc = "![Four possible arcs][four_arcs]")]
    #[cfg_attr(
        feature = "doc-images",
        embed_doc_image::embed_doc_image("four_arcs", "docs/img/example_four_arcs.svg")
    )]
    #[cfg_attr(
        not(feature = "doc-images"),
        doc = "**Doc images not enabled**. Compile docs with `cargo doc --features 'doc-images'` and Rust version >= 1.54."
    )]
    pub fn from_start_stop_radius(
        start: [f64; 2],
        stop: [f64; 2],
        radius: f64,
        positive: bool,
        large_arc: bool,
    ) -> crate::error::Result<Self> {
        fn compute_sweep_angle(
            start: [f64; 2],
            stop: [f64; 2],
            center: [f64; 2],
            positive: bool,
        ) -> (f64, f64) {
            let start_angle = (start[1] - center[1]).atan2(start[0] - center[0]);
            let stop_angle = (stop[1] - center[1]).atan2(stop[0] - center[0]);
            let mut delta = stop_angle - start_angle;
            if positive {
                if delta < 0.0 {
                    delta += TAU;
                }
            } else {
                if delta > 0.0 {
                    delta -= TAU;
                }
            }
            (start_angle, delta)
        }

        compare_variables!(radius > 0.0)?;

        let (x1, y1) = (start[0], start[1]);
        let (x2, y2) = (stop[0], stop[1]);

        // Vector from start to stop
        let dx = x2 - x1;
        let dy = y2 - y1;
        let distance_start_stop = (dx * dx + dy * dy).sqrt();

        // If the distance is zero, start and stop are identical
        if distance_start_stop == 0.0 {
            return Err(crate::error::ErrorType::PointsIdentical.into());
        }

        //
        let mut radius = radius;
        if 0.5 * distance_start_stop > radius {
            radius = 0.5 * distance_start_stop
        }

        // --- Midpoint ---
        let mx = (x1 + x2) * 0.5;
        let my = (y1 + y2) * 0.5;

        // --- Distance from midpoint to center ---
        let half_d = distance_start_stop * 0.5;
        let h = (radius * radius - half_d * half_d).sqrt();

        // --- Normalized perpendicular vector ---
        let inv_d = 1.0 / distance_start_stop;
        let px = -dy * inv_d;
        let py = dx * inv_d;

        // Two possible centers
        let c1 = [mx + px * h, my + py * h];
        let c2 = [mx - px * h, my - py * h];

        // Choose correct center based on orientation
        let (start_angle_1, sweep_angle_1) = compute_sweep_angle(start, stop, c1, positive);
        let (start_angle_2, sweep_angle_2) = compute_sweep_angle(start, stop, c2, positive);

        let is_large1 = sweep_angle_1.abs() > std::f64::consts::PI;
        let is_large2 = sweep_angle_2.abs() > std::f64::consts::PI;

        let (center, start_angle, sweep_angle) =
            match (is_large1 == large_arc, is_large2 == large_arc) {
                (true, false) => (c1, start_angle_1, sweep_angle_1),
                (false, true) => (c2, start_angle_2, sweep_angle_2),
                (true, true) => {
                    // This can happen numerically near PI — pick one consistently
                    (c1, start_angle_1, sweep_angle_1)
                }
                (false, false) => {
                    // Numerical issue near π — pick the closer one
                    let target = if large_arc {
                        std::f64::consts::PI * 1.5
                    } else {
                        std::f64::consts::PI * 0.5
                    };

                    // Calculate the deviation from the "target" arc
                    let err1 = (sweep_angle_1.abs() - target).abs();
                    let err2 = (sweep_angle_2.abs() - target).abs();

                    if err1 < err2 {
                        (c1, start_angle_1, sweep_angle_1)
                    } else {
                        (c2, start_angle_2, sweep_angle_2)
                    }
                }
            };

        Self::from_center_radius_start_sweep_angle(center, radius, start_angle, sweep_angle)
    }

    /**
    Creates an [`ArcSegment`] from its `start`, `stop` and `center`.

    TODO: Explain that this constructor is corrective!

    With those three input parameters, the `radius` of the arc is overdefined,
    since it can be calculated both from the distance `start` - `center` and
    `stop` - `center`. Hence, this function calculates both distances and
    compares them. If they are (approximately) equal, this function then
    forwards to [`ArcSegment::from_start_stop_radius`]. See the docstring of
    that function for an explanation of the `positive` and `large_arc`
    arguments.

    # Examples

    ```
    use planar_geo::segment::ArcSegment;
    use std::f64::consts::PI;

    let start = [0.0, 2.0];
    let stop = [2.0, 0.0];
    let radius: f64 = 2.0;

    let arc = ArcSegment::from_start_stop_center([0.0, 2.0], [2.0, 0.0], [0.0, 0.0], true, false).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.radius(), 2.0);
    approx::assert_abs_diff_eq!(arc.sweep_angle(), 0.5 * PI, epsilon = 1e-3);

    // Radii start-center and stop-center are not equal
    assert!(ArcSegment::from_start_stop_center([0.0, 2.0], [3.0, 0.0], [0.0, 0.0], true, false).is_err());
    ```
     */
    pub fn from_start_stop_center(
        start: [f64; 2],
        stop: [f64; 2],
        center: [f64; 2],
        positive: bool,
        large_arc: bool,
    ) -> crate::error::Result<Self> {
        let radius_start_center =
            ((start[0] - center[0]).powi(2) + (start[1] - center[1]).powi(2)).sqrt();
        let radius_stop_center =
            ((stop[0] - center[0]).powi(2) + (stop[1] - center[1]).powi(2)).sqrt();

        if relative_ne!(
            radius_start_center,
            radius_stop_center,
            epsilon = crate::DEFAULT_EPSILON,
            max_relative = crate::DEFAULT_MAX_RELATIVE
        ) {
            // We already know this will fail, we are just interested in the
            // produced error message.
            compare_variables!(radius_start_center == radius_stop_center)?;
        }

        return ArcSegment::from_start_stop_radius(
            start,
            stop,
            radius_start_center,
            positive,
            large_arc,
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

    let arc = ArcSegment::fillet([0.0, 1.0], [1.0, 1.0], [1.0, 0.0], 0.5).expect("valid_arc");
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
    ) -> crate::error::Result<Self> {
        // The code within this function is based on:
        // https://stackoverflow.com/questions/24771828/algorithm-for-creating-rounded-corners-in-a-polygon

        // Check if two neighboring points are identical. In that case,
        // the points are also collinear
        if start == corner || corner == stop {
            return Err(crate::error::ErrorType::Collinear([start, corner, stop]).into());
        }

        // Check if the radius is zero or negative
        compare_variables!(radius > 0.0)?;

        let start_angle = (corner[1] - start[1]).atan2(corner[0] - start[0]);
        let stop_angle = (corner[1] - stop[1]).atan2(corner[0] - stop[0]);

        let pos_sweep_angle: f64;
        let neg_sweep_angle: f64;

        // Values for a counter-clockwise arc
        if stop_angle < start_angle {
            pos_sweep_angle = stop_angle + TAU - start_angle;
        } else {
            pos_sweep_angle = stop_angle - start_angle;
        }

        // Values for a clockwise arc
        if stop_angle > start_angle {
            neg_sweep_angle = start_angle + TAU - stop_angle;
        } else {
            neg_sweep_angle = start_angle - stop_angle;
        }

        // Angle between the sides (absolute value must be smaller than 180 degree)
        let mut sweep_angle: f64;
        let dir: f64;
        if pos_sweep_angle.abs() < PI {
            sweep_angle = pos_sweep_angle;
            dir = 0.5;
        } else {
            sweep_angle = neg_sweep_angle;
            dir = -0.5;
        }
        sweep_angle = sweep_angle % TAU;

        // Segment length between fillet start point and cornerent point
        let mut segment_len = radius / (sweep_angle / 2.0).tan().abs();

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
            segment_len * (sweep_angle / 2.0).tan()
        } else {
            radius
        };

        // Distance between cornerent point and intersections of the normals at fillet
        // start and end point
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
        let xm = (start_angle + sweep_angle * dir).cos() * radius + xc;
        let ym = (start_angle + sweep_angle * dir).sin() * radius + yc;
        let middle_fillet = [xm, ym];

        return ArcSegment::from_start_middle_stop(start_fillet, middle_fillet, stop_fillet);
    }

    /**
    Returns whether `self` is a full circle (its `sweep_angle` covers a full
    2*pi rad)
     */
    pub fn is_circle(&self) -> bool {
        return self.sweep_angle().abs() >= TAU;
    }

    /**
    Returns if `angle` is covered by `self`.
    ```
    use planar_geo::prelude::*;

    let center = [0.0, 1.0];
    let radius = 1.0;
    let start_angle = -4.0*std::f64::consts::PI; // This is mathematically identical to 0.0
    let sweep_angle = -2.5*std::f64::consts::PI; // This is mathematically identical to -pi/2
    let arc = ArcSegment::from_center_radius_start_sweep_angle(center, radius, start_angle, sweep_angle).expect("valid_arc");

    assert!(arc.covers_angle(-2.25*std::f64::consts::PI));
    assert!(!arc.covers_angle(-2.75*std::f64::consts::PI));

    // A circle contains all angles
    let start_angle = -4.0*std::f64::consts::PI; // This is mathematically identical to 0.0
    let sweep_angle = -2.0*std::f64::consts::PI; // This is mathematically identical to -pi/2
    let arc = ArcSegment::from_center_radius_start_sweep_angle(center, radius, start_angle, sweep_angle).expect("valid_arc");

    assert!(arc.covers_angle(0.0));
    assert!(arc.covers_angle(12.0));
    ```
     */
    pub fn covers_angle(&self, angle: f64) -> bool {
        let angle = angle % TAU;
        let start_angle = self.start_angle();
        let sweep_angle = self.sweep_angle();
        let stop_angle = start_angle + sweep_angle;

        // Arc is counter-clockwise
        if sweep_angle > 0.0 {
            return (angle >= start_angle && angle <= stop_angle)
                || (angle + TAU >= start_angle && angle + TAU <= stop_angle)
                || (angle - TAU >= start_angle && angle - TAU <= stop_angle);
        } else {
            return (angle >= stop_angle && angle <= start_angle)
                || (angle + TAU >= stop_angle && angle + TAU <= start_angle)
                || (angle - TAU >= stop_angle && angle - TAU <= start_angle);
        }
    }

    /**
    Returns the `start` of `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the sweep angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
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

    // Clockwise arc => the sweep angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
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

    // Clockwise arc => the sweep angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
    assert_eq!(arc.number_points(), 3);
    ```
     */
    pub fn number_points(&self) -> usize {
        return 3;
    }

    /**
    Returns the length of `self`. This is the (absolute) product of
    [`ArcSegment::radius`] and [`ArcSegment::sweep_angle`] and represents the
    part of the circumference covered by `self`.

    # Examples

    ```
    use planar_geo::prelude::*;

    // Clockwise arc => the sweep angle is negative
    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
    approx::assert_abs_diff_eq!((arc.radius() * arc.sweep_angle()).abs(), arc.length());
    ```
     */
    pub fn length(&self) -> f64 {
        return (self.radius() * self.sweep_angle()).abs();
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
        return CentroidData::from(self).into();
    }

    /**
    Reverses the direction of `self` (i.e., exchanges start and stop point)

    ```
    use planar_geo::prelude::*;

    let start = [3.0, 1.0];
    let middle = [2.0 + 1.0/2.0_f64.sqrt(), 1.0 - 1.0/2.0_f64.sqrt()];
    let stop = [2.0, 0.0];

    let mut arc = ArcSegment::from_start_middle_stop(start, middle, stop).expect("valid_arc");
    let start = arc.start();
    let stop = arc.stop();
    arc.reverse();
    assert_eq!(arc.start(), stop);
    assert_eq!(arc.stop(), start);
    ```
     */
    pub fn reverse(&mut self) {
        self.start_angle = self.stop_angle();
        self.sweep_angle = -self.sweep_angle;
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
        let angle = self.start_angle() + normalized * self.sweep_angle;
        return [
            self.center[0] + self.radius * angle.cos(),
            self.center[1] + self.radius * angle.sin(),
        ];
    }

    /**
    Returns the points of a polygon polysegment which approximates `self`. The number
    of points is defined by the
    [`SegmentPolygonizer`](crate::segment::SegmentPolygonizer) (see its
    docstring). The points are regularily distributed over the segment, which
    means that two subsequent points always have the same euclidian distance
    from each other.

    # Examples

    ```
    use planar_geo::prelude::*;

    let arc = ArcSegment::circle([0.0, 0.0], 2.0).unwrap();

    let mut iter = arc.polygonize(SegmentPolygonizer::InnerSegments(4));

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
            super::SegmentPolygonizer::InnerSegments(segs) => segs,
            super::SegmentPolygonizer::MaximumSegmentLength(max_len) => {
                if max_len <= 0.0 {
                    0
                } else {
                    let max_angle = 2.0 * (0.5 * max_len / radius).clamp(-1.0, 1.0).asin();
                    (self.sweep_angle().abs() / max_angle).ceil() as usize
                }
            }
            super::SegmentPolygonizer::MaximumAngle(max_angle) => {
                if max_angle <= 0.0 {
                    0
                } else {
                    (self.sweep_angle().abs() / max_angle).ceil() as usize
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

    This is a shorthand for `self.polygonize( Polygonizer::InnerSegments(1))`.

    # Examples

    ```
    use planar_geo::segment::ArcSegment;

    let arc = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 1.0, -1.0).expect("valid_arc");

    let mut iter = arc.points();

    assert_eq!(iter.next(), Some(arc.segment_point(0.0)));
    assert_eq!(iter.next(), Some(arc.segment_point(1.0)));
    assert!(iter.next().is_none());
    ```
     */
    pub fn points<'a>(&'a self) -> super::PolygonPointsIterator<'a> {
        return self.polygonize(super::SegmentPolygonizer::InnerSegments(1));
    }

    /**
    Switches start and end / stop points of `self`.

    # Examples

    ```
    use std::f64::consts::PI;
    use approx;
    use planar_geo::prelude::*;

    let mut arc = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, PI).expect("valid_arc");
    approx::assert_abs_diff_eq!(arc.start(), [2.0, 0.0], epsilon = 1e-3);
    approx::assert_abs_diff_eq!(arc.stop(), [-2.0, 0.0], epsilon = 1e-3);

    arc.invert();
    approx::assert_abs_diff_eq!(arc.start(), [-2.0, 0.0], epsilon = 1e-3);
    approx::assert_abs_diff_eq!(arc.stop(), [2.0, 0.0], epsilon = 1e-3);

    arc.invert();
    approx::assert_abs_diff_eq!(arc.start(), [2.0, 0.0], epsilon = 1e-3);
    approx::assert_abs_diff_eq!(arc.stop(), [-2.0, 0.0], epsilon = 1e-3);
    ```
     */
    pub fn invert(&mut self) {
        self.start_angle = self.stop_angle();
        self.sweep_angle = -self.sweep_angle;
    }

    /**
    Returns whether `self` and `other` lay on the same circle (center and radius
    are equal).

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let a1 = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, -1.0, 1.0).expect("valid_arc");
    let a2 = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, 1.0).expect("valid_arc");
    let a3 = ArcSegment::from_center_radius_start_sweep_angle([1.0, 0.0], 2.0, 0.0, 1.5 * PI).expect("valid_arc");

    assert!(a1.same_circle(&a2, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!a1.same_circle(&a3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    assert!(!a2.same_circle(&a3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn same_circle(&self, other: &ArcSegment, epsilon: f64, max_relative: f64) -> bool {
        relative_eq!(
            self.center,
            other.center,
            epsilon = epsilon,
            max_relative = max_relative
        ) && relative_eq!(
            self.radius,
            other.radius,
            epsilon = epsilon,
            max_relative = max_relative
        )
    }

    /**
    Returns whether `line_segment` is a tangent of `self`.

    A [`LineSegment`](super::LineSegment) is a tangent of an [`ArcSegment`] if
    they touch in a single point and the line segment is perpendicular to the
    straight line connecting this point and the center of the arc.

    # Examples

    ```
    use std::f64::consts::PI;
    use planar_geo::prelude::*;

    let arc = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, PI).expect("valid_arc");

    // Horizontal line segment which touches the arc at [0.0, 2.0]
    let ls = LineSegment::new([-2.0, 2.0], [2.0, 2.0]).expect("valid_arc");
    assert!(arc.is_tangent(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));

    // The underlying line is a tangent, but point [0.0, 2.0] is not part of the line segment
    let ls = LineSegment::new([-2.0, 2.0], [-4.0, 2.0]).expect("valid_arc");
    assert!(!arc.is_tangent(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));

    // Line segment intersects the arc twice
    let ls = LineSegment::new([-2.0, 1.9], [2.0, 1.9]).expect("valid_arc");
    assert!(!arc.is_tangent(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));

    // Vertical line segment which touches the arc at [2.0, 0.0]
    let ls = LineSegment::new([2.0, 0.0], [2.0, 1.0]).expect("valid_arc");
    assert!(arc.is_tangent(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));

    // Line segment touches the arc at [2.0, 0.0], but is not perpendicular to [0.0, 0.0] - [2.0, 0.0]
    let ls = LineSegment::new([2.0, 0.0], [3.0, 1.0]).expect("valid_arc");
    assert!(!arc.is_tangent(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ```
     */
    pub fn is_tangent(
        &self,
        line_segment: &super::LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> bool {
        return line_segment.is_tangent(self, epsilon, max_relative);
    }

    /// Returns whether `self` and `other` are touching.
    ///
    /// The two segments are touching if they are intersecting but not
    /// dividing each other. If `other` is a
    /// [`LineSegment`](super::LineSegment), then it must be a tangent of
    /// `self`. If `other` is an [`ArcSegment`], `self` and `other` must either
    /// lay on the same circle ([`ArcSegment::same_circle`] is true) or have a
    /// single intersection point and their tangents at this point must be
    /// identical.
    #[doc = ""]
    #[cfg_attr(feature = "doc-images", doc = "![Touching and dividing][touching]")]
    #[cfg_attr(
        feature = "doc-images",
        embed_doc_image::embed_doc_image("touching", "docs/img/touching.svg")
    )]
    #[cfg_attr(
        not(feature = "doc-images"),
        doc = "**Doc images not enabled**. Compile docs with `cargo doc --features 'doc-images'` and Rust version >= 1.54."
    )]
    ///
    /// # Examples
    ///
    /// ```
    /// use planar_geo::prelude::*;
    ///
    /// // Circle around [0, 0] with radius 2
    /// let c1 = ArcSegment::circle([0.0, 0.0], 2.0).unwrap();
    ///
    /// // Arc segment around [0, 0] with radius 2
    /// let a1 = ArcSegment::from_center_radius_start_sweep_angle([0.0, 0.0], 2.0, 0.0, 1.0).expect("valid_arc");
    ///
    /// // Circle around [0, 4] with radius 2
    /// let c2 = ArcSegment::circle([0.0, 4.0], 2.0).unwrap();
    ///
    /// // Circle around [0, 3] with radius 2
    /// let c3 = ArcSegment::circle([0.0, 3.0], 2.0).unwrap();
    ///
    /// // Line segment from [-2, -2] to [2, -2]
    /// let ls = LineSegment::new([-2.0, -2.0], [2.0, -2.0]).expect("valid_arc");
    ///
    /// // c1 and a1 lay on the same circle -> Touching
    /// assert!(c1.touches_segment(&a1, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ///
    /// // c1 and c2 touch in a single point
    /// assert!(c1.touches_segment(&c2, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ///
    /// // c1 and c3 overlap (two intersections)
    /// assert!(!c1.touches_segment(&c3, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    ///
    /// // ls is a tangent of c1
    /// assert!(c1.touches_segment(&ls, DEFAULT_EPSILON, DEFAULT_MAX_RELATIVE));
    /// ```
    pub fn touches_segment<'a, T: Into<super::SegmentRef<'a>>>(
        &self,
        other: T,
        epsilon: f64,
        max_relative: f64,
    ) -> bool {
        return self
            .touches_and_intersections(other.into(), epsilon, max_relative)
            .0;
    }

    pub(crate) fn touches_and_intersections(
        &self,
        other: super::SegmentRef<'_>,
        epsilon: f64,
        max_relative: f64,
    ) -> (bool, PrimitiveIntersections) {
        fn ep_as(s: &super::ArcSegment, pt: [f64; 2], epsilon: f64, max_relative: f64) -> bool {
            return relative_eq!(
                pt,
                s.start(),
                epsilon = epsilon,
                max_relative = max_relative
            ) || relative_eq!(pt, s.stop(), epsilon = epsilon, max_relative = max_relative);
        }

        match other {
            super::SegmentRef::LineSegment(line_segment) => {
                return line_segment.touches_and_intersections(self.into(), epsilon, max_relative);
            }
            super::SegmentRef::ArcSegment(arc_segment) => {
                if std::ptr::eq(self, arc_segment) {
                    return (true, PrimitiveIntersections::Zero);
                }

                let intersections =
                    self.intersections_arc_segment(arc_segment, epsilon, max_relative);
                let touches = match intersections {
                    PrimitiveIntersections::Zero => false,
                    PrimitiveIntersections::One(i) => {
                        if ep_as(self, i, epsilon, max_relative)
                            || ep_as(arc_segment, i, epsilon, max_relative)
                        {
                            true
                        } else {
                            // Arc segment are touching if their tangent is the
                            // same. This is the case if the vectors i - center are
                            // colinear.
                            let a = [i[0] - self.center[0], i[1] - self.center[1]];
                            let b = [i[0] - arc_segment.center[0], i[1] - arc_segment.center[1]];
                            let cross = a[0] * b[1] - a[1] * b[0];
                            relative_eq!(cross, 0.0, epsilon = epsilon, max_relative = max_relative)
                        }
                    }
                    PrimitiveIntersections::Two([i1, i2]) => {
                        // If center and radius are identical, the arc segments are part
                        // of the same circle and therefore touch
                        if self.same_circle(&arc_segment, epsilon, max_relative) {
                            true
                        } else {
                            // Are the intersections end points?
                            (ep_as(self, i1, epsilon, max_relative)
                                || ep_as(arc_segment, i1, epsilon, max_relative))
                                && (ep_as(self, i2, epsilon, max_relative)
                                    || ep_as(arc_segment, i2, epsilon, max_relative))
                        }
                    }
                };
                return (touches, intersections);
            }
        }
    }

    /// Intersects a (world-space) infinite line with this circle.
    ///
    /// This function uses a modified version of the circle-line intersection
    /// algorithm from cp-algorithms.com:
    /// https://cp-algorithms.com/geometry/circle-line-intersection.html
    ///
    /// The input line is assumed to be in WORLD coordinates:
    ///
    /// a x + b y + c = 0
    ///
    /// where (a, b) is NOT required to be normalized.
    ///
    /// The circle has center (ox, oy) and radius r.
    ///
    /// ------------------------------------------------------------
    /// COORDINATE MODEL
    /// ------------------------------------------------------------
    ///
    /// This function internally translates the problem into the
    /// circle-centered coordinate system:
    ///
    /// x' = x - ox
    /// y' = y - oy
    ///
    /// After translation, the circle becomes:
    ///
    /// x'^2 + y'^2 = r^2
    ///
    /// and the line becomes:
    ///
    /// a x' + b y' + c' = 0
    /// where c' = c + a·ox + b·oy
    ///
    /// IMPORTANT:
    /// - The translation is purely internal.
    /// - All returned points are converted back to WORLD space.
    /// - No caller-visible coordinate transform is required.
    ///
    /// ------------------------------------------------------------
    /// GEOMETRIC METHOD
    /// ------------------------------------------------------------
    ///
    /// 1. Project origin onto the line in circle-centered space.
    /// 2. Compute distance from circle center to line.
    /// 3. If distance > r → no intersection.
    /// 4. If distance == r → tangent (1 point).
    /// 5. If distance < r → 2 intersection points.
    ///
    /// Direction vector along the line is perpendicular to (a, b).
    pub(crate) fn intersections_line_circle(
        &self,
        a: f64,
        b: f64,
        c: f64,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        let [ox, oy] = self.center();
        let r = self.radius();

        // Translate line into circle-centered coordinate system
        let c = c + a * ox + b * oy;

        let normal_len_sq = a * a + b * b;

        // Projection of origin onto the line (in local frame)
        let x0 = -a * c / normal_len_sq;
        let y0 = -b * c / normal_len_sq;

        // squared distance from origin to line
        let dist_to_line_sq = c * c / normal_len_sq;

        // remaining squared radius after projection
        let radial_discriminant = r * r - dist_to_line_sq;

        let mut intersections = PrimitiveIntersections::Zero;

        // Tangent case (1 intersection)
        if relative_eq!(
            radial_discriminant,
            0.0,
            epsilon = epsilon,
            max_relative = max_relative
        ) {
            let pt = [x0 + ox, y0 + oy];
            if self.covers_point(pt, epsilon, max_relative) {
                intersections.push(pt);
            }
        } else if radial_discriminant > 0.0 {
            // Secant case (2 intersections)

            let m = (radial_discriminant / normal_len_sq).sqrt();

            // Points in circle-centered frame
            let xa = x0 + b * m;
            let ya = y0 - a * m;
            let xb = x0 - b * m;
            let yb = y0 + a * m;

            // Convert to world coordinates BEFORE angle test
            let pa = [xa + ox, ya + oy];
            let pb = [xb + ox, yb + oy];

            // Angles are defined relative to circle center in WORLD space
            let angle_a = (pa[1] - oy).atan2(pa[0] - ox);
            let angle_b = (pb[1] - oy).atan2(pb[0] - ox);

            if self.covers_angle(angle_a) {
                intersections.push([xa + ox, ya + oy]);
            }
            if self.covers_angle(angle_b) {
                intersections.push([xb + ox, yb + oy]);
            }
        }
        return intersections;
    }

    /// Split arc in y-monotonic segments (where each y-value appears at most
    /// once)
    pub(crate) fn split_y_monotonic(&self) -> [Option<ArcSegment>; 3] {
        // Cover special case of a circle
        let start_angle = if self.is_circle() {
            -FRAC_PI_2
        } else {
            self.start_angle()
        };

        let mut angles = [Some(start_angle), None, None, None];
        let mut free_slot = 1;
        if self.is_positive() {
            if self.start_angle() > FRAC_PI_2 {
                if self.covers_angle(3.0 * FRAC_PI_2) {
                    angles[free_slot] = Some(3.0 * FRAC_PI_2);
                    free_slot += 1;
                }
                if self.covers_angle(FRAC_PI_2) {
                    angles[free_slot] = Some(FRAC_PI_2);
                    free_slot += 1;
                }
            } else {
                if self.covers_angle(FRAC_PI_2) {
                    angles[free_slot] = Some(FRAC_PI_2);
                    free_slot += 1;
                }
                if self.covers_angle(3.0 * FRAC_PI_2) {
                    angles[free_slot] = Some(3.0 * FRAC_PI_2);
                    free_slot += 1;
                }
            }
        } else {
            if self.start_angle() < FRAC_PI_2 {
                if self.covers_angle(3.0 * FRAC_PI_2) {
                    angles[free_slot] = Some(3.0 * FRAC_PI_2);
                    free_slot += 1;
                }
                if self.covers_angle(FRAC_PI_2) {
                    angles[free_slot] = Some(FRAC_PI_2);
                    free_slot += 1;
                }
            } else {
                if self.covers_angle(FRAC_PI_2) {
                    angles[free_slot] = Some(FRAC_PI_2);
                    free_slot += 1;
                }
                if self.covers_angle(3.0 * FRAC_PI_2) {
                    angles[free_slot] = Some(3.0 * FRAC_PI_2);
                    free_slot += 1;
                }
            }
        }
        if self.start_angle().rem_euclid(TAU) != self.stop_angle().rem_euclid(TAU) {
            angles[free_slot] = Some(self.stop_angle().rem_euclid(TAU));
        }

        let mut arcs = [None, None, None];
        free_slot = 0;
        for w in angles.windows(2) {
            // Check if both angles are "Some"
            let start_angle: f64 = match w[0] {
                Some(a) => a,
                None => break,
            };
            let mut sweep_angle = match w[1] {
                Some(a) => {
                    let delta = (a - start_angle).rem_euclid(TAU);
                    if delta > PI {
                        (start_angle - a).rem_euclid(TAU)
                    } else {
                        delta
                    }
                }
                None => break,
            };
            if !self.is_positive() {
                sweep_angle = -sweep_angle;
            }

            if let Ok(arc) = ArcSegment::from_center_radius_start_sweep_angle(
                self.center(),
                self.radius(),
                start_angle,
                sweep_angle,
            ) {
                arcs[free_slot] = Some(arc);
                free_slot += 1;
            }
        }

        return arcs;
    }
}

impl crate::primitive::private::Sealed for ArcSegment {}

impl Primitive for ArcSegment {
    fn covers_point(&self, point: [f64; 2], epsilon: f64, max_relative: f64) -> bool {
        if relative_eq!(
            self.start(),
            point,
            epsilon = epsilon,
            max_relative = max_relative
        ) {
            return true;
        }
        if relative_eq!(
            self.stop(),
            point,
            epsilon = epsilon,
            max_relative = max_relative
        ) {
            return true;
        }

        // Quick first check: If point p is outside the segment bounding box, it can't
        // be part of the arc segment.
        if !BoundingBox::from(self).approx_covers_point(point, epsilon) {
            return false;
        }

        let shifted_pt = [point[0] - self.center()[0], point[1] - self.center()[1]];

        if relative_eq!(
            self.radius().powi(2),
            shifted_pt[0].powi(2) + shifted_pt[1].powi(2),
            epsilon = epsilon,
            max_relative = max_relative
        ) {
            let angle = shifted_pt[1].atan2(shifted_pt[0]);
            return self.covers_angle(angle);
        } else {
            return false;
        }
    }

    fn covers_arc_segment(&self, arc_segment: &Self, epsilon: f64, max_relative: f64) -> bool {
        if std::ptr::eq(self, arc_segment) {
            return true;
        }

        // Special treatment of a full circle: other is contained if center and
        // radius are identical to that of self.
        if self.is_circle() {
            return relative_eq!(
                self.center(),
                arc_segment.center(),
                epsilon = epsilon,
                max_relative = max_relative
            ) && relative_eq!(
                self.radius(),
                arc_segment.radius(),
                epsilon = epsilon,
                max_relative = max_relative
            );
        }

        match self.intersections_primitive(arc_segment, epsilon, max_relative) {
            // Deal with special case where self and other are identical
            PrimitiveIntersections::Zero => false,
            PrimitiveIntersections::One(_) => false,
            PrimitiveIntersections::Two([pt1, pt2]) => {
                let start = arc_segment.start();
                let stop = arc_segment.stop();

                if relative_ne!(start, pt1, epsilon = epsilon, max_relative = max_relative)
                    && relative_ne!(stop, pt1, epsilon = epsilon, max_relative = max_relative)
                {
                    false
                } else if relative_ne!(start, pt2, epsilon = epsilon, max_relative = max_relative)
                    && relative_ne!(stop, pt2, epsilon = epsilon, max_relative = max_relative)
                {
                    false
                } else {
                    // Middle angle of other must be on self
                    self.covers_angle(arc_segment.start_angle() + 0.5 * arc_segment.sweep_angle())
                }
            }
        }
    }

    fn covers_line_segment(
        &self,
        _line_segment: &crate::segment::LineSegment,
        _epsilon: f64,
        _max_relative: f64,
    ) -> bool {
        return false;
    }

    fn covers_line(&self, _line: &crate::line::Line, _epsilon: f64, _max_relative: f64) -> bool {
        return false;
    }

    fn intersections_line(
        &self,
        line: &crate::line::Line,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        line.intersections_arc_segment(self, epsilon, max_relative)
    }

    fn intersections_line_segment(
        &self,
        line_segment: &crate::segment::LineSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        line_segment.intersections_arc_segment(self, epsilon, max_relative)
    }

    /// This function uses a modified version of the circle-circle intersection
    /// algorithm from cp-algorithms.com:
    /// https://cp-algorithms.com/geometry/circle-circle-intersection.html
    fn intersections_arc_segment(
        &self,
        arc_segment: &ArcSegment,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        // If the segments are identical, they do not intersect by definition
        if std::ptr::eq(self, arc_segment) {
            return PrimitiveIntersections::Zero;
        }

        /*
        Compute intersection of two circles via radical axis reduction.

        We derive a line (radical axis) in WORLD coordinates:

            a x + b y + c = 0

        where:
            a = 2 (C1x - C2x)
            b = 2 (C1y - C2y)
            c = C2·C2 - C1·C1 + r1² - r2²

        This line represents all points with equal power w.r.t. both circles.

        IMPORTANT INVARIANT:
        - No coordinate translation is performed here.
        - This line is already in WORLD space.
        - Any local coordinate handling is done inside `intersections_line_circle`.
        */
        let radius_arc1 = self.radius();
        let radius_arc2 = arc_segment.radius();
        let center_arc1 = self.center();
        let center_arc2 = arc_segment.center();

        // Circles are identical
        if relative_eq!(
            radius_arc1,
            radius_arc2,
            epsilon = epsilon,
            max_relative = max_relative
        ) && relative_eq!(
            center_arc1,
            center_arc2,
            epsilon = epsilon,
            max_relative = max_relative
        ) {
            let mut intersections = PrimitiveIntersections::Zero;
            if self.covers_angle(arc_segment.start_angle()) {
                intersections.push(arc_segment.start());
            }
            if self.covers_angle(arc_segment.stop_angle()) {
                intersections.push(arc_segment.stop());
            }
            if intersections.len() == 2 {
                return intersections;
            }
            if arc_segment.covers_angle(self.start_angle()) {
                intersections.push(self.start());
            }
            if intersections.len() == 2 {
                return intersections;
            }
            if arc_segment.covers_angle(self.stop_angle()) {
                intersections.push(self.stop());
            }
            return intersections;
        }

        // Radical axis (WORLD coordinates)
        let a = 2.0 * (center_arc1[0] - center_arc2[0]);
        let b = 2.0 * (center_arc1[1] - center_arc2[1]);
        let c = center_arc2[0].powi(2) + center_arc2[1].powi(2)
            - center_arc1[0].powi(2)
            - center_arc1[1].powi(2)
            + radius_arc1.powi(2)
            - radius_arc2.powi(2);

        /*
        General case:

        1. Intersect radical axis (line) with THIS circle
        2. This returns points in WORLD coordinates
        3. Filter by whether points lie inside arc_segment
        */
        let mut intersections = PrimitiveIntersections::Zero;
        for i in self.intersections_line_circle(a, b, c, epsilon, max_relative) {
            if arc_segment.covers_point(i, epsilon, max_relative) {
                intersections.push(i);
            }
        }
        return intersections;
    }

    fn intersections_primitive<T: Primitive>(
        &self,
        other: &T,
        epsilon: f64,
        max_relative: f64,
    ) -> PrimitiveIntersections {
        other.intersections_arc_segment(self, epsilon, max_relative)
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

        // sweep angle is simply inverted
        self.sweep_angle = -self.sweep_angle;

        // Calculate the angle between start and stop
        let line_angle = (stop[1] - start[1]).atan2(stop[0] - start[0]);
        let diff = line_angle - self.start_angle;
        self.start_angle = (self.start_angle + 2.0 * diff).rem_euclid(TAU);
    }
}

impl ToBoundingBox for ArcSegment {
    fn bounding_box(&self) -> BoundingBox {
        let c = self.center();

        let mut xmin_pt = c.clone();
        xmin_pt.translate([-self.radius(), 0.0]);

        let mut xmax_pt = c.clone();
        xmax_pt.translate([self.radius(), 0.0]);

        let mut ymin_pt = c.clone();
        ymin_pt.translate([0.0, -self.radius()]);

        let mut ymax_pt = c.clone();
        ymax_pt.translate([0.0, self.radius()]);

        let start = self.start();
        let stop = self.stop();

        let xmin = if self.covers_angle(PI) {
            xmin_pt[0]
        } else {
            // Smaller of the two values
            if start[0] < stop[0] {
                start[0]
            } else {
                stop[0]
            }
        };

        let xmax = if self.covers_angle(0.0) {
            xmax_pt[0]
        } else {
            // Larger of the two values
            if start[0] > stop[0] {
                start[0]
            } else {
                stop[0]
            }
        };

        let ymin = if self.covers_angle(-FRAC_PI_2) {
            ymin_pt[1]
        } else {
            // Smaller of the two values
            if start[1] < stop[1] {
                start[1]
            } else {
                stop[1]
            }
        };

        let ymax = if self.covers_angle(FRAC_PI_2) {
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

impl approx::AbsDiffEq for ArcSegment {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        return self.start().abs_diff_eq(&other.start(), epsilon)
            && self
                .segment_point(0.5)
                .abs_diff_eq(&other.segment_point(0.5), epsilon)
            && self.stop().abs_diff_eq(&other.stop(), epsilon);
    }
}

impl approx::RelativeEq for ArcSegment {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        return self
            .start()
            .relative_eq(&other.start(), epsilon, max_relative)
            && self.segment_point(0.5).relative_eq(
                &other.segment_point(0.5),
                epsilon,
                max_relative,
            )
            && self
                .stop()
                .relative_eq(&other.stop(), epsilon, max_relative);
    }
}

impl approx::UlpsEq for ArcSegment {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: f64, max_ulps: u32) -> bool {
        return self.start().ulps_eq(&other.start(), epsilon, max_ulps)
            && self
                .segment_point(0.5)
                .ulps_eq(&other.segment_point(0.5), epsilon, max_ulps)
            && self.stop().ulps_eq(&other.stop(), epsilon, max_ulps);
    }
}

impl From<&ArcSegment> for CentroidData {
    fn from(value: &ArcSegment) -> Self {
        let inv_3 = 1.0 / 3.0;

        let center = value.center();
        let angle = value.sweep_angle();
        let radius = value.radius();

        // centroid and area for shape 1
        let x = (value.start()[0] + value.stop()[0]) * inv_3;
        let y = (value.start()[1] + value.stop()[1]) * inv_3;
        let area = 0.5 * (value.start()[0] * value.stop()[1] - value.stop()[0] * value.start()[1]);
        let data_shape_1 = CentroidData { area, x, y };

        // centroid and area for shape 2
        // Centroid of a circular segment is calculated by moving the segment so that:
        // a) Center equals origin
        // b) circular segment is symmetric around the y-axis
        // Calculate the centroid and reverse the operation.
        let half_angle = 0.5 * angle;
        let area = half_angle * radius.powi(2);
        let centroid_transformed = [2.0 / 3.0 * radius * half_angle.sin() / half_angle, 0.0];

        let middle_angle = value.start_angle() + 0.5 * value.sweep_angle();
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
Calculates the sweep angle from a `start_angle` and a `stop_angle`. The arc
between a start and a stop angle can either be mathematically `positive`
(i.e. counter-clockwise) or negative (not `positive`, i.e. clockwise).
 */
pub fn calculate_sweep_angle(start_angle: f64, stop_angle: f64, positive: bool) -> f64 {
    let sweep_angle: f64;
    if positive {
        if start_angle < stop_angle {
            sweep_angle = stop_angle - start_angle;
        } else {
            sweep_angle = TAU - (start_angle - stop_angle);
        }
    } else {
        if start_angle > stop_angle {
            sweep_angle = stop_angle - start_angle;
        } else {
            sweep_angle = (stop_angle - start_angle) - TAU;
        }
    }
    return sweep_angle;
}
