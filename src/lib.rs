#![cfg_attr(docsrs, doc = include_str!("../README.md"))]
#![cfg_attr(not(docsrs), doc = include_str!("../README_local.md"))]
#![deny(missing_docs)]

/**
A resonable value for the absolute tolerance.

Comparing floating-point numbers is often necessary within this crate, e.g.
when deciding whether two geometries intersect or not. Since floating point
numbers do not have arbitrary precision, most real numbers cannot be represented
exactly (e.g. 0.1). Hence, comparisons within this crate are done using the
[`ulps_eq`](approx::UlpsEq::ulps_eq) function from the [`approxim`](approx)
crate, which requires specifying an absolute tolerance and a maximum units in
last place (ULPs) deviation. The former is mainly relevant when comparing
numbers close to zero, while the latter gets important with (absolute) big
numbers.

The following links taken directly from the [`approxim`](approx) crate contain
more information regarding the behaviour of floating point numbers, particulary
when comparing them:
- [Comparing Floating Point Numbers, 2012 Edition](https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/)
- [The Floating Point Guide - Comparison](https://floating-point-gui.de/errors/comparison/)
- [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)

This constant and [`DEFAULT_MAX_ULPS`] provide sane defaults which will work
well in most cases and are recommeded to be used in the various functions
requiring an `epsilon` and a [`max_ulps`] argument. The reason they need to be
provided as explicit arguments and are not simply used by default is that expert
users can modify the comparison behaviour of e.g. the intersection algorithms
if required for a particular use case.

The specified value is the squere root of the machine precision
(`f64::EPSILON.sqrt()`).
*/
pub const DEFAULT_EPSILON: f64 = 0.000000014901161193847656_f64; //

/**
A resonable default value for the maximum units in last place (ULPs) tolerance.
The (maximum) ULPs tolerance is used in a lot of functions across this crate,
see the docstring of [`DEFAULT_EPSILON`] for more.
 */
pub const DEFAULT_MAX_ULPS: u32 = 4;

// Base
pub mod error;
pub mod prelude;

// Primitives
pub mod line;
pub mod primitive;
pub mod segment;

// Composite types
pub mod composite;
pub mod contour;
pub mod segment_chain;
pub mod shape;

#[cfg(feature = "visualize")]
pub mod visualize;

// =============================================================================
// Shared code

/**
Affine transformations for geometric types.

All geometric types within this crate as well as the basic "point" type
`[f64;2`] implement this trait to allow for easy affine transformations. All
examples in the docstrings of the individual trait methods are for the point
type, because all other geometric types are based on it. Hence, their trait
implementation just delegates to `impl Transformation for [f64; 2]` for all
their points.
 */
pub trait Transformation {
    /**
    Translates `self` by the given `shift`.

    ```
    use planar_geo::Transformation;

    let mut pt = [1.0, 1.0];
    pt.translate([1.0, -1.0]);
    assert_eq!(pt, [2.0, 0.0]);
    ```
     */
    fn translate(&mut self, shift: [f64; 2]);

    /**
    Rotates `self` around the `center` by the given `angle` (in rad).

    ```
    use std::f64::consts::PI;
    use planar_geo::Transformation;

    let mut pt = [1.0, 1.0];
    pt.rotate([1.0, -1.0], PI);
    approx::assert_abs_diff_eq!(pt, [1.0, -3.0], epsilon = 1e-15);
    ```
     */
    fn rotate(&mut self, center: [f64; 2], angle: f64);

    /**
    Scales `self` by `factor` with respect to the origin `[0.0, 0.0]`.

    ```
    use planar_geo::Transformation;

    let mut pt = [1.0, 1.0];
    pt.scale(2.0);
    assert_eq!(pt, [2.0, 2.0]);
    ```
     */
    fn scale(&mut self, factor: f64);

    /**
    Mirrors `self` about a line defined by two points.

    ```
    use planar_geo::Transformation;

    let mut pt = [1.0, 1.0];
    pt.line_reflection([0.0, 0.0], [0.0, 2.0]); // Vertical line
    approx::assert_abs_diff_eq!(pt, [-1.0, 1.0], epsilon = 1e-15);
    pt.line_reflection([0.0, -1.0], [2.0, -1.0]); // Horizontal line
    approx::assert_abs_diff_eq!(pt, [-1.0, -3.0], epsilon = 1e-15);
    ```
     */
    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> ();

    /**
    Mirrors `self` about a `point`. This operation is equivalent to a rotation
    around the point with the angle PI.

    ```
    use planar_geo::prelude::*;

    let mut pt = [1.0, 1.0];
    pt.point_reflection([0.0, 0.0]);
    approx::assert_abs_diff_eq!(pt, [-1.0, -1.0], epsilon = 1e-15);
    pt.point_reflection([0.0, 2.0]);
    approx::assert_abs_diff_eq!(pt, [1.0, 5.0], epsilon = 1e-15);
    ```
     */
    fn point_reflection(&mut self, point: [f64; 2]) -> () {
        self.rotate(point, std::f64::consts::PI);
    }
}

impl Transformation for [f64; 2] {
    fn translate(&mut self, shift: [f64; 2]) {
        *self = [self[0] + shift[0], self[1] + shift[1]];
    }

    fn rotate(&mut self, center: [f64; 2], angle: f64) {
        let t = crate::Rotation2::new(angle);
        let pt = [self[0] - center[0], self[1] - center[1]];
        let mut p = t * pt;
        p.translate([center[0], center[1]]);
        *self = p;
    }

    fn scale(&mut self, factor: f64) {
        *self = [self[0] * factor, self[1] * factor];
    }

    fn line_reflection(&mut self, start: [f64; 2], stop: [f64; 2]) -> () {
        // Treatment of special case: vertical line
        if start[0] == stop[0] {
            *self = [-self[0] + 2.0 * start[0], self[1]];

        // Treatment of special case: Horizontal line
        } else if start[1] == stop[1] {
            *self = [self[0], -self[1] + 2.0 * start[1]];
        // All other cases
        } else {
            // Solve the line equation
            let m = (stop[1] - start[1]) / (stop[0] - start[0]);
            let c = (stop[0] * start[1] - start[0] * stop[1]) / (stop[0] - start[0]);

            let d = (self[0] + (self[1] - c) * m) / (1.0 + m.powi(2));
            *self = [2.0 * d - self[0], 2.0 * d * m - self[1] + 2.0 * c];
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CentroidData {
    pub(crate) area: f64,
    pub(crate) x: f64,
    pub(crate) y: f64,
}

impl CentroidData {
    /**
    Calculates the common centroid of self and other.
     */
    pub(crate) fn union(&self, other: &Self) -> Self {
        let area = self.area + other.area;
        if area == 0.0 {
            let x = 0.0;
            let y = 0.0;
            return Self { area, x, y };
        } else {
            let x = (self.area * self.x + other.area * other.x) / area;
            let y = (self.area * self.y + other.area * other.y) / area;
            return Self { area, x, y };
        }
    }

    /**
    Calculates the common centroid of self and other, but subtract other from
    self (i.e. the original contour of other is subtracted from self)
     */
    pub(crate) fn subtract(&self, other: &Self) -> Self {
        let area = self.area - other.area;
        if area == 0.0 {
            let x = 0.0;
            let y = 0.0;
            return Self { area, x, y };
        } else {
            let x = (self.area * self.x - other.area * other.x) / area;
            let y = (self.area * self.y - other.area * other.y) / area;
            return Self { area, x, y };
        }
    }
}

impl From<CentroidData> for [f64; 2] {
    fn from(value: CentroidData) -> Self {
        return [value.x, value.y];
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Rotation2 {
    sin: f64,
    cos: f64,
}

impl Rotation2 {
    pub(crate) fn new(angle: f64) -> Self {
        return Self {
            sin: angle.sin(),
            cos: angle.cos(),
        };
    }
}

impl std::ops::Mul<[f64; 2]> for Rotation2 {
    type Output = [f64; 2];

    fn mul(self, rhs: [f64; 2]) -> Self::Output {
        return [
            rhs[0] * self.cos - rhs[1] * self.sin,
            rhs[0] * self.sin + rhs[1] * self.cos,
        ];
    }
}
