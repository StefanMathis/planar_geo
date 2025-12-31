/*!
This module reexports commonly used types, traits and functions from this crate
for ease of use. For example, a user of this crate can simply write:

```
use planar_geo::prelude::{ArcSegment, Composite, Primitive};
```

instead of

```
use planar_geo::segment::ArcSegment;
use planar_geo::composite::Composite;
use planar_geo::primitive::Primitive;
```
 */

pub use crate::Transformation;
pub use crate::composite::*;
pub use crate::contour::{ArrowHeadSize, Contour};
pub use crate::line::*;
pub use crate::primitive::{Primitive, PrimitiveIntersections};
pub use crate::segment::arc_segment::*;
pub use crate::segment::line_segment::*;
pub use crate::segment::*;
pub use crate::segment_chain::*;
pub use crate::shape::Shape;
pub use crate::{DEFAULT_EPSILON, DEFAULT_MAX_ULPS};

#[cfg(feature = "visualize")]
pub use crate::visualize::*;
pub use approx;

///! Reexport of useful common functionality
pub use bounding_box::*;
