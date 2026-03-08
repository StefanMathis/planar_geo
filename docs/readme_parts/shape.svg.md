_This image was created with examples/type_overview.rs_

The "segment" types are so-called [`Primitive`]s: Simple straight
([`LineSegment`]) or arc ([`ArcSegment`]) connections between two points. They
form the basis for the [`Composite`] types [`SegmentChain`], [`Contour`] (a
closed segment chain) and [`Shape`] (composed of one or more contours).

For these types, this crate offers the following features:
- Property calculation (e.g. length, surface area, centroids),
- [`Transformation`] (scaling, shifting, rotation, mirroring),
- Intersection calculation between all combinations of the aforementioned types,
see the [`Primitive`] and [`Composite`] traits.

If the corresponding features are activated, it is also possible to serialize
and deserialize (using [serde]) and to visualize (using [gtk-rs]) these types.
See the [Features](#features) section for more.

One distinct feature of this library is the treatment of arcs: Arcs are not
approximated as polylines (i.e. a connected series of line segments), but
instead being represented as "true" arcs (see e.g. the examples for surface
area calculation below). In fact, this is the reason why this library was
written in the first place.

Since the "point" type is defined using the floating-point type `f64`, a lot of
operations (i.e. intersection calculation) are prone to rounding-errors. These
operations therefore require specifying an absolute tolerance `epsilon` and a
maximum units in last place tolerance `max_ulps`, which are used as inputs
for [`ulps_eq`] (from the [approxim] crate) to e.g. determine whether two points
are approximately equal. It is recommended to use the "default" tolerances
[`DEFAULT_EPSILON`] and [`DEFAULT_MAX_ULPS`] unless there is a good reason to
use other values.

The following paragraphs will provide some examples for the aforementioned
features.

## Construction and property calculation

The following code snippet shows how to construct the shape shown in the image
below and calculate some of its properties, e.g. centroid and surface area. The
image itself has been created by running `examples/shape.rs`.