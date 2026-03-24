These intersections have been calculated using the generic [`intersections`]
interface, which is available for any of the geometric types defined in this
crate:

```rust
use planar_geo::prelude::*;

let e = DEFAULT_EPSILON;
let m = DEFAULT_MAX_ULPS;

let vertices = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let contour = Contour::new(Polysegment::from_points(vertices));
let vertices = &[[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
let hole = Contour::new(Polysegment::from_points(vertices));
let shape = Shape::new(vec![contour, hole]).expect("valid inputs");

let vertices = &[[2.0, 1.0], [2.0, 0.5], [0.0, 0.5]];
let polysegment = Polysegment::from_points(vertices);

let intersections: Vec<Intersection> = polysegment.intersections(&shape, e, m);
assert_eq!(intersections.len(), 4);
```

# Features

## Serialization and deserialization

When the `serde` feature is enabled, all geometric objects can be serialized and
deserialized using the [serde](https://crates.io/crates/serde) crate.

## Visualization

When the `cairo` feature is enabled, all geometric objects have a `draw`
method which can be used to draw them onto a 
[cairo](https://gtk-rs.org/gtk-rs-core/stable/latest/docs/cairo/) [`Context`].
See the [module-level documentation](draw)
for more. All images used in the documentation were created using this
functionality.

## Embedded images in documentation

When building the documentation with `doc-images` enabled integrates images
into the docstrings of items using the 
[embed-doc-image](https://crates.io/crates/embed-doc-image) crate.

# Documentation

When building the documentation, it is recommended to enable all features with
`cargo doc --all-features`; otherwise the generated documentation will not have
any images and will miss the items hidden behind feature flags.