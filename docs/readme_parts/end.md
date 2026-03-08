# Features

## Serialization and deserialization

When the `serde` feature is enabled, all geometric objects can be serialized and
deserialized using the [serde](https://crates.io/crates/serde) crate.

## Visualization

When the `cairo` feature is enabled, all geometric objects have a `draw`
method which can be used to draw them onto a 
[cairo](https://gtk-rs.org/gtk-rs-core/stable/latest/docs/cairo/) [`Context`].
See the [module-level documentation](visualize)
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