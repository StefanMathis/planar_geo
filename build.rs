use std::fs;

fn main() {
    // Skip README generation on docs.rs
    if std::env::var("DOCS_RS").as_deref() == Ok("1") {
        return;
    }

    /*
    Compose README.md from the building blocks in docs/readme_parts,
    interleaving image links in between. Finally, all {{VERSION}} placeholders
    are replaced by the actual version read from Cargo.toml.
     */

    let mut readme =
        fs::read_to_string("docs/readme_parts/links.md").expect("Failed to read template");
    readme.push_str(
        &fs::read_to_string("docs/readme_parts/type_overview.svg.md")
            .expect("Failed to read template"),
    );
    readme.push_str("\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/images/type_overview.svg \"Overview geometric types of the planar_geo crate\")\n\n");

    readme.push_str(
        &fs::read_to_string("docs/readme_parts/shape.svg.md").expect("Failed to read template"),
    );
    readme.push_str("\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/images/shape.svg \"Example shape\")\n\n");

    readme.push_str(
        &fs::read_to_string("docs/readme_parts/intersection_segments.svg.md")
            .expect("Failed to read template"),
    );
    readme.push_str("\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/images/intersection_segments.svg \"Point and segment intersection\")\n\n");

    readme.push_str(
        &fs::read_to_string("docs/readme_parts/intersection_composites.svg.md")
            .expect("Failed to read template"),
    );
    readme.push_str("\n\n![](https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/images/intersection_composites.svg \"Intersection between contours and a segment chain\")\n\n");

    readme.push_str(
        &fs::read_to_string("docs/readme_parts/end.md").expect("Failed to read template"),
    );

    let readme = readme.replace(
        "{{VERSION}}",
        &std::env::var("CARGO_PKG_VERSION").expect("version is available in build.rs"),
    );
    let _ = fs::write("README.md", readme);
}
