fn main() {
    // If building for docs.rs, DO NOT create the README files from the template
    if let Ok(env) = std::env::var("DOCS_RS") {
        if &env == "1" {
            return ();
        }
    }

    let mut readme = std::fs::read_to_string("README.template.md").unwrap();
    readme = readme.replace(
        "{{VERSION}}",
        std::env::var("CARGO_PKG_VERSION")
            .expect("version is available in build.rs")
            .as_str(),
    );

    // Generate README_local.md using local images
    let mut local = readme.replace("{{type_overview.svg}}", "docs/type_overview.svg");
    local = local.replace("{{shape.svg}}", "docs/shape.svg");
    local = local.replace(
        "{{intersection_segments.svg}}",
        "docs/intersection_segments.svg",
    );
    local = local.replace(
        "{{intersection_composites.svg}}",
        "docs/intersection_composites.svg",
    );
    std::fs::write("README_local.md", local).unwrap();

    // Generate README,md using online hosted images
    let mut docsrs = readme.replace(
        "{{type_overview.svg}}",
        "https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/type_overview.svg",
    );
    docsrs = docsrs.replace(
        "{{shape.svg}}",
        "https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/shape.svg",
    );
    docsrs = docsrs.replace(
        "{{intersection_segments.svg}}",
        "https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/intersection_segments.svg",
    );
    docsrs = docsrs.replace(
        "{{intersection_composites.svg}}",
        "https://raw.githubusercontent.com/StefanMathis/planar_geo/refs/heads/main/docs/intersection_composites.svg",
    );
    std::fs::write("README.md", docsrs).unwrap();
}
