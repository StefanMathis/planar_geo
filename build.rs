use std::fs;

fn main() {
    // Skip README generation on docs.rs
    if std::env::var("DOCS_RS").as_deref() == Ok("1") {
        return;
    }

    /*
    Compose README.md from docs/links.md, docs/main.md and (if available) end.md
    All {{VERSION}} placeholders in links.md are replaced by the actual version read from Cargo.toml.
     */
    let mut readme = fs::read_to_string("docs/links.md")
        .expect("Failed to read links.md")
        .replace(
            "{{VERSION}}",
            &std::env::var("CARGO_PKG_VERSION")
                .expect("version is available when running build.rs"),
        );
    readme.push_str(&fs::read_to_string("docs/main.md").expect("Failed to read main.md"));
    if let Ok(end) = fs::read_to_string("docs/end.md") {
        readme.push_str(&end);
    }

    let _ = fs::write("README.md", readme);
}
