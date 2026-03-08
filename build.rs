use std::{fs, path::Path};

struct Image {
    placeholder: String,
    filename: String,
    remote_url: String,
}

impl Image {
    fn new(filename: &str) -> Self {
        let placeholder = format!("{{{{{}}}}}", filename); // e.g., "{{relative_permeability.svg}}"
        let remote_url = format!(
            "https://raw.githubusercontent.com/StefanMathis/stem_material/refs/heads/main/images/{}",
            filename
        );

        Self {
            placeholder,
            filename: filename.to_string(),
            remote_url,
        }
    }
}

fn main() {
    // Skip README generation on docs.rs
    if std::env::var("DOCS_RS").as_deref() == Ok("1") {
        return;
    }

    let readme_template =
        fs::read_to_string("README.template.md").expect("Failed to read README template");

    let versioned_readme = readme_template.replace(
        "{{VERSION}}",
        &std::env::var("CARGO_PKG_VERSION").expect("version is available in build.rs"),
    );

    let images = [
        Image::new("type_overview.svg"),
        Image::new("shape.svg"),
        Image::new("intersection_segments.svg"),
        Image::new("intersection_composites.svg"),
    ];

    let out_dir = std::env::var("CARGO_MANIFEST_DIR").expect("manifest is available in build.rs");
    let doc_dir = format!("{}/target/doc/stem_material", out_dir);
    let _ = fs::create_dir_all(&doc_dir);

    // Generate README_local.md
    let mut local_readme = versioned_readme.clone();
    for img in &images {
        let local_path = format!("docs/{}", img.filename);
        local_readme = local_readme.replace(&img.placeholder, &local_path);

        let img_src = Path::new(&out_dir).join(&local_path);
        let img_dst = Path::new(&doc_dir).join(&local_path);
        let _ = fs::copy(img_src, img_dst);
    }
    let _ = fs::write("README_local.md", local_readme);

    // Generate README.md for docs.rs with remote images
    let mut docsrs_readme = versioned_readme.clone();
    for img in &images {
        docsrs_readme = docsrs_readme.replace(&img.placeholder, &img.remote_url);
    }
    let _ = fs::write("README.md", docsrs_readme);
}
