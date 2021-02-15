use nalgebra::{Point3, Vector3};

mod graphics_utils;
mod renderer;

// Image
const ASPECT_RATIO: f64 = 16. / 9.;
const IMAGE_WIDTH: usize = 400;
const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;

fn main() {
    let viewport_height = 2.0;
    let viewport_width = ASPECT_RATIO * viewport_height;
    let focal_length = 1.0;
    let mut scene = renderer::RenderedScene::new(
        ASPECT_RATIO,
        IMAGE_HEIGHT,
        IMAGE_WIDTH,
        viewport_height,
        viewport_width,
        focal_length,
        graphics_utils::ColorRGB::new(0., 0., 0.),
    );

    println!("Image height: {}", IMAGE_HEIGHT);

    // Render
    let plane = renderer::Plane {
        top_left_corner: Point3::new(0., 0., 0.),
        bottom_right_corner: Point3::new(0., 0., 0.),
    };
    scene.render_object(&plane);

    let world_sphere = renderer::Sphere {
        center: Point3::new(0., -100.5, -1.),
        radius: 100.,
    };
    scene.render_object(&world_sphere);

    // let sphere = renderer::Sphere {
    //     center: Point3::new(0., 0., -1.),
    //     radius: 0.5,
    // };
    // scene.render_object(&sphere);

    // Export
    let filename = "./test.png";
    match scene.export_to_png(filename) {
        Ok(()) => println!("File exported to {}", filename),
        Err(e) => println!("{}", e),
    }
}
