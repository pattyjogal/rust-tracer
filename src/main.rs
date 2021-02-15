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

    // Render
    for j in (0..IMAGE_HEIGHT).rev() {
        for i in 0..IMAGE_WIDTH {
            scene.add_pixel(i, j);
        }
    }

    // Export
    let filename = "./test.png";
    match scene.export_to_png(filename) {
        Ok(()) => println!("File exported to {}", filename),
        Err(e) => println!("{}", e),
    }
}
