use nalgebra::{Point3, Vector3};

mod camera;
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
    println!("Image height: {}", IMAGE_HEIGHT);

    // Render
    let camera = camera::Camera::new(
        Point3::new(0., 0., 0.),
        Point3::new(0., 0., -1.),
        Vector3::new(0., 1., 0.),
        90.,
        ASPECT_RATIO,
        camera::CameraMode::Perspective,
    );

    let red_mat = renderer::Material::new(0.3, 0.7, graphics_utils::ColorRGB::new(1., 0., 0.));
    let teal_mat = renderer::Material::new(0.3, 0.7, graphics_utils::ColorRGB::new(0., 0.5, 0.5));

    let light = renderer::PointLight {
        point: Point3::new(1., 1., 1.),
        color: graphics_utils::ColorRGB::new(1., 1., 1.),
    };

    let objects: Vec<Box<dyn renderer::Renderable>> = vec![
        // Box::new(renderer::Sphere {
        //     center: Point3::new(0., -100.5, -1.),
        //     radius: 100.,
        // }),
        Box::new(renderer::Plane {
            point: Point3::new(0., 0.5, 0.),
            normal: Vector3::new(0., 1., 0.),
            material: teal_mat,
        }),
        Box::new(renderer::Sphere {
            center: Point3::new(-1., 0., -1.),
            radius: 0.5,
            material: red_mat,
        }),
        Box::new(renderer::Sphere {
            center: Point3::new(0., 0., -2.),
            radius: 0.5,
            material: red_mat,
        }),
        // Box::new(renderer::Sphere {s
        //     center: Point3::new(1., 0., -1.),
        //     radius: 0.5,
        // }),
    ];

    let mut scene = renderer::RenderedScene::new(
        ASPECT_RATIO,
        IMAGE_HEIGHT,
        IMAGE_WIDTH,
        viewport_height,
        viewport_width,
        focal_length,
        graphics_utils::ColorRGB::new(0., 0., 0.),
        objects,
        camera,
        1,
        light,
    );
    scene.render();

    // Export
    let filename = "./test.png";
    match scene.export_to_png(filename) {
        Ok(()) => println!("File exported to {}", filename),
        Err(e) => println!("{}", e),
    }
}
