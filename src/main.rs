use nalgebra::{Point3, Vector3};

mod camera;
mod graphics_utils;
mod renderer;

// Image
const ASPECT_RATIO: f64 = 16. / 9.;
const IMAGE_WIDTH: usize = 400;
const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;

fn main() {
    // Render
    let camera = camera::Camera::new(
        Point3::new(0., 0.25, 0.),
        Point3::new(0., 0., 100000.),
        Vector3::new(0., -0.5, -0.5),
        45.,
        ASPECT_RATIO,
        camera::CameraMode::Perspective,
    );

    let red_mat = renderer::Material::new(
        0.03,
        0.97,
        graphics_utils::ColorRGB::new(0.99215686274, 0.65490196078, 0.99215686274),
    );
    let gray_mat = renderer::Material::new(
        0.03,
        0.97,
        graphics_utils::ColorRGB::new(131. / 255., 149. / 255., 167. / 255.),
    );
    let off_gray_mat =
        renderer::Material::new(0.03, 0.97, graphics_utils::ColorRGB::new(0.25, 0.25, 0.35));

    let light = renderer::PointLight {
        point: Point3::new(20.0, 20.3, -15.5),
        color: graphics_utils::ColorRGB::new(0.6, 0.6, 0.6),
    };

    let objects: Vec<Box<dyn renderer::Renderable>> = vec![
        Box::new(renderer::Plane {
            point: Point3::new(0., 0., 0.),
            normal: Vector3::new(0., 1., 0.),
            material: gray_mat,
        }),
        // Large pink pyramid
        Box::new(renderer::Triangle {
            p0: Point3::new(-0.1, 0.0, 2.1),
            p1: Point3::new(0.2, 0.0, 2.5),
            p2: Point3::new(-0.07, 2.0, 2.15),
            material: red_mat,
        }),
        Box::new(renderer::Triangle {
            p0: Point3::new(-0.1, 0.0, 2.1),
            p1: Point3::new(-0.3, 0.0, 2.5),
            p2: Point3::new(-0.07, 2.0, 2.15),
            material: red_mat,
        }),
        Box::new(renderer::Triangle {
            p1: Point3::new(0.2, 0.0, 2.5),
            p0: Point3::new(-0.3, 0.0, 2.5),
            p2: Point3::new(-0.07, 2.0, 2.15),
            material: red_mat,
        }),
        // Small red sphere
        Box::new(renderer::Sphere {
            center: Point3::new(0.65, 0.25, 2.0),
            radius: 0.25,
            material: renderer::Material::new(
                0.07,
                0.93,
                graphics_utils::ColorRGB::new(255. / 255., 99. / 255., 72. / 255.),
            ),
        }),
        // Large orange sphere
        Box::new(renderer::Sphere {
            center: Point3::new(-1.3, 0.25, 3.7),
            radius: 0.5,
            material: renderer::Material::new(
                0.01,
                0.99,
                graphics_utils::ColorRGB::new(255. / 255., 195. / 255., 18. / 255.),
            ),
        }),
        // Small navy blue pyramid
        Box::new(renderer::Triangle {
            p0: Point3::new(-0.3, 0.0, 1.4),
            p1: Point3::new(-0.4, 0.0, 1.5),
            p2: Point3::new(-0.35, 0.25, 1.45),
            material: off_gray_mat,
        }),
        Box::new(renderer::Triangle {
            p0: Point3::new(-0.3, 0.0, 1.4),
            p1: Point3::new(-0.165, 0.0, 1.44),
            p2: Point3::new(-0.35, 0.25, 1.45),
            material: off_gray_mat,
        }),
    ];

    let mut scene = renderer::RenderedScene::new(
        IMAGE_HEIGHT,
        IMAGE_WIDTH,
        graphics_utils::ColorRGB::new(0., 0., 0.),
        objects,
        camera,
        8,
        light,
    );
    scene.render();

    // Export
    let filename = format!("./test.png");
    match scene.export_to_png(&filename) {
        Ok(()) => println!("File exported to {}", filename),
        Err(e) => println!("{}", e),
    }
}
