use nalgebra::{Point3, Vector3};

use ordered_float::NotNan;
use rand::Rng;
use std::rc::Rc;

mod camera;
mod graphics_utils;
mod material;
mod renderer;

// Image
// const ASPECT_RATIO: f64 = 1.;
const ASPECT_RATIO: f64 = 16. / 9.;
const IMAGE_WIDTH: usize = 400;
const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;

fn generate_spheres_in_box(
    n: usize,
    radius: f64,
    start: Point3<f64>,
    end: Point3<f64>,
) -> Vec<Rc<dyn renderer::Renderable>> {
    let mut ret: Vec<Rc<dyn renderer::Renderable>> = vec![];
    let mut rng = rand::thread_rng();
    let mat = Rc::new(material::Lambertian::new(
        graphics_utils::ColorRGB::new(1., rng.gen::<f64>(), rng.gen::<f64>()),
    ));

    for _ in 0..n {
        ret.push(Rc::new(renderer::Sphere {
            radius,
            center: random_point_in_box(start, end),
            material: mat.clone(),
        }))
    }

    ret
}

fn generate_mesh(filename: &str) -> Vec<Rc<dyn renderer::Renderable>> {
    let mesh = renderer::Mesh::from_obj(filename).expect("Could not load file");
    vec![Rc::new(mesh)]
}

fn random_point_in_box(start: Point3<f64>, end: Point3<f64>) -> Point3<f64> {
    let mut rng = rand::thread_rng();
    Point3::new(
        rng.gen_range(start.x..end.x),
        rng.gen_range(start.y..end.y),
        rng.gen_range(start.z..end.z),
    )
}

fn main() {
    // Render
    let camera = camera::Camera::new(
        Point3::new(0., 1.0, -1.),
        Point3::new(0., 0., 3.),
        Vector3::new(0., -1., 0.),
        45.,
        ASPECT_RATIO,
        camera::CameraMode::Perspective,
    );

    let diffuse_mat = Rc::new(material::Lambertian::new(graphics_utils::ColorRGB::new(
        0.7,
        0.3,
        0.3,
    )));
    let ground_mat = Rc::new(material::Lambertian::new(graphics_utils::ColorRGB::new(
        0.8,
        0.8,
        0.0,
    )));
    let glass_mat = Rc::new(material::Dielectric::new(1.5));
    let metal_mat = Rc::new(material::Metal::new(graphics_utils::ColorRGB::new(
        0.75, 0.75, 0.75,
    )));

    let light = renderer::PointLight {
        point: Point3::new(2., 4., -3.),
        color: graphics_utils::ColorRGB::new(0.9, 0.9, 0.9),
    };

    let objects: Vec<Rc<dyn renderer::Renderable>> = vec![
        // Rc::new(renderer::Plane {
        //     point: Point3::new(0., 0., 0.),
        //     normal: Vector3::new(0., 1., 0.),
        //     material: gray_mat,
        // }),
        // Large pink pyramid
        // Rc::new(renderer::Triangle {
        //     p0: Point3::new(-0.1, 0.0, 2.1).map(|x| NotNan::new(x).unwrap()),
        //     p1: Point3::new(0.2, 0.0, 2.5).map(|x| NotNan::new(x).unwrap()),
        //     p2: Point3::new(-0.07, 2.0, 2.15).map(|x| NotNan::new(x).unwrap()),
        //     material: red_mat,
        // }),
        // Rc::new(renderer::Triangle {
        //     p0: Point3::new(-0.1, 0.0, 2.1).map(|x| NotNan::new(x).unwrap()),
        //     p1: Point3::new(-0.3, 0.0, 2.5).map(|x| NotNan::new(x).unwrap()),
        //     p2: Point3::new(-0.07, 2.0, 2.15).map(|x| NotNan::new(x).unwrap()),
        //     material: red_mat,
        // }),
        // Rc::new(renderer::Triangle {
        //     p1: Point3::new(0.2, 0.0, 2.5).map(|x| NotNan::new(x).unwrap()),
        //     p0: Point3::new(-0.3, 0.0, 2.5).map(|x| NotNan::new(x).unwrap()),
        //     p2: Point3::new(-0.07, 2.0, 2.15).map(|x| NotNan::new(x).unwrap()),
        //     material: red_mat,
        // }),
        Rc::new(renderer::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.5,
            material: glass_mat.clone(),
        }),
        Rc::new(renderer::Sphere {
            center: Point3::new(-1., 0., 0.),
            radius: 0.5,
            material: metal_mat.clone(),
        }),
        Rc::new(renderer::Sphere {
            center: Point3::new(1., 0., 0.),
            radius: 0.5,
            material: metal_mat.clone(),
        }),
        Rc::new(renderer::Sphere {
            center: Point3::new(0., -100.5, 0.),
            radius: 100.,
            material: diffuse_mat.clone(),
        }),
        Rc::new(renderer::Sphere {
            center: Point3::new(0., 0., 4.),
            radius: 1.0,
            material: diffuse_mat.clone(),
        }),
        // Rc::new(renderer::Sphere {
        //     center: Point3::new(1., 0., 3.),
        //     radius: 0.5,
        //     material: material::Material::new(
        //         0.01,
        //         0.99,
        //         graphics_utils::ColorRGB::new(255. / 255., 195. / 255., 18. / 255.),
        //     ),
        // }),
        // Rc::new(renderer::Sphere {
        //     center: Point3::new(-1., 0., 3.),
        //     radius: 0.5,
        //     material: Rc::new(material::Lambertian::new(
        //         graphics_utils::ColorRGB::new(255. / 255., 195. / 255., 18. / 255.),
        //     )),
        // }),
        // Rc::new(renderer::Sphere {
        //     center: Point3::new(0., 0., -3.),
        //     radius: 0.5,
        //     material: material::Material::new(
        //         0.01,
        //         0.99,
        //         graphics_utils::ColorRGB::new(255. / 255., 195. / 255., 18. / 255.),
        //     ),
        // }),
    ];

    let mut scene = renderer::RenderedScene::new(
        IMAGE_HEIGHT,
        IMAGE_WIDTH,
        graphics_utils::ColorRGB::new(0., 0., 0.),
        objects,
        // generate_spheres_in_box(100000, 0.05, Point3::new(-2., -2., 5.), Point3::new(2., 2., 6.)),
        // generate_mesh("./teapot.obj"),
        camera,
        4,
        light,
        renderer::HitDetectionMode::Naive,
    );
    scene.render();

    // Export
    let filename = format!("./test.png");
    match scene.export_to_png(&filename) {
        Ok(()) => println!("File exported to {}", filename),
        Err(e) => println!("{}", e),
    }
}
