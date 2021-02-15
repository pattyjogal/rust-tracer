use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::path::Path;
use std::vec::Vec;

use nalgebra::{Point3, Vector3};
use png;

use super::graphics_utils::{ColorRGB, Ray, Hit, Hittable};

pub struct RenderedScene {
  aspect_ratio: f64,
  image_height: usize,
  image_width: usize,
  viewport_height: f64,
  viewport_width: f64,
  focal_length: f64,
  pixel_data: Vec<ColorRGB>,
  origin: Vector3<f64>,
  horizontal: Vector3<f64>,
  vertical: Vector3<f64>,
  lower_left_corner: Vector3<f64>,
}

impl RenderedScene {
  pub fn new(
    aspect_ratio: f64,
    image_height: usize,
    image_width: usize,
    viewport_height: f64,
    viewport_width: f64,
    focal_length: f64,
    default_color: ColorRGB,
  ) -> Self {
    let origin = Vector3::new(0., 0., 0.);
    let horizontal = Vector3::new(viewport_width, 0., 0.);
    let vertical = Vector3::new(0., viewport_height, 0.);
    let lower_left_corner =
      origin - (horizontal / 2.) - (vertical / 2.) - Vector3::new(0., 0., focal_length);

    println!("{}\n{}\n{}", horizontal, vertical, lower_left_corner);
    let pixel_data: Vec<ColorRGB> = vec![default_color; image_height * image_width];

    RenderedScene {
      aspect_ratio,
      image_height,
      image_width,
      viewport_height,
      viewport_width,
      focal_length,
      pixel_data,
      origin,
      horizontal,
      vertical,
      lower_left_corner,
    }
  }

  pub fn render_object(&mut self, object: &impl Colorable) {
    for j in (0..self.image_height) {
      for i in 0..self.image_width {
        self.render_pixel(i, j, object);
      }
    }
  }

  fn render_pixel(&mut self, x: usize, y: usize, object: &impl Colorable) {
    let u = (x as f64) / (self.image_width as f64 - 1.);
    let v = (y as f64) / (self.image_height as f64 - 1.);
    let ray = Ray {
      origin: Point3::from(self.origin),
      direction: self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
    };

    if object.check_ray_hit(&ray) {
      self.pixel_data[x + y * self.image_width] = object.color_at_ray_hit(&ray);
    }
  }

  pub fn export_to_png(&self, filename: &str) -> Result<(), io::Error> {
    let path = Path::new(filename);
    let file = File::create(path)?;
    let ref mut w = BufWriter::new(&file);

    let mut encoder = png::Encoder::new(w, self.image_width as u32, self.image_height as u32);
    encoder.set_color(png::ColorType::RGB);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;

    let data = self
      .pixel_data
      .iter()
      .flat_map(ColorRGB::as_slice)
      .map(|n| (n * 255.) as u8)
      .collect::<Vec<u8>>();
    writer.write_image_data(&data[..])?;

    Ok(())
  }
}

pub trait Colorable {
  fn color_at_ray_hit(&self, ray: &Ray) -> ColorRGB;
}

/// Plane
pub struct Plane {
  pub top_left_corner: Point3<f64>,
  pub bottom_right_corner: Point3<f64>,
}

impl Hittable for Plane {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit: &Hit) -> bool {
    true
  }
}

impl Colorable for Plane {
  fn color_at_ray_hit(&self, ray: &Ray) -> ColorRGB {
    let unit_direction = ray.direction / (ray.direction.norm() as f64);
    let t = 0.5 * (unit_direction.y + 1.0);

    (1.0 - t) * ColorRGB::new(1., 1., 1.) + t * ColorRGB::new(0.5, 0.7, 1.0)
  }
}

/// Sphere
pub struct Sphere {
  pub center: Point3<f64>,
  pub radius: f64,
}

impl Sphere {
  fn calc_discriminant(&self, ray: &Ray) -> f64 {
    let translated_origin = ray.origin - self.center;
    let a = ray.direction.dot(&ray.direction);
    let b = 2.0 * translated_origin.dot(&ray.direction);
    let c = translated_origin.dot(&translated_origin) - self.radius * self.radius;

    b * b - 4. * a * c
  }
}

impl Hittable for Sphere {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit: &Hit) -> bool {
    self.calc_discriminant(ray) > 0.
  }
}

impl Colorable for Sphere {
  fn color_at_ray_hit(&self, ray: &Ray) -> ColorRGB {
    let translated_origin = ray.origin - self.center;
    let discriminant = self.calc_discriminant(&ray);
    let a = ray.direction.dot(&ray.direction);
    let b = 2.0 * translated_origin.dot(&ray.direction);
    let t = (-b - discriminant.sqrt()) / (2.0 * a);

    let z_scaled = ray.index(t) - Vector3::new(0., 0., -1.);
    let n = z_scaled / z_scaled.coords.norm();
    
    0.5 * ColorRGB::new(n.x + 1., n.y + 1., n.z + 1.)
  }
}
