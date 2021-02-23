use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::path::Path;
use std::vec::Vec;

use nalgebra::{Point3, Vector3, distance};
use png;

use super::camera::Camera;
use super::graphics_utils::{unit_vector, ColorRGB, Hit, Hittable, Ray};

const EPSILON: f64 = 0.001;

pub struct RenderedScene {
  image_height: usize,
  image_width: usize,
  pixel_data: Vec<ColorRGB>,
  objects: Vec<Box<dyn Renderable>>,
  camera: Camera,
  mj_fine_grid_size: usize,
  light: PointLight,
}

impl RenderedScene {
  pub fn new(
    image_height: usize,
    image_width: usize,
    viewport_height: f64,
    viewport_width: f64,
    focal_length: f64,
    default_color: ColorRGB,
    objects: Vec<Box<dyn Renderable>>,
    camera: Camera,
    mj_fine_grid_size: usize,
    light: PointLight,
  ) -> Self {
    let origin = Vector3::new(0., 0., 0.);
    let horizontal = Vector3::new(viewport_width, 0., 0.);
    let vertical = Vector3::new(0., viewport_height, 0.);
    let lower_left_corner =
      origin - (horizontal / 2.) - (vertical / 2.) - Vector3::new(0., 0., focal_length);

    let pixel_data: Vec<ColorRGB> = vec![default_color; image_height * image_width];

    RenderedScene {
      image_height,
      image_width,
      pixel_data,
      objects,
      camera,
      mj_fine_grid_size,
      light,
    }
  }

  pub fn render(&mut self) {
    for j in 0..self.image_height {
      for i in 0..self.image_width {
        self.render_pixel(i, j);
      }
    }
  }

  fn render_pixel(&mut self, x: usize, y: usize) {
    let u = (x as f64) / (self.image_width as f64 - 1.);
    let v = (y as f64) / (self.image_height as f64 - 1.);

    // Multi-Jittered sampling
    for i in 0..self.mj_fine_grid_size {
      for j in 0..self.mj_fine_grid_size {
        let i_offset = (i as f64 / self.mj_fine_grid_size as f64) / (self.image_width as f64 - 1.);
        let j_offset = (j as f64 / self.mj_fine_grid_size as f64) / (self.image_height as f64 - 1.);
        let ray = self.camera.get_ray(u + i_offset, v + j_offset);
        match self.hit_objects(&ray, EPSILON, std::f64::INFINITY) {
          Some((hit, object)) => {
            // TODO: Only works if default color is black
            // let color = 0.5 * object.color_at_ray_hit(&ray);
            let color = object.calculate_shade_at_hit(
              object.material().color,
              &self.light,
              &hit,
              &self.objects,
            );
            let pixel_val = self.pixel_data[x + y * self.image_width];
            self.pixel_data[x + y * self.image_width] =
              pixel_val + color / self.mj_fine_grid_size.pow(2) as f64
          }
          None => {}
        }
      }
    }
  }

  fn hit_objects(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(Hit, &Box<dyn Renderable>)> {
    let mut possible_obj_hit: Option<(Hit, &Box<dyn Renderable>)> = None;
    let mut closest_so_far = t_max;

    for object in &self.objects {
      match object.check_ray_hit(ray, t_min, closest_so_far) {
        Some(hit) => {
          closest_so_far = hit.t;
          possible_obj_hit = Some((hit, object));
        }
        None => {}
      }
    }

    possible_obj_hit
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
  fn color_at_ray_hit(&self, hit: &Ray) -> ColorRGB;
  fn material(&self) -> &Material;
}

pub trait Renderable: Colorable + Hittable {
  fn calculate_shade_at_hit(
    &self,
    ambient_color: ColorRGB,
    light: &PointLight,
    hit: &Hit,
    objects: &Vec<Box<dyn Renderable>>,
  ) -> ColorRGB {
    // Calculate shading
    let ambient = self.material().k_ambient * ambient_color;
    let diffuse = self.material().k_diffuse
      * Vector3::from(light.point - hit.point).dot(&hit.normal)
      * light.color;
    let shaded_color = self.material().color.component_mul(&(ambient + diffuse));

    // Calculate shadow
    let ray_to_light = Ray {
      origin: hit.point,
      direction: Vector3::from(light.point - hit.point),
    };
    for object in objects {
      match object.check_ray_hit(&ray_to_light, 0.00015, 1.0) {
        Some(hit) => return shaded_color - ColorRGB::new(0.3, 0.3, 0.3),
        None => {}
      }
    }

    shaded_color
  }
}

/// Plane
pub struct Plane {
  pub point: Point3<f64>,
  pub normal: Vector3<f64>,
  pub material: Material,
}

impl Renderable for Plane {}

impl Hittable for Plane {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
    let Ray { origin, direction } = ray;
    let t = (self.point - origin).dot(&self.normal) / direction.dot(&self.normal);

    if t.is_nan() || t < t_min || t > t_max {
      None
    } else {
      Some(Hit::new(ray, t, self.normal))
    }
  }
}

impl Colorable for Plane {
  fn color_at_ray_hit(&self, ray: &Ray) -> ColorRGB {
    return ColorRGB::new(0.8, 0.8, 0.8);
  }

  fn material(&self) -> &Material {
    &self.material
  }
}

/// Sphere
pub struct Sphere {
  pub center: Point3<f64>,
  pub radius: f64,
  pub material: Material,
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

impl Renderable for Sphere {}

impl Hittable for Sphere {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
    let discriminant = self.calc_discriminant(ray);
    if discriminant < 0. {
      return None;
    }
    let sqrtd = discriminant.sqrt();

    if distance(&ray.origin, &self.center) < self.radius {
      return None;
    }

    // Find nearest root in acceptable range
    let translated_origin = ray.origin - self.center;
    let a = ray.direction.dot(&ray.direction);
    let hb = translated_origin.dot(&ray.direction);

    let mut root = (-hb - sqrtd) / a;
    if root < t_min || t_max < root {
      root = (-hb + sqrtd) / a;
      if root < t_min || t_max < root {
        return None;
      }
    }

    let point = ray.index(root);
    let outward_normal = (point - self.center) / self.radius;

    // TODO: Maybe move Hit to Ray?
    Some(Hit::new(ray, root, outward_normal))
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

  fn material(&self) -> &Material {
    &self.material
  }
}

/// Triangle
pub struct Triangle {
  pub p0: Point3<f64>,
  pub p1: Point3<f64>,
  pub p2: Point3<f64>,
  pub material: Material,
}

impl Renderable for Triangle {}

impl Hittable for Triangle {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
    let e1 = Vector3::from(self.p1 - self.p0);
    let e2 = Vector3::from(self.p2 - self.p0);
    let h = ray.direction.cross(&e2);
    let a = e1.dot(&h);

    if a > -EPSILON && a < EPSILON {
      return None;
    }

    let f = 1.0 / a;
    let s = ray.origin - self.p0;
    let u = f * s.dot(&h);
    if u < 0. || u > 1. {
      return None;
    }

    let q = s.cross(&e1);
    let v = f * ray.direction.dot(&q);
    if v < 0. || u + v > 1. {
      return None;
    }

    let t = f * e2.dot(&q);
    if t > EPSILON {
      return Some(Hit::new(&ray, t, e1.cross(&e2)));
    }

    None
  }
}

impl Colorable for Triangle {
  fn color_at_ray_hit(&self, ray: &Ray) -> ColorRGB {
    return ColorRGB::new(0.8, 0.8, 0.8);
  }

  fn material(&self) -> &Material {
    &self.material
  }
}

// Point Light
pub struct PointLight {
  pub point: Point3<f64>,
  pub color: ColorRGB,
}

#[derive(Clone, Copy)]
pub struct Material {
  k_diffuse: f64,
  k_ambient: f64,
  color: ColorRGB,
}

impl Material {
  pub fn new(k_diffuse: f64, k_ambient: f64, color: ColorRGB) -> Material {
    Material {
      k_diffuse,
      k_ambient,
      color,
    }
  }
}
