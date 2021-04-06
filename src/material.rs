use super::graphics_utils::{ColorRGB, Ray, Hit, unit_vector};
use nalgebra::Vector3;
use rand::Rng;

/// A representation of the object's shading parameters
pub trait Material {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(ColorRGB, Ray)>;
    fn albedo(&self) -> ColorRGB;
}

pub struct Lambertian {
  /// The color of the material
  pub albedo: ColorRGB,
}

impl Lambertian {
  pub fn new(albedo: ColorRGB) -> Self {
    Self {
      albedo,
    }
  }
}

impl Material for Lambertian {
  fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(ColorRGB, Ray)> {
    // let mut scatter_dir = hit.normal;// + random_unit_vector();
    let mut scatter_dir = hit.normal + random_unit_vector();
    if !hit.front_face {
      scatter_dir = - scatter_dir;
    }

    if near_zero(&scatter_dir) {
      scatter_dir = hit.normal;
    }

    Some((self.albedo, Ray { origin: hit.point, direction: scatter_dir }))
  }

  fn albedo(&self) -> ColorRGB {
    self.albedo
  }
}

pub struct Metal {
  /// The tint of the metal
  pub albedo: ColorRGB,
}

impl Metal {
  pub fn new(albedo: ColorRGB) -> Self {
    Self {
      albedo,
    }
  }
}

impl Material for Metal {
  fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(ColorRGB, Ray)> {
    let reflected = reflect(&unit_vector(ray.direction), &hit.normal);
    let scattered = Ray { origin: hit.point, direction: reflected };

    if reflected.dot(&hit.normal) > 0. {
      Some((self.albedo, scattered))
    } else {
      None
    }
  }

  fn albedo(&self) -> ColorRGB {
    self.albedo
  }
}

fn random_unit_vector() -> Vector3<f64> {
  let mut rng = rand::thread_rng();
  let theta = rng.gen_range(0.0..(2. * std::f64::consts::PI));
  let z: f64 = rng.gen_range(-1.0..1.);
  let unit_z = (1. - z.powf(2.)).sqrt();

  Vector3::new(unit_z * theta.cos(), unit_z * theta.sin(), z)
}

fn near_zero(vec: &Vector3<f64>) -> bool {
  let eps: f64 = 1e-6;
  vec[0].abs() < eps && vec[1].abs() < eps && vec[2].abs() < eps
}

fn reflect(v: &Vector3<f64>, n: &Vector3<f64>) -> Vector3<f64> {
  v - 2. * v.dot(n) * unit_vector(*n)
}