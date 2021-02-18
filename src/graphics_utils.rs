use nalgebra::{Point3, Vector3};

pub type ColorRGB = Vector3<f64>;

pub struct Ray {
  pub origin: Point3<f64>,
  pub direction: Vector3<f64>,
}

impl Ray {
  pub fn index(&self, t: f64) -> Point3<f64> {
    self.origin + t * self.direction
  }
}

pub struct Hit {
  pub t: f64,
  pub point: Point3<f64>,
  pub normal: Vector3<f64>,
  pub front_face: bool,
}

impl Hit {
  pub fn new(ray: &Ray, t: f64, outward_normal: Vector3<f64>) -> Hit {
    let front_face = ray.direction.dot(&outward_normal) < 0.;

    Hit {
      t,
      point: ray.index(t),
      normal: if front_face { outward_normal } else { -outward_normal },
      front_face,
    }
  }
}

pub trait Hittable {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn test_indexing() {
    let ray = Ray {
      origin: Point3::new(0., 0., 0.),
      direction: Vector3::new(0., 1., 0.),
    };

    assert_eq!(ray.index(1.), Point3::new(0., 1., 0.));
  }
}
