use nalgebra::{Point3, Vector3};

pub type ColorRGB = Vector3<f64>;

pub struct Ray {
  pub origin: Point3<f64>,
  pub direction: Vector3<f64>,
}

impl Ray {
  fn index(&self, t: f64) -> Point3<f64> {
    self.origin + t * self.direction
  }

  pub fn color(&self) -> ColorRGB {
    let unit_direction = self.direction / (self.direction.norm() as f64);
    let t = 0.5 * (unit_direction.y + 1.0);

    (1.0 - t) * ColorRGB::new(1., 1., 1.) + t * ColorRGB::new(0.5, 0.7, 1.0)
  }
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
