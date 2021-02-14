use nalgebra::Vector3;
use nalgebra::Point3;

struct Ray {
  origin: Point3,
  direction: Vector3,
}

impl Index<usize> for Ray {
  type Output = Point3;

  fn index(&self, t: f64) -> &Self::Output {
    self.origin + t * self.direction;
  }
}