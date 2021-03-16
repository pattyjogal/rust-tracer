use nalgebra::{Point3, Vector3};
use float_ord::FloatOrd;

use std::mem::swap;

pub type ColorRGB = Vector3<f64>;

const N_DIMENSIONS: usize = 3;

/// A representation of a ray with a point and origin
pub struct Ray {
  pub origin: Point3<f64>,
  pub direction: Vector3<f64>,
}

impl Ray {
  /// Allow for indexing into the parameterized ray equation
  ///
  /// # Arguments
  /// `t` - The time along the ray; `0` at the origin
  ///
  /// # Returns
  /// The point at time `t` on the ray
  pub fn index(&self, t: f64) -> Point3<f64> {
    self.origin + t * self.direction
  }
}

pub struct Hit {
  /// The time along the ray where the hit occurred
  pub t: f64,
  /// The point of contact with the ray and the object
  pub point: Point3<f64>,
  /// The normal vector to the hit
  pub normal: Vector3<f64>,
  /// Whether or not the normal is facing out
  pub front_face: bool,
}

impl Hit {
  /// Constructs a hit from ray-intersection information
  ///
  /// # Arguments
  /// `ray` - A reference to the ray that generated the hit
  /// `t` - The time along the ray when the hit occurred
  /// `outward_normal` - The outward-facing normal vector of the intersected surface
  ///
  /// # Returns
  /// The constructed hit object
  pub fn new(ray: &Ray, t: f64, outward_normal: Vector3<f64>) -> Hit {
    let front_face = ray.direction.dot(&outward_normal) < 0.;

    Hit {
      t,
      point: ray.index(t),
      normal: if front_face {
        outward_normal
      } else {
        -outward_normal
      },
      front_face,
    }
  }
}

/// Represents an object that can be hit by a ray
pub trait Hittable {
  /// Determines if a ray hit the hittable object between the range of `t`s
  ///
  /// # Arguments
  /// `ray` - A reference to the ray used to intersect
  /// `t_min` - The minimum time along the ray to check
  /// `t_max` - The maximum time along the ray to check
  ///
  /// # Returns
  /// An `Option` of a possible closest hit
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;

  /// Defines a bounding box for the hittable object, if applicable
  ///
  /// # Arguments
  /// `t_start` - The start time for the render (used in animation)
  /// `t_end` - The end time for the render
  ///
  /// # Returns
  /// An `Option` of a AABB, if the Hittable would produce one
  fn get_bounding_box(&self, t_start: f64, t_end: f64) -> Option<AxisAlignedBoundingBox>;
}

#[derive(Copy, Clone)]
pub struct AxisAlignedBoundingBox {
  start: Point3<f64>,
  end: Point3<f64>,
}

impl AxisAlignedBoundingBox {
  pub fn new(start: Point3<f64>, end: Point3<f64>) -> AxisAlignedBoundingBox {
    AxisAlignedBoundingBox { start, end }
  }

  pub fn check_ray_hit(&self, ray: &Ray, mut t_min: f64, mut t_max: f64) -> bool {
    for i in 0..N_DIMENSIONS {
      let inverse_direction = 1. / ray.direction[i];
      let mut t0 = (self.start[i] - ray.origin[i]) * inverse_direction;
      let mut t1 = (self.end[i] - ray.origin[i]) * inverse_direction;
      if inverse_direction < 0. {
        swap(&mut t0, &mut t1);
      }

      t_min = if t0 > t_min { t0 } else { t_min };
      t_max = if t1 < t_max { t1 } else { t_max };

      if t_max <= t_min {
        return false;
      }
    }

    true
  }
}

pub fn compute_surrounding_box(first: &AxisAlignedBoundingBox, second: &AxisAlignedBoundingBox) -> AxisAlignedBoundingBox {
  let small = Point3::new(
    first.start.x.min(second.start.x),
    first.start.y.min(second.start.y),
    first.start.z.min(second.start.z),
  );

  let big = Point3::new(
    first.start.x.max(second.start.x),
    first.start.y.max(second.start.y),
    first.start.z.max(second.start.z),
  );

  AxisAlignedBoundingBox::new(small, big)
}

/// Computes a vector with magnitude 1 from another vector
///
/// # Arguments
/// `vec` - The vector to convert to a unit vector
///
/// # Returns
/// A unit vector
pub fn unit_vector(vec: Vector3<f64>) -> Vector3<f64> {
  vec / vec.norm()
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
