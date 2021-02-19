use nalgebra::{Point3, Vector3};

use super::graphics_utils::{unit_vector, Ray};

pub enum CameraMode {
  Perspective,
  Orthographic,
}

pub struct Camera {
  origin: Point3<f64>,
  horizontal: Vector3<f64>,
  vertical: Vector3<f64>,
  lower_left_corner: Point3<f64>,
  mode: CameraMode,
}

impl Camera {
  pub fn new(
    look_from: Point3<f64>,
    look_at: Point3<f64>,
    up: Vector3<f64>,
    fov: f64,
    aspect_ratio: f64,
    mode: CameraMode,
  ) -> Camera {
    let theta = fov * std::f64::consts::PI / 180.;
    let h = (theta / 2.).tan();
    let viewport_height = 2. * h;
    let viewport_width = aspect_ratio * viewport_height;

    let w = unit_vector(look_from - look_at);
    let u = unit_vector(up.cross(&w));
    let v = w.cross(&u);

    let origin = look_from;
    let horizontal = viewport_width * u;
    let vertical = viewport_height * v;

    Camera {
      origin,
      horizontal,
      vertical,
      lower_left_corner: origin - horizontal / 2. - vertical / 2. - w,
      mode,
    }
  }

  pub fn get_ray(&self, x: f64, y: f64) -> Ray {
    let pixel_world_coords = self.lower_left_corner + x * self.horizontal + y * self.vertical;
    match &self.mode {
      CameraMode::Perspective => Ray {
        origin: self.origin,
        direction: pixel_world_coords - self.origin,
      },
      CameraMode::Orthographic => Ray {
        origin: Point3::from(pixel_world_coords),
        direction: Vector3::new(0., 0., -1.),
      }
    }
  }
}
