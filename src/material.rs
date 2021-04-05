use super::graphics_utils::ColorRGB;

/// A representation of the object's shading parameters
#[derive(Clone, Copy)]
pub struct Material {
  /// Diffuse constant term from Phong equation
  pub k_diffuse: f64,
  /// Ambient constant term from Phone equation
  pub k_ambient: f64,
  /// The color of the material
  pub color: ColorRGB,
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