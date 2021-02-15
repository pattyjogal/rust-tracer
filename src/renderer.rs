use std::vec::Vec;
use std::path::Path;
use std::fs::File;
use std::io;
use std::io::{BufWriter};

use png;
use nalgebra::{Point3, Vector3};

use super::graphics_utils::{ColorRGB, Ray};

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
      origin - horizontal / 2. - vertical / 2. - Vector3::new(0., 0., focal_length);

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

  pub fn add_pixel(&mut self, x: usize, y: usize) {
    let u = (x as f64) / (self.image_width as f64 - 1.);
    let v = (y as f64) / (self.image_width as f64 - 1.);
    let ray = Ray {
      origin: Point3::from(self.origin),
      direction: self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
    };

    self.pixel_data[x + y * self.image_width] = ray.color();
  }

  pub fn export_to_png(&self, filename: &str) -> Result<(), io::Error> {
    let path = Path::new(filename);
    let file = File::create(path)?;
    let ref mut w = BufWriter::new(&file);

    let mut encoder = png::Encoder::new(w, self.image_width as u32, self.image_height as u32);
    encoder.set_color(png::ColorType::RGB);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;

    let data = self.pixel_data.iter().flat_map(ColorRGB::as_slice).map(|n| (n * 255.) as u8).collect::<Vec<u8>>();
    writer.write_image_data(&data[..])?;

    Ok(())
  }
}
