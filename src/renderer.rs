use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::path::Path;
use std::rc::Rc;
use std::time::Instant;
use std::vec::Vec;

use nalgebra::{distance, Point3, Vector3};
use png;
use rand::Rng;

use super::camera::Camera;
use super::graphics_utils::{
  compute_surrounding_box, AxisAlignedBoundingBox, ColorRGB, Hit, Hittable, Ray,
};

/// A small value to offset some ray origins in calculations
const EPSILON: f64 = 0.001;

pub enum HitDetectionMode {
  Naive,
  Bvh,
}

pub enum NodeType<'a> {
  Node(&'a BVHNode),
  Primitive,
}

/// A representation of a 3D space
pub struct RenderedScene {
  /// The output screen's height in pixels
  image_height: usize,
  /// The output screen's width in pixels
  image_width: usize,
  /// A single-dimensional vector of pixel colors, flattened on rows
  pixel_data: Vec<ColorRGB>,
  /// A vector containing all the renderable objects in the scene
  objects: Vec<Rc<dyn Renderable>>,
  /// The camera used to look into the scene
  camera: Camera,
  /// The size of the fine multi-jittered grid (where the coarse grid is the square root of this number)
  mj_fine_grid_size: usize,
  /// A point light to illuminate the scene
  light: PointLight,
  hit_detection_mode: HitDetectionMode,
  root_bvh_node: Option<BVHNode>,
}

impl RenderedScene {
  /// Constructs a rendered scene
  ///
  /// # Arguments
  /// `image_height` - The height of the output image in pixels
  /// `image_width` - The width of the output image in pixels
  /// `default_color` - The color to show in the background when no object is hit by a ray
  /// `objects` - A collection of objects to display in the scene
  /// `camera` - The camera used to generate the render of the scene
  /// `mj_find_grid_size` - The side length of the fine grid used in Multi-Jittered Sampling
  /// `light` - A point light to illuminate the scene
  ///
  /// # Returns
  /// An initialized scene
  pub fn new(
    image_height: usize,
    image_width: usize,
    default_color: ColorRGB,
    objects: Vec<Rc<dyn Renderable>>,
    camera: Camera,
    mj_fine_grid_size: usize,
    light: PointLight,
    hit_detection_mode: HitDetectionMode,
  ) -> Self {
    let pixel_data: Vec<ColorRGB> = vec![default_color; image_height * image_width];
    let root_bvh_node = match &hit_detection_mode {
      HitDetectionMode::Bvh => Some(BVHNode::new(&objects, 0, objects.len())),
      _ => None,
    };

    RenderedScene {
      image_height,
      image_width,
      pixel_data,
      objects,
      camera,
      mj_fine_grid_size,
      light,
      hit_detection_mode,
      root_bvh_node,
    }
  }

  /// Computes the colors at each pixel, storing them for exporting
  pub fn render(&mut self) {
    let now = Instant::now();
    for i in 0..self.image_width {
      for j in 0..self.image_height {
        self.render_pixel(i, j);
      }
    }
  }

  /// Determines the color of a single pixel and stores it in the scene's rendered output
  ///
  /// # Arguments
  /// `x` - The horizontal coordinate of the pixel, where left is `0`
  /// `y` - The vertical coordinate of the pixel, where bottom is `0`
  fn render_pixel(&mut self, x: usize, y: usize) {
    // Adjusted viewplane world-coordinates
    let u = (x as f64) / (self.image_width as f64 - 1.);
    let v = (y as f64) / (self.image_height as f64 - 1.);

    // Sets to determine constraints of MJS
    let mut used_rows = HashSet::<usize>::new();
    let mut used_cols = HashSet::<usize>::new();
    let mut used_squares = HashSet::<(usize, usize)>::new();

    // Multi-Jittered sampling
    for i in 0..self.mj_fine_grid_size {
      if used_rows.contains(&i) {
        continue;
      }

      for j in 0..self.mj_fine_grid_size {
        if used_cols.contains(&j) {
          continue;
        }

        // The current MJ subgrid index pair
        let square = (
          i / (self.mj_fine_grid_size as f64).sqrt() as usize,
          j / (self.mj_fine_grid_size as f64).sqrt() as usize,
        );
        if used_squares.contains(&square) {
          continue;
        }

        used_rows.insert(i);
        used_cols.insert(j);
        used_squares.insert(square);

        // The subpixel offset for the chosen MJ square
        let i_offset = (i as f64 / self.mj_fine_grid_size as f64) / (self.image_width as f64 - 1.);
        let j_offset = (j as f64 / self.mj_fine_grid_size as f64) / (self.image_height as f64 - 1.);

        let ray = self.camera.get_ray(u + i_offset, v + j_offset);
        let now = Instant::now();
        let hit_objects = self.hit_objects(&ray, EPSILON, std::f64::INFINITY);
        match hit_objects {
          Some((hit, object)) => {
            // TODO: Only works for now if default color is black
            let color = object.calculate_shade_at_hit(
              &object.material().color,
              &self.light,
              &hit,
              &self.objects,
            );
            let pixel_val = self.pixel_data[x + y * self.image_width];
            self.pixel_data[x + y * self.image_width] =
              pixel_val + color / self.mj_fine_grid_size as f64
          }
          None => {}
        }
      }
    }
  }

  /// Given a ray and a range of times along the ray, attempts to find the closest object that ray hits, if any
  ///
  /// # Arguments
  /// `ray` - A reference to the ray used to hit the objects
  /// `t_min` - The minimum time along the ray to check
  /// `t_max` - The maxumum time along the ray to check
  ///
  /// # Returns
  /// A possible pairing of a hit record and the closest object that was hit
  fn hit_objects(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(Hit, &Rc<dyn Renderable>)> {
    match self.hit_detection_mode {
      HitDetectionMode::Naive => {
        let mut possible_obj_hit: Option<(Hit, &Rc<dyn Renderable>)> = None;
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
      HitDetectionMode::Bvh => match self
        .root_bvh_node
        .as_ref()
        .unwrap()
        .get_intersected_object(ray, t_min, t_max)
      {
        Some(object) => Some((Hit::new(ray, 0., Vector3::new(0., 0., 0.)), object)),
        None => None,
      },
    }
  }

  /// Writes the stored pixel data out to a PNG image
  ///
  /// Must be called after `RenderedScene::render()`
  ///
  /// # Arguments
  /// `filename` - The filepath to write the image to
  ///
  /// # Returns
  /// The result of the image-wiriting operation
  ///
  /// # Errors
  /// Produces an error type if there were any I/O issues while
  /// attempting to write to the given file
  pub fn export_to_png(&self, filename: &str) -> Result<(), io::Error> {
    println!("Exporting...");
    let path = Path::new(filename);
    let file = File::create(path)?;
    let ref mut w = BufWriter::new(&file);

    // Create the PNG encoder
    let mut encoder = png::Encoder::new(w, self.image_width as u32, self.image_height as u32);
    encoder.set_color(png::ColorType::RGB);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;

    // Translate the color data to byte values
    let data = self
      .pixel_data
      .iter()
      .flat_map(ColorRGB::as_slice)
      .map(|n| (n * 255.) as u8)
      .collect::<Vec<u8>>();
    writer.write_image_data(&data[..])?;

    Ok(())
  }

  fn get_bounding_box(&self, start_time: f64, end_time: f64) -> Option<AxisAlignedBoundingBox> {
    let mut possible_output_box = None;

    for object in &self.objects {
      match object.get_bounding_box(start_time, end_time) {
        Some(iter_box) => {
          possible_output_box = match possible_output_box {
            Some(output_box) => Some(compute_surrounding_box(&output_box, &iter_box)),
            None => Some(iter_box),
          }
        }
        None => {}
      }
    }

    possible_output_box
  }
}

/// Objects that can produce a color when hit by a ray
pub trait Colorable {
  /// A getter for the object's material
  ///
  /// # Returns
  /// A reference to the object's material
  fn material(&self) -> &Material;
}

/// Objects that can be shown in the scene
///
/// Must be both `Colorable` and `Hittable`, as items need to be
/// detected and properly colored in to be shown
pub trait Renderable: Colorable + Hittable {
  /// Given that a ray has intersected this object, compute the exact shade that the
  /// object's material will show. Factors in shadows and shading.
  ///
  /// # Arguments
  /// `ambient_color` - The color used for the ambient term in the Phong equation
  /// `light` - The light to use for checking shading and shadows
  /// `hit` - The ray hit record that produced this color
  /// `objects` - The objects in the scene that can product a shadow against this object
  ///
  /// # Returns
  /// The color shown at this hit point
  fn calculate_shade_at_hit(
    &self,
    ambient_color: &ColorRGB,
    light: &PointLight,
    hit: &Hit,
    objects: &Vec<Rc<dyn Renderable>>,
  ) -> ColorRGB {
    // Calculate shading
    let ambient = self.material().k_ambient * ambient_color;
    let diffuse = self.material().k_diffuse
      * Vector3::from(light.point - hit.point).dot(&hit.normal)
      * light.color;
    let shaded_color = self.material().color.component_mul(&(ambient + diffuse));

    // Calculate shadow
    // let ray_to_light = Ray {
    //   origin: hit.point,
    //   direction: Vector3::from(light.point - hit.point),
    // };

    // TODO: Maybe this slows it a bit?
    // for object in objects {
    //   match object.check_ray_hit(&ray_to_light, 0.015, 1.0) {
    //     Some(_hit) => return shaded_color - ColorRGB::new(0.3, 0.3, 0.3),
    //     None => {}
    //   }
    // }

    shaded_color
  }

  fn node_type(&self) -> NodeType;
}

/// Plane
pub struct Plane {
  /// A point that lies on the plane
  pub point: Point3<f64>,
  /// A normal vector for the plane
  pub normal: Vector3<f64>,
  /// The material used to shade the plane
  pub material: Material,
}

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

  fn get_bounding_box(&self, _: f64, _: f64) -> Option<AxisAlignedBoundingBox> {
    None
  }
}

impl Colorable for Plane {
  fn material(&self) -> &Material {
    &self.material
  }
}

impl Renderable for Plane {
  fn node_type(&self) -> NodeType {
    NodeType::Primitive
  }
}

/// Sphere
pub struct Sphere {
  /// The center point of the sphere
  pub center: Point3<f64>,
  /// The radius of the sphere
  pub radius: f64,
  /// The material used to shade the sphere
  pub material: Material,
}

impl Sphere {
  /// Helper to calculate the discriminant for root-finding in the sphere
  ///
  /// # Arguments
  /// `ray` - The reference ray to use for calculating the discriminant
  ///
  /// # Returns
  /// The discriminant
  fn calc_discriminant(&self, ray: &Ray) -> f64 {
    let translated_origin = ray.origin - self.center;
    let a = ray.direction.dot(&ray.direction);
    let b = 2.0 * translated_origin.dot(&ray.direction);
    let c = translated_origin.dot(&translated_origin) - self.radius * self.radius;

    b * b - 4. * a * c
  }
}

impl Hittable for Sphere {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
    let discriminant = self.calc_discriminant(ray);
    if discriminant < 0. {
      return None;
    }
    let sqrtd = discriminant.sqrt();

    // Disallow a hit starting from within the sphere
    if distance(&ray.origin, &self.center) < self.radius {
      return None;
    }

    // Find nearest root in acceptable range
    let translated_origin = ray.origin - self.center;
    let a = ray.direction.dot(&ray.direction);
    let hb = translated_origin.dot(&ray.direction);

    // Check both roots to see if a hit occurred at all
    let mut root = (-hb - sqrtd) / a;
    if root < t_min || t_max < root {
      root = (-hb + sqrtd) / a;
      if root < t_min || t_max < root {
        return None;
      }
    }

    let point = ray.index(root);
    let outward_normal = (point - self.center) / self.radius;

    Some(Hit::new(ray, root, outward_normal))
  }

  fn get_bounding_box(&self, _: f64, _: f64) -> Option<AxisAlignedBoundingBox> {
    Some(AxisAlignedBoundingBox::new(
      self.center - Vector3::new(self.radius, self.radius, self.radius),
      self.center + Vector3::new(self.radius, self.radius, self.radius),
    ))
  }
}

impl Colorable for Sphere {
  fn material(&self) -> &Material {
    &self.material
  }
}

impl Renderable for Sphere {
  fn node_type(&self) -> NodeType {
    NodeType::Primitive
  }
}

/// Triangle
pub struct Triangle {
  /// The first vertex
  pub p0: Point3<f64>,
  /// The second vertex
  pub p1: Point3<f64>,
  /// The third vertex
  pub p2: Point3<f64>,
  /// The material used to shade the triangle
  pub material: Material,
}

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
    if t > t_min && t < t_max {
      return Some(Hit::new(&ray, t, e1.cross(&e2)));
    }

    None
  }
  fn get_bounding_box(&self, _: f64, _: f64) -> Option<AxisAlignedBoundingBox> {
    None
  }
}

impl Colorable for Triangle {
  fn material(&self) -> &Material {
    &self.material
  }
}

impl Renderable for Triangle {
  fn node_type(&self) -> NodeType {
    NodeType::Primitive
  }
}

/// A light that shines isotropically from a point
pub struct PointLight {
  /// The location of the light
  pub point: Point3<f64>,
  /// The light's color
  pub color: ColorRGB,
}

/// A representation of the object's shading parameters
#[derive(Clone, Copy)]
pub struct Material {
  /// Diffuse constant term from Phong equation
  k_diffuse: f64,
  /// Ambient constant term from Phone equation
  k_ambient: f64,
  /// The color of the material
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

struct BVHNode {
  bounding_box: AxisAlignedBoundingBox,
  left_child: Rc<dyn Renderable>,
  right_child: Rc<dyn Renderable>,
  dummy_material: Material,
}

fn compare_boxes(a: &Rc<dyn Renderable>, b: &Rc<dyn Renderable>, axis: usize) -> bool {
  match (a.get_bounding_box(0., 0.), b.get_bounding_box(0., 0.)) {
    (Some(a_box), Some(b_box)) => a_box.start[axis] < b_box.start[axis],
    _ => false,
  }
}

impl BVHNode {
  fn new(raw_objects: &Vec<Rc<dyn Renderable>>, start: usize, end: usize) -> BVHNode {
    let mut rng = rand::thread_rng();
    let axis = rng.gen_range(0..=2);
    let comparator = |a, b| compare_boxes(a, b, axis);
    let object_span = end - start;
    let mut objects = raw_objects.to_vec();
    let left;
    let right;

    match object_span {
      1 => {
        left = objects[start].clone();
        right = objects[start].clone();
      }
      2 => {
        if comparator(&objects[start], &objects[start + 1]) {
          left = objects[start].clone();
          right = objects[start + 1].clone();
        } else {
          left = objects[start + 1].clone();
          right = objects[start].clone();
        }
      }
      _ => {
        objects.sort_by(|a, b| {
          if compare_boxes(a, b, axis) {
            std::cmp::Ordering::Less
          } else {
            std::cmp::Ordering::Greater
          }
        });
        let mid = start + object_span / 2;
        left = Rc::new(BVHNode::new(&objects, start, mid));
        right = Rc::new(BVHNode::new(&objects, mid, end));
      }
    }

    let left_box = left
      .get_bounding_box(0., 0.)
      .expect("Attempted to make BVH with unboundable object");
    let right_box = right
      .get_bounding_box(0., 0.)
      .expect("Attempted to make BVH with unboundable object");

    BVHNode {
      bounding_box: compute_surrounding_box(&left_box, &right_box),
      left_child: left,
      right_child: right,
      dummy_material: Material {
        color: ColorRGB::new(0., 0., 0.),
        k_ambient: 0.,
        k_diffuse: 0.,
      },
    }
  }

  fn get_intersected_object(
    &self,
    ray: &Ray,
    t_min: f64,
    t_max: f64,
  ) -> Option<&Rc<dyn Renderable>> {
    if !self.bounding_box.check_ray_hit(ray, t_min, t_max) {
      return None;
    }

    let left_res = match self.left_child.node_type() {
      NodeType::Node(node) => node.get_intersected_object(ray, t_min, t_max),
      NodeType::Primitive => match self.left_child.check_ray_hit(ray, t_min, t_max) {
        Some(_) => Some(&self.left_child),
        None => None,
      },
    };

    if left_res.is_some() {
      return left_res;
    }

    let right_res = match self.right_child.node_type() {
      NodeType::Node(node) => node.get_intersected_object(ray, t_min, t_max),
      NodeType::Primitive => match self.right_child.check_ray_hit(ray, t_min, t_max) {
        Some(_) => Some(&self.right_child),
        None => None,
      },
    };

    if right_res.is_some() {
      return right_res;
    }

    None
  }

  fn get_bounding_box(&self, _: f64, _: f64) -> Option<AxisAlignedBoundingBox> {
    Some(self.bounding_box)
  }
}

impl Hittable for BVHNode {
  fn check_ray_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
    if !self.bounding_box.check_ray_hit(ray, t_min, t_max) {
      return None;
    }

    if let Some(hit) = self.left_child.check_ray_hit(ray, t_min, t_max) {
      return Some(hit);
    }

    if let Some(hit) = self.right_child.check_ray_hit(ray, t_min, t_max) {
      return Some(hit);
    }

    None
  }

  fn get_bounding_box(&self, _: f64, _: f64) -> Option<AxisAlignedBoundingBox> {
    Some(self.bounding_box)
  }
}

// TODO: This is dumb, but because I can't convert Vec<&Renderable> -> Vec<&Hittable>,
// I have to make this into a "colorable" object
impl Colorable for BVHNode {
  fn material(&self) -> &Material {
    &self.dummy_material
  }
}

impl Renderable for BVHNode {
  fn node_type(&self) -> NodeType {
    NodeType::Node(self)
  }
}

impl std::fmt::Debug for BVHNode {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
    f.write_str(&format!("[\n{:?}\n]", self.bounding_box))
  }
}
