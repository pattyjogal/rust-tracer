# A Basic Ray-Tracer (Ray-Caster?) in Rust

A simple ray-tracing renderer written in Rust.

Currently supports:
- Objects
  - Triangles
  - Spheres
  - Planes
  - Meshes
- Light Effects
  - Hard Shadows
  - Diffuse Shading
  - Specular Shading
  - Smooth Mesh Shading
- Anti-Aliasing
  - Multi-Jittered Sampling
- Camera Projections
  - Perspective
  - Orthographic

## Building and Running
As long as you have Rust installed w/ the `cargo` command-line tool, simply run

```
cargo run
```

from the base directory.

## Sources
Some functionality adapted from [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html) by Peter Shirley.