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

## Gallery
<img width="400" height="225" alt="balls" src="https://github.com/user-attachments/assets/aca5c01c-043e-4ee7-b949-deebefa1ad6e" />
<img width="400" height="225" alt="trippy" src="https://github.com/user-attachments/assets/37a3a8c7-243a-4370-a169-ff463fed362c" />


## Building and Running
As long as you have Rust installed w/ the `cargo` command-line tool, simply run

```
cargo run
```

from the base directory.

## Sources
Some functionality adapted from [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html) by Peter Shirley.
