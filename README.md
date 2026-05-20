# Rusty Ray Tracer

A simple path tracer written in Rust that supports loading `.obj` files and rendering scenes into `.exr` image files.

## Features

- **Path tracing with Next Event Estimation (NEE):** Stratified sampling splitting direct and indirect illumination for faster convergence.
- **OBJ file loading:** Import triangulated `.obj` scenes with paired `.mtl` material files.
- **Bounding Volume Hierarchy (BVH):** Acceleration structure built using the Surface Area Heuristic (SAH) for efficient ray-scene intersection.
- **Cosine-weighted importance sampling:** Hemisphere sampling weighted by cosine for indirect illumination.
- **Adaptive path termination (Russian Roulette):** Configurable fixed or probabilistic path length control.
- **Tile-based parallel rendering:** Multi-threaded rendering using Rayon with configurable tile size and channel bound.
- **EXR output:** High dynamic range output in OpenEXR format.
- **TOML configuration:** All camera and rendering parameters are defined in an external config file.
- **Configurable logging:** Adjustable log verbosity via CLI.

## Requirements

- Rust toolchain (stable, edition 2021+)
- A triangulated `.obj` file with an accompanying `.mtl` file

### Material Conventions

The `.mtl` file must contain at least one material whose name includes `light` or `Light`. This material will be treated as an emissive light source. The ambient color (`Ka`) of that material defines the light emission, and the diffuse color (`Kd`) defines its reflectance.

All other materials are treated as diffuse (Lambertian) surfaces.

## Installation

```bash
git clone <repository-url>
cd rusty-ray-tracer
cargo build --release
```

## Usage

```bash
cargo run --release -- [OPTIONS]
```

### Command-Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--object` | `-o` | Path to input `.obj` file (only triangulated meshes supported; paired `.mtl` must contain a material named `light` or `Light`) | `./cornell-box.obj` |
| `--no-graphical-coordinates` | `-g` | Assume mathematical axis orientation for coordinates in .obj file (vertical z-axis) | `false` |
| `--config` | `-c` | Path to TOML rendering configuration file | `./rendering.toml` |
| `--format` | `-f` | Output image format (`exr` or `png`) | `exr` |
| `--output-name` | `-n` | Output file name (without extension) | `rendering` |
| `--log` | `-l` | Log level (`TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, `OFF`) | `OFF` |

### Example

```bash
cargo run --release -- -o ./scenes/cornell-box.obj -c ./scene.toml -l INFO
```

This renders the scene and outputs `rendering.exr` in the working directory.

## Configuration

All camera and rendering parameters are defined in a TOML configuration file (mathematical axis orientation).

### Example `scene.toml`

```toml
[camera]
position = [0.0, 10.0, 3.7]
view_direction = [0.0, -10.0, -1.0]
focal_distance = 2.0
window_size = [2.0, 1.125]
from_graphical_coordinates = true

[render]
width = 2560
height = 1440
rays_per_pixel = 100
tile_size = 110
channel_bound = 12
backside_intersections = "ignore"

[render.path_length]
mode = "adaptive"
min_bounces = 1
termination_probability = 0.2
```

### Camera Parameters

| Parameter | Description |
|-----------|-------------|
| `position` | Camera position in 3D space `[x, y, z]` |
| `view_direction` | Direction the camera looks towards `[x, y, z]` |
| `focal_distance` | Distance from camera point to the view window plane |
| `window_size` | Width and height of the view window `[w, h]` |


### Render Parameters

| Parameter | Description |
|-----------|-------------|
| `width` | Horizontal resolution in pixels |
| `height` | Vertical resolution in pixels |
| `rays_per_pixel` | Number of samples per pixel |
| `tile_size` | Tile edge length for parallel work distribution |
| `channel_bound` | Bounded channel capacity (number of CPU cores is a good choice) |
| `backside_intersections` | `"ignore"` or `"end_ray"` |

### Path Length

The `[render.path_length]` section controls ray termination. Two modes are available:

**Fixed** - rays are terminated after a set number of bounces:

```toml
[render.path_length]
mode = "fixed"
max_bounces = 5
```

**Adaptive (Russian Roulette)** - rays are terminated probabilistically after a minimum number of bounces:

```toml
[render.path_length]
mode = "adaptive"
min_bounces = 1
termination_probability = 0.2
```

## Project Structure

```
.
├── src/
│   ├── main.rs              # Entry point, CLI, scene setup, rendering orchestration
│   └── rtcore/
│       ├── mod.rs           # Core: Ray, Scene, ColorStack, ViewPoint, path tracing logic
│       ├── objects.rs       # Triangle, Material, BVH (BBVTNode), Intersection
│       └── utilities.rs     # Random number generation helpers
├── scene.toml               # Scene and render configuration
├── Cargo.toml
└── README.md
```

## Dependencies

| Crate | Purpose |
|-------|---------|
| `clap` | Command-line argument parsing |
| `serde` | Serialization/deserialization for config |
| `toml` | TOML config file parsing |
| `nalgebra` | Linear algebra (vectors, math) |
| `tobj` | OBJ/MTL file loading |
| `exr` | OpenEXR image writing |
| `image` | PNG image writing (optional) |
| `rayon` | Data-parallel tile rendering |
| `crossbeam` | Bounded channels for producer/consumer pattern |
| `rand` | Random number generation |
| `tracing` / `tracing-subscriber` | Structured logging |
| `indicatif` | Progress bar |

## How It Works

1. **Configuration loading:** The TOML config file is parsed to obtain camera and render parameters.
2. **Scene loading:** The `.obj` file is parsed into triangles; materials are classified as diffuse or emissive.
3. **BVH construction:** All triangles are organized into an axis-aligned bounding box tree using SAH splits.
4. **Rendering:** The image is divided into tiles. Each tile is processed in parallel — for every pixel, `rays_per_pixel` rays are cast through random sub-pixel positions.
5. **Path tracing:** Each primary ray is intersected with the BVH. At each hit point:
   - **Direct lighting** is estimated by sampling a random point on a random light source and testing visibility (shadow ray).
   - **Indirect lighting** is estimated by spawning a cosine-weighted random ray into the hemisphere and recursing.
   - Paths are terminated by Russian Roulette or a fixed bounce limit.
6. **Output:** The accumulated radiance per pixel is written to a `rendering.exr` file.

## Limitations

- Only triangulated meshes are supported.
- Specular and glossy materials are defined but not yet fully implemented (fallback to diffuse BRDF).
- Normal interpolation (`face_normal` for non-flat shading) has a known issue and falls back to flat normals when vertex normals differ.

<!-- ## License

*Add your license here.* -->
