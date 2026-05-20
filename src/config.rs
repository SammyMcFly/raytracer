use serde::Deserialize;

#[derive(Deserialize, Debug)]
pub struct Config {
    pub camera: CameraConfig,
    pub render: RenderConfig,
}

#[derive(Deserialize, Debug)]
pub struct CameraConfig {
    pub position: [f64; 3],
    pub view_direction: [f64; 3],
    pub focal_distance: f64,
    pub window_size: [f64; 2],
}

#[derive(Deserialize, Debug)]
pub struct RenderConfig {
    pub width: usize,
    pub height: usize,
    pub rays_per_pixel: u64,
    pub tile_size: usize,
    pub channel_bound: usize,
    pub backside_intersections: String,       // "ignore" or "end_ray"
    pub path_length: PathLengthConfig,
}

#[derive(Deserialize, Debug, Clone)]
#[serde(tag = "mode")]
pub enum PathLengthConfig {
    #[serde(rename = "fixed")]
    Fixed { max_bounces: u64 },
    #[serde(rename = "adaptive")]
    Adaptive { min_bounces: u64, termination_probability: f64 },
}
