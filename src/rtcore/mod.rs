//! Core of Rusty Ray Tracer
//!
//! All the most important structs, methods and functions that Rusty Ray Tracer builds on are defined here
//!
//!
//!
use exr::prelude::WritableImage;
use tracing::{debug, error, trace}; // debug, error, info, span, trace, warn};
use indicatif::{MultiProgress, ProgressBar, HumanDuration};
use indicatif::style::ProgressStyle;

use nalgebra::{Vector3, Matrix3};
use nalgebra::Rotation3;

use std::sync::mpsc;
use std::sync::Arc;
use std::thread;

use std::time::Instant;
use std::ops::{Add, AddAssign};
use core::panic;
use std::f64::consts;


pub mod utilities;
pub mod objectslib;


#[derive(Debug, Clone, Copy)]
pub enum RayColor {
    Red,
    Green,
    Blue,
}

impl RayColor {
    fn new_rndm() -> RayColor {
        match utilities::random3() {
            0 => { RayColor::Red },
            1 => { RayColor::Green },
            2 => { RayColor::Blue },
            _ => {
                error!("Invalid random value during RayColor init");
                panic!("Invalid random value during RayColor init");
            },
        }
    }
}

/// Struct representation of a Ray that can be cast into a [[Scene]]
#[derive(Debug)]
pub struct Ray {
    /// Base of Ray vector
    base: Vector3<f64>, // Redundant info if ray belongs to scene
    /// Normalized direction of Ray
    direction: Vector3<f64>, // scaled with scaling factor t
    color: RayColor,
}

impl Ray {
    /// Construct new Ray
    fn new(base: Vector3<f64>, direction: Vector3<f64>, color: RayColor) -> Self {
        Self {
            base,
            direction: direction/direction.norm(),
            color,
        }
    }
    /// Construct new Ray from base into a random direction defined by dir
    pub fn new_into_hemisphere(base: Vector3<f64>, dir: Vector3<f64>, color: RayColor) -> Ray {
        let azimuthal = utilities::random()*2.0*consts::PI;
        // polar angle is not equally distributed over its range (distribution: sin(theta))
        // according to the inversion method the function that needs to be applied to uniform samples x between 0 and 1
        // so that they adhere to a sine distribution is arccos(1-x)
        let polar = (1.0-utilities::random()).acos();
        let mut rand_dir_in_hemisphere = match dir.dot(&Vector3::z_axis()) {
            1.0 | -1.0 => {
                let perp = Vector3::new(1.0, 0.0, 0.0);
                let rot_perp = Rotation3::new(polar*perp);
                let rot_para = Rotation3::new(azimuthal*dir/dir.norm());
                rot_para*(rot_perp*dir)
            },
            _ => {
                let perp = dir.cross(&Vector3::z_axis());
                let rot_perp = Rotation3::new(polar*perp/perp.norm());
                let rot_para = Rotation3::new(azimuthal*dir/dir.norm());
                rot_para*(rot_perp*dir)
            },
        };
        rand_dir_in_hemisphere = rand_dir_in_hemisphere/rand_dir_in_hemisphere.norm();
        if dir.dot(&rand_dir_in_hemisphere) < 0.0 {
            error!("Created ray is not in hemisphere");
            todo!("Created ray is not in hemisphere");
        }
        trace!("{}", rand_dir_in_hemisphere.norm());
        Ray {
            base,
            direction: rand_dir_in_hemisphere, // /rand_dir_in_hemisphere.norm(),
            color,
        }
    }
    fn intersect_triangle(&self, p0: Vector3<f64>, p1: Vector3<f64>, p2: Vector3<f64>) -> Option<(f64, f64, f64)> {
        let s = self.base-p0;
        let e1 = p1-p0;
        let e2 = p2-p0;
        let factor = 1.0/self.direction.cross(&e2).dot(&e1);
        let t = factor*s.cross(&e1).dot(&e2);
        let b1 = factor*self.direction.cross(&e2).dot(&s);
        let b2 = factor*s.cross(&e1).dot(&self.direction);
        if b1 >= 0.0 && b2 >= 0.0 && b1+b2 <= 1.0 && t > 0.0 {
            // intersection
            return Some((t, b1, b2));
        }
        None
    }
    /// Check if a ray is absorbed depending on its color
    pub fn is_absorbed(&self, refl: &objectslib::Reflectance) -> bool {
        let random_val =  utilities::random();
        match self.color {
            RayColor::Red => { random_val >= refl.r },
            RayColor::Green => { random_val >= refl.g },
            RayColor::Blue => { random_val >= refl.b },
        }
    }
}


pub struct ColorStack {
    accumulated_color: Vec<Vec<Vector3<f64>>>,
    // height: Vec<Vec<Vector3<f64>>>,
}

impl ColorStack {
    fn new(dim_x: u32, dim_y: u32) -> Self {
        ColorStack {
            accumulated_color: vec![vec![Vector3::zeros(); dim_x as usize]; dim_y as usize],
            // height: vec![vec![Vector3::zeros(); dim_x as usize]; dim_y as usize],
        }
    }
    // fn inc_height(&mut self, x: usize, y: usize, color: RayColor) {
    //     match color {
    //         RayColor::Red => self.height[y][x][0] += 1.0,
    //         RayColor::Green => self.height[y][x][1] += 1.0,
    //         RayColor::Blue => self.height[y][x][2] += 1.0,
    //     }
    // }
}

impl Add for ColorStack {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.accumulated_color.is_empty() || other.accumulated_color.is_empty() {
            panic!("Cannot add empty instances of ColorStack")
        }
        if self.accumulated_color.len() != other.accumulated_color.len()
        || self.accumulated_color[0].len() != other.accumulated_color[0].len() {
            panic!("Cannot add two instances of ColorStack with different dimensions")
        }
        let mut result = Self::new(self.accumulated_color[0].len() as u32, self.accumulated_color.len() as u32);
        for i in 0..self.accumulated_color.len() {
            for j in 0..self.accumulated_color[0].len() {
                result.accumulated_color[i][j] = self.accumulated_color[i][j] + other.accumulated_color[i][j];
                // result.height[i][j] = self.height[i][j] + other.height[i][j];
            }
        }
        result
    }
}

impl AddAssign for ColorStack {
    fn add_assign(&mut self, other: Self) {
        if self.accumulated_color.is_empty() || other.accumulated_color.is_empty() {
            panic!("Cannot add empty instances of ColorStack")
        }
        if self.accumulated_color.len() != other.accumulated_color.len()
        || self.accumulated_color[0].len() != other.accumulated_color[0].len() {
            panic!("Cannot add two instances of ColorStack with different dimensions")
        }
        let mut result = Self::new(self.accumulated_color[0].len() as u32, self.accumulated_color.len() as u32);
        for i in 0..self.accumulated_color.len() {
            for j in 0..self.accumulated_color[0].len() {
                result.accumulated_color[i][j] = self.accumulated_color[i][j] + other.accumulated_color[i][j];
                // result.height[i][j] = self.height[i][j] + other.height[i][j];
            }
        }
        *self = result;
    }
}

impl std::convert::Into<image::ImageBuffer<image::Rgb::<u8>, Vec<u8>>> for ColorStack {
    fn into(self) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>> {
        if !self.accumulated_color.is_empty() {
            // Create a new image buffer with resolution: res_x * res_y
            let mut imgbuf = image::ImageBuffer::new(
                self.accumulated_color[0].len().try_into().unwrap(),
                self.accumulated_color.len().try_into().unwrap());

            // Populate the image bugger with the
            for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
                // let height_x_y_r = color_stack.height[y as usize][x as usize][0];
                // let height_x_y_g = color_stack.height[y as usize][x as usize][1];
                // let height_x_y_b = color_stack.height[y as usize][x as usize][2];
                // let r = (color_stack.accumulated_color[y as usize][x as usize][0]/height_x_y_r) as u8;
                // let g = (color_stack.accumulated_color[y as usize][x as usize][1]/height_x_y_g) as u8;
                // let b = (color_stack.accumulated_color[y as usize][x as usize][2]/height_x_y_b) as u8;
                // map colors
                let scale = 30.0;
                let r = self.accumulated_color[y as usize][x as usize][0]*scale;
                let r = if r <= 255.0 {r as u8} else {255};
                let g = self.accumulated_color[y as usize][x as usize][1]*scale;
                let g = if g <= 255.0 {g as u8} else {255};
                let b = self.accumulated_color[y as usize][x as usize][2]*scale;
                let b = if b <= 255.0 {b as u8} else {255};
                *pixel = image::Rgb([r, g, b]);
            }
            imgbuf
        } else {
            image::ImageBuffer::default()
        }
    }
}

// std::convert::Into<exr::prelude::Image<exr::prelude::Layer<exr::prelude::SpecificChannels<(f32, f32, f32), (exr::prelude::ChannelDescription, exr::prelude::ChannelDescription, exr::prelude::ChannelDescription)>>>> for
impl ColorStack {
    pub fn into_exr(self) -> exr::prelude::Image<exr::prelude::Layer<exr::prelude::SpecificChannels<impl Fn(exr::prelude::Vec2<usize>) -> (f32, f32, f32), (exr::prelude::ChannelDescription, exr::prelude::ChannelDescription, exr::prelude::ChannelDescription)>>>
    {
        // Create a new image with resolution: res_x * res_y
        let colorstack = self.accumulated_color.clone();
        let channels = exr::prelude::SpecificChannels::rgb(move |pos: exr::prelude::Vec2<usize>| {
            (colorstack[pos.1][pos.0][0] as f32,
                colorstack[pos.1][pos.0][1] as f32,
                colorstack[pos.1][pos.0][2] as f32)
        });
        let image = exr::prelude::Image::from_layer(
            exr::prelude::Layer::new(
                (self.accumulated_color[0].len(), self.accumulated_color.len()),
                exr::prelude::LayerAttributes::named("main-rgb-layer"),
                exr::prelude::Encoding::UNCOMPRESSED,
                channels
            )
        );
        image
    }
}


#[derive(Debug, Clone)]
pub struct ViewWindow {
    pub direction: Vector3<f64>, // surface normal in kartesian coodinates
    pub size: (f64, f64),
}

impl ViewWindow {
    pub fn from(direction: Vector3<f64>, size: (f64, f64)) -> Self {
        Self {
            direction: direction/direction.norm(),
            size,
        }
    }
    fn get_vector_from_center(&self, offset_horizontal: f64, offset_vertical: f64) -> Vector3<f64> {
        let horizontal = self.direction.cross(&Vector3::z_axis());
        let vertical = horizontal.cross(&self.direction);

        horizontal/horizontal.norm()*offset_horizontal + vertical/vertical.norm()*offset_vertical
    }
}

#[derive(Debug, Clone)]
pub struct ViewPoint {
    pub point: Vector3<f64>,
    pub window: ViewWindow,
    pub distance_to_window_plane: f64,
}

impl ViewPoint {
    /// Create ray from view_point onto window with random, uniformly distributed position of intersection (with window)
    ///
    /// Output: Tuple containing:
    /// Ray,
    /// x position on ViewWindow (left is 0, right is width of window),
    /// y position on ViewWindow (bottom is 0, top is height of window)
    fn create_rndm_ray_through_pixel(&self, pixel_x: usize, pixel_y: usize, pixel_width: f64, pixel_height: f64) -> Ray {
        let rndm_pos = utilities::random_2d_pos(
            0.0..=pixel_width,
            0.0..=pixel_height);
        // calculate pixel position (left-upper-corner)
        let pixel_offset_from_center_x = -self.window.size.0/2.0+(pixel_x as f64)*pixel_width;
        let pixel_offset_from_center_y = self.window.size.1/2.0-(pixel_y as f64)*pixel_height;
        let vec_center_window_to_rndm_pos = self.window.get_vector_from_center(
            pixel_offset_from_center_x+rndm_pos.0,
            pixel_offset_from_center_y-rndm_pos.1);

        let direction_from_point_to_center_of_window = self.window.direction
            *self.distance_to_window_plane;

        let direction_from_point_to_rndm_pos_on_window =
            direction_from_point_to_center_of_window + vec_center_window_to_rndm_pos;

        Ray::new(self.point, direction_from_point_to_rndm_pos_on_window, RayColor::new_rndm())
    }
}

// implement polar coordinates (azimuthal angle, polar angle)?
// #[derive(Debug)]
// enum CoordinateSystem {
//     Kartesian,
//     Cylindric,
//     Spherical,
// }

#[derive(Debug, Clone)]
pub struct Scene {
    pub view_point: ViewPoint,
    // light_sources: Vec<i8>,
    pub models: Vec<tobj::Model>,
    // coordinate_system: CoordinateSystem,
}

impl Scene {
    fn eval_intersection(
        &self,
        ray: &Ray,
        model: &tobj::Model,
        _triangle_indices: (usize, usize, usize),
        normal_index: usize,
        t_b1_b2: (f64, f64, f64)
    ) -> Vector3<f64> {
        // handle absorption randomly
        // cast new ray
        if model.name == "Plane1" {
            let reflectance = objectslib::Reflectance::new(0.3, 0.3, 0.3);
            // let face_normal = Vector3::from_vec(model.mesh.normals[0..=2]
            //     .iter()
            //     .map(|val| { *val as f64 } )
            //     .collect());
            // let ray_refl = Ray {
            //     // intersection point
            //     base: ray.base+ray.direction*t_b1_b2.0,
            //     // reflected ray direction (omega_in = -ray:direction)
            //     direction: -2.0*ray.direction.dot(&face_normal)*face_normal+ray.direction,
            //     color: ray.color,
            // };
            // match self.trace_ray(&ray_refl){
            //     Some(res) => {
            //         res/2
            //     },
            //     None => {
            //         Vector3::from([100, 100, 100])
            //     }
            // }
            // check for absorption
            if ray.is_absorbed(&reflectance) { // absorbed
                Vector3::new(0.0, 0.0, 0.0)
            } else { // not absorbed
                let face_normal = Vector3::from_vec(model.mesh.normals[0..=2]
                    .iter()
                    .map(|val| { *val as f64 } )
                    .collect());
                let ray_refl =  Ray::new_into_hemisphere(
                    ray.base+ray.direction*t_b1_b2.0,
                    face_normal,
                    ray.color);
                match self.trace_ray(&ray_refl){
                    Some(res) => {
                        res/consts::PI
                    },
                    None => {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                }
            }
        } else if model.name == "Plane2" || model.name == "Plane3" {
            let reflectance = objectslib::Reflectance::new(0.3, 0.3, 0.3);
            // check for absorption
            if ray.is_absorbed(&reflectance) { // absorbed
                Vector3::new(0.0, 0.0, 0.0)
            } else { // not absorbed
                let face_normal = Vector3::from_vec(model.mesh.normals[0..=2]
                    .iter()
                    .map(|val| { *val as f64 } )
                    .collect());
                let ray_refl =  Ray::new_into_hemisphere(
                    ray.base+ray.direction*t_b1_b2.0,
                    face_normal,
                    ray.color);
                match self.trace_ray(&ray_refl){
                    Some(res) => {
                        res/consts::PI
                    },
                    None => {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                }
            }
            // Vector3::from([100, 100, 150])
        } else if model.name == "Ls" {
            // let reflectance = Vector3::new(0.3, 0.3, 0.3);
            match ray.color {
                RayColor::Red => Vector3::from([1.0, 0.0, 0.0]),
                RayColor::Green => Vector3::from([0.0, 1.0, 0.0]),
                RayColor::Blue => Vector3::from([0.0, 0.0, 0.9]),
            }
        } else if model.name == "Cube" {
            let reflectance = objectslib::Reflectance::new(0.9, 0.1, 0.1);
            // check for absorption
            if ray.is_absorbed(&reflectance) { // absorbed
                Vector3::new(0.0, 0.0, 0.0)
            } else { // not absorbed
                let face_normal = Vector3::from_vec(model.mesh.normals[3*normal_index..=(3*normal_index+2)]
                    .iter()
                    .map(|val| { *val as f64 } )
                    .collect());
                let ray_refl =  Ray::new_into_hemisphere(
                    ray.base+ray.direction*t_b1_b2.0,
                    face_normal,
                    ray.color);
                trace!("{}, {}, {:?}, {}", normal_index, face_normal, ray_refl, face_normal.dot(&ray_refl.direction));
                match self.trace_ray(&ray_refl){
                    Some(res) => {
                        res/consts::PI
                    },
                    None => {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                }
            }
            // Vector3::from([255, 0, 0])
        } else { // pro forma
            error!("Could not evaluate intersection. Hit object not found");
            Vector3::from([0.0, 0.0, 0.0])
        }
    }
    fn trace_ray(&self, ray: &Ray) -> Option<Vector3<f64>> {
        // determine first intersection with an object
        // init first intersection
        let mut first_model: Option<&tobj::Model> = None;
        let mut triangle_indices: Option<(usize, usize, usize)> = None;
        let mut normal_index: Option<usize> = None;
        let mut t_b1_b2: Option<(f64, f64, f64)> = None;
        // go through all models
        for model in &self.models {
            if !&model.mesh.face_arities.is_empty() {
                panic!("Found something else than a triangle");
            }

            // go through all triangles (only triangles are supported)
            for (i, ind_tri) in model.mesh.indices.chunks(3).enumerate() {
                // extract vertices
                let p0 = Vector3::from([
                    (model.mesh.positions[3*ind_tri[0] as usize]) as f64,
                    (model.mesh.positions[(3*ind_tri[0] as usize)+1]) as f64,
                    (model.mesh.positions[(3*ind_tri[0] as usize)+2]) as f64]);
                let p1 = Vector3::from([
                    (model.mesh.positions[3*ind_tri[1] as usize]) as f64,
                    (model.mesh.positions[(3*ind_tri[1] as usize)+1]) as f64,
                    (model.mesh.positions[(3*ind_tri[1] as usize)+2]) as f64]);
                let p2 = Vector3::from([
                    (model.mesh.positions[3*ind_tri[2] as usize]) as f64,
                    (model.mesh.positions[(3*ind_tri[2] as usize)+1]) as f64,
                    (model.mesh.positions[(3*ind_tri[2] as usize)+2]) as f64]);

                // intersect ray with current triangle
                if let Some((t, b1, b2)) = ray.intersect_triangle(p0, p1, p2) {
                    match t_b1_b2 {
                        Some(prev) => {
                            if t < prev.0 {
                                t_b1_b2 = Some((t, b1, b2));
                                first_model = Some(model);
                                triangle_indices = Some((
                                    ind_tri[0] as usize,
                                    ind_tri[1] as usize,
                                    ind_tri[2] as usize));
                                normal_index = Some(model.mesh.normal_indices[3*i] as usize);
                            }
                        },
                        None => {
                            t_b1_b2 = Some((t, b1, b2));
                            first_model = Some(model);
                            triangle_indices = Some((
                                ind_tri[0] as usize,
                                ind_tri[1] as usize,
                                ind_tri[2] as usize));
                            normal_index = Some(model.mesh.normal_indices[3*i] as usize);
                        },
                    }
                }
            }
        }
        // evaluate first intersection
        if let Some(fm) = first_model {
            return Some(self.eval_intersection(
                ray,
                fm,
                triangle_indices.unwrap(),
                normal_index.unwrap(),
                t_b1_b2.unwrap()));
        }
        None
    }
    fn trace_rays_in_new_thread(
        transmitter: mpsc::Sender<ColorStack>,
        scene: Arc<Self>,
        // scene: Self,
        // window_size: (f64, f64),
        res_x: u32,
        res_y: u32,
        rays_numbered: std::ops::Range<u64>,
        bar: Arc<ProgressBar>,
    ) -> thread::JoinHandle<()> {
        thread::spawn(move || {
            let mut color_stack = ColorStack::new(res_x, res_y);
            // total number of pixels
            let n_pixels = (res_x as u64)*(res_y as u64);
            // calc width and height of a pixel
            let pixel_width = scene.view_point.window.size.0/(res_x as f64);
            let pixel_height = scene.view_point.window.size.1/(res_y as f64);
            for i in rays_numbered {
                // determine associated pixel (start top left in lines to bottom right)
                let pixel_x = (i%res_x as u64) as usize;
                let pixel_y = ((i/res_x as u64)%(res_y as u64)) as usize;

                let ray = scene.view_point
                    .create_rndm_ray_through_pixel(pixel_x, pixel_y, pixel_width, pixel_height);

                // intersect ray with scene objects and eval intersections
                if let Some(ray_data) = scene.trace_ray(&ray) {
                    color_stack.accumulated_color[pixel_y][pixel_x] += ray_data;
                }
                // color_stack.inc_height(pixel_x, pixel_y, ray.color);

                if i%n_pixels == n_pixels-1 {
                    bar.inc(1);
                }
            }
            bar.abandon();
            // send result
            transmitter.send(color_stack).unwrap();
        })
    }
    pub fn look(&self, res_x: u32, res_y: u32, n_rays_per_pixel: u64, n_threads: u8) -> ColorStack {
        // set starting time
        // let started = Instant::now();

        // create pointer to scene
        // Note: Copying scene for independent computing does not seem worth at all
        // (no computation time improvement at this stage)
        let scene = Arc::new(self.clone());

        // init progress bar (for optimum performance remove bar)
        let multi_p_bar = MultiProgress::new();
        let bar_style = ProgressStyle::with_template(
                "[{elapsed_precise}]{wide_bar:.cyan/blue} {pos}/{len} rays [{eta_precise}]\n{msg}").unwrap();
        let bar = Arc::new(multi_p_bar.add(ProgressBar::new(n_rays_per_pixel)));
        bar.set_style(bar_style.clone());

        // trace rays in thread
        let receiver = {
            let (transmitter, receiver) = mpsc::channel();
            // Number of rays
            let n_rays_total: u64 = (res_x as u64)*(res_y as u64)*n_rays_per_pixel;
            let n_rays_per_thread = n_rays_total/(n_threads as u64);
            let n_rays_per_thread_rest = n_rays_total%(n_threads as u64);
            for i in 0..n_threads {
                // calculate number of rays computed by this thread
                let rays_numbered = if (i as u64) < n_rays_per_thread_rest {
                    ((n_rays_per_thread+1)*(i as u64))..((n_rays_per_thread+1)*((i+1) as u64))
                } else {
                    (n_rays_per_thread*(i as u64)+n_rays_per_thread_rest)..(n_rays_per_thread*((i+1) as u64)
                        +n_rays_per_thread_rest)
                };

                let transmitter = transmitter.clone();
                let scene = Arc::clone(&scene);
                let bar = Arc::clone(&bar);

                Self::trace_rays_in_new_thread(
                    transmitter,
                    scene,
                    res_x,
                    res_y,
                    rays_numbered,
                    bar);
            }
            receiver
        };

        // initialize accumulated color matrix that will later be used to populate the image buffer
        let mut color_stack = ColorStack::new(res_x, res_y);

        // receive messages until all transmitters are dropped (and processes are finished)
        for rcv_msg in receiver.iter() {
            color_stack += rcv_msg;
        }

        // multi_p_bar..clear().unwrap();
        // println!("Finished evaluating {} rays in {}", n_rays, HumanDuration(started.elapsed()));
        println!("Finished... ");

        color_stack
    }
}
