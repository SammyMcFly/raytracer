//! Core of Rusty Ray Tracer
//!
//! All the most important structs, methods and functions that Rusty Ray Tracer builds on are defined here
//!
//! This ray tracer builds on the concept of path tracing with next event estimation:
//! - stratification of reflected ray domain into light sources and indirect lighting (whole hemisphere)
//! -
//!
use objects::{BBVTNode, Intersection, Reflectance, Triangle};
use tracing::{debug, error, info, trace}; // debug, error, info, span, trace, warn};
use indicatif::{MultiProgress, ProgressBar, HumanDuration};
use indicatif::style::ProgressStyle;

use nalgebra::{Vector3, Matrix3};
use nalgebra::Rotation3;

use std::mem::swap;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;

use std::time::Instant;
use std::ops::{Add, AddAssign};
use core::panic;
use std::f64::consts;


pub mod utilities;
pub mod objects;




type Color = Vector3<f64>;

pub enum SampleWeighting {
    Uniform,
    Cosine,
}

/// Struct representation of a Ray that can be cast into a [[Scene]]
#[derive(Debug)]
pub struct Ray {
    /// Base of Ray vector
    base: Vector3<f64>, // Redundant info if ray belongs to scene
    /// Normalized direction of Ray
    direction: Vector3<f64>, // scaled with scaling factor t
    // color: RayColor,
    n_bounces: u64,
}

impl Ray {
    /// Construct new Ray
    fn new(base: Vector3<f64>, direction: Vector3<f64>, n_bounces: u64) -> Self {
        Self {
            base,
            direction: direction/direction.norm(),
            // color,
            n_bounces,
        }
    }
    /// Construct new Ray from base into hemisphere defined by dir
    ///
    /// Direction is random and weighted according to "weighting"
    pub fn new_into_hemisphere(base: Vector3<f64>, dir: Vector3<f64>,
            n_bounces: u64, weighting: SampleWeighting) -> Ray {
        info!("Creating new ray with random direction inside the hemisphere");
        let azimuthal = utilities::random()*2.0*consts::PI;
        // polar angle is not equally distributed over its range (distribution: sin(theta))
        // according to the inversion method the function that needs to be applied to uniform samples x between 0 and 1
        // so that they adhere to a sine distribution is arccos(1-x)
        let polar = match weighting {
            SampleWeighting::Uniform => (1.0-utilities::random()).acos(),
            SampleWeighting::Cosine => (utilities::random().sqrt()).asin(),
        };
        trace!("Angles of ray creation: azimuthal = {azimuthal}, polar = {polar}");
        trace!("{}", dir);
        let mut rand_dir_in_hemisphere = match dir.dot(&Vector3::new(0.0, 0.0, 1.0)) {
            1.0 | -1.0 => {
                let perp = Vector3::new(1.0, 0.0, 0.0);
                let rot_perp = Rotation3::new(polar*perp);
                let rot_para = Rotation3::new(azimuthal*dir/dir.norm());
                trace!("{}", rot_para*(rot_perp*dir));
                rot_para*(rot_perp*dir)
            },
            _ => {
                let perp = dir.cross(&Vector3::z_axis());
                let rot_perp = Rotation3::new(polar*perp/perp.norm());
                let rot_para = Rotation3::new(azimuthal*dir/dir.norm());
                trace!("{}", rot_para*(rot_perp*dir));
                rot_para*(rot_perp*dir)
            },
        };
        rand_dir_in_hemisphere = rand_dir_in_hemisphere/rand_dir_in_hemisphere.norm();
        if dir.dot(&rand_dir_in_hemisphere) < 0.0 {
            error!("Created ray is not in hemisphere");
            panic!("Created ray is not in hemisphere");
        }
        trace!("Length of vector {}", rand_dir_in_hemisphere.norm());
        Self::new(
            base,
            rand_dir_in_hemisphere,
            n_bounces,
        )
    }
    /// Intersect ray with triangle
    ///
    /// Return:
    /// - scaling of ray t to intersection (>0)
    /// - scaling of vector between p0 and p1 (to intersection)
    /// - scaling of vector between p0 and p2 (to intersection)
    fn intersect_triangle(&self, p0: &Vector3<f64>, p1: &Vector3<f64>, p2: &Vector3<f64>) -> Option<(f64, f64, f64)> {
        info!("Intersecting ray with triangle");
        let s = self.base-p0;
        let e1 = p1-p0;
        let e2 = p2-p0;
        let factor = 1.0/self.direction.cross(&e2).dot(&e1);
        let t = factor*s.cross(&e1).dot(&e2);
        let b1 = factor*self.direction.cross(&e2).dot(&s);
        let b2 = factor*s.cross(&e1).dot(&self.direction);
        if b1 >= 0.0 && b2 >= 0.0 && b1+b2 <= 1.0 && t > 0.0 {
            // intersection
            info!("Found intersection between ray and triangle");
            return Some((t, b1, b2));
        }
        info!("No intersection between ray and triangle found");
        None
    }
    /// Intersect ray with (cubic, axis-aligned) bounding volume
    ///
    /// See: https://www.pbr-book.org/3ed-2018/Shapes/Basic_Shape_Interface#RayndashBoundsIntersections
    ///
    /// Return:
    /// - scaling of ray t to intersection (>0)
    /// - scaling of vector between p0 and p1 (to intersection)
    /// - scaling of vector between p0 and p2 (to intersection)
    fn intersect_bv(&self, bounds_min: Vector3<f64>, bounds_max: Vector3<f64>) -> bool {
        let mut t_0: Option<f64> = None;
        let mut t_1: Option<f64> = None;
        for (i, (b_min_i, b_max_i)) in bounds_min.iter().zip(&bounds_max).enumerate() {
            // inverse ray direction
            let ray_dir_inv = 1.0/self.direction[i];
            let mut t_near = (b_min_i-self.base[i])*ray_dir_inv;
            let mut t_far = (b_max_i-self.base[i])*ray_dir_inv;
            if t_near > t_far {
                swap(&mut t_near, &mut t_far);
            }
            // update t_0
            if let Some(val) = t_0 {
                t_0 = if t_near > val { Some(t_near) } else { Some(val) };
            } else {
                t_0 = Some(t_near);
            }
            // update t_1
            if let Some(val) = t_1 {
                t_1 = if t_far < val { Some(t_far) } else { Some(val) };
            } else {
                t_1 = Some(t_far);
            }
            if t_0 > t_1 {
                return false;
            }
        }
        t_1.unwrap() > 0.0
    }
}

/// 2D-Array that stores a [[Color]] for every pixel
#[derive(Debug)]
pub struct ColorStack {
    accumulated_color: Vec<Vec<Color>>,
    // height: Vec<Vec<Vector3<f64>>>,
}

impl ColorStack {
    fn new(dim_x: u32, dim_y: u32) -> Self {
        ColorStack {
            accumulated_color: vec![vec![Vector3::zeros(); dim_x as usize]; dim_y as usize],
            // height: vec![vec![Vector3::zeros(); dim_x as usize]; dim_y as usize],
        }
    }
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
        info!("Writing ColorStack to file");
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
        info!("Create random ray through pixel with index ({pixel_x},{pixel_y})");
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

        Ray::new(self.point, direction_from_point_to_rndm_pos_on_window, 0)
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
    view_point: ViewPoint,
    // /// Collection of all 3D veritces in the scene represented by Vector3<f64>
    // vertices: Vec<Arc<Vector3<f64>>>,
    // /// Collection of all 3D normal vectors of surfaces represented by Vector3<f64>
    // normals: Vec<Arc<Vector3<f64>>>,
    // /// Collection of all surface materials
    // materials: Vec<Arc<objects::Material>>,
    /// Collection of all light sources
    light_sources: Vec<objects::Triangle>,
    // /// Collection of all triangle in the scene, they index 3 vertices, normals and materials
    // // triangles: Vec<objectslib::Triangle>,
    /// Collection of all light sources
    // light_sources: Vec<i8>,
    /// Binary bounding volume tree, which is searched for every intersection test
    ///
    /// Has all triangles of the scene as leafs
    bbv_tree_root: objects::BBVTNode,
    /// Maximum bounces a ray does in the scene until it is terminated
    max_bounces: u64,

    // coordinate_system: CoordinateSystem,
}

impl Scene {
    /// Create Scene
    ///
    /// models is expected to contain models consisting only of triangles and having a mdl.mesh.material_id
    pub fn from(
            view_point: ViewPoint,
            models: Vec<tobj::Model>,
            materials: Vec<tobj::Material>,
            max_bounces: u64,
    ) -> Self {
        info!("Creating scene");
        let mut materials_extracted = Vec::new();
        for mat in materials {
            let [r, g, b] = mat.diffuse.unwrap();
            if mat.name.contains("light_source") {
                let light_emission = Reflectance::from_vec(
                    Vec::from(mat.ambient.unwrap()).iter().map(|a| *a as f64).collect());
                materials_extracted.push(Arc::new(objects::Material::DiffuseLightSource {
                    refl: Reflectance::new(r as f64, g as f64, b as f64),
                    light_emission,
                }));
            } else {
                materials_extracted.push(Arc::new(objects::Material::Diffuse {
                    refl: Reflectance::new(r as f64, g as f64, b as f64),
                }));
            }
        }

        let mut vertices = Vec::new();
        let mut normals = Vec::new();
        let mut triangles = Vec::new();
        let mut light_sources = Vec::new();
        for mdl in models {
            let mut vertices_i = Vec::new();
            let mut normals_i = Vec::new();
            for vertex in mdl.mesh.positions.chunks(3) {
                vertices_i.push(Arc::new(Vector3::new(vertex[0] as f64, vertex[1] as f64, vertex[2] as f64)));
            }
            for normal in mdl.mesh.normals.chunks(3) {
                normals_i.push(Arc::new(Vector3::new(normal[0] as f64, normal[1] as f64, normal[2] as f64)));
            }
            for (vertex_indices, normal_indices) in mdl.mesh.indices.chunks(3)
            .zip(mdl.mesh.normal_indices.chunks(3)){
                let tr = objects::Triangle::new(
                    vertices_i[vertex_indices[0] as usize].clone(),
                    vertices_i[vertex_indices[1] as usize].clone(),
                    vertices_i[vertex_indices[2] as usize].clone(),
                    normals_i[normal_indices[0] as usize].clone(),
                    normals_i[normal_indices[1] as usize].clone(),
                    normals_i[normal_indices[2] as usize].clone(),
                    materials_extracted[mdl.mesh.material_id.unwrap()].clone(),
                );
                triangles.push(tr.clone());
                if materials_extracted[mdl.mesh.material_id.unwrap()].is_light_source() {
                    light_sources.push(tr);
                }

            }
            vertices.extend(vertices_i);
            normals.extend(normals_i);
        }

        // build bounding volume hierarchy
        let bbv_tree_root = BBVTNode::build_bbv_tree(triangles);

        debug!("Outer bounds: {:?}, {:?}", bbv_tree_root.bounds_min, bbv_tree_root.bounds_max);

        Self {
            view_point,
            // vertices,
            // normals,
            // materials: materials_extracted,
            light_sources,
            bbv_tree_root,
            max_bounces,
        }
    }
    /// Currently does not support normal vector interpolation
    fn get_random_light_source(&self) -> (&objects::Triangle, f64) { // , reference_point: Vector3<f64>, triangle: objects::Triangle
        // todo check that light source is not selected from the same light source
        // // select light source, which face the reference point and are located in the hemisphere defined by
        // // the reference point's normal vector
        // let mut ls_candidates: Vec<Triangle> = Vec::new();
        // for ls in &self.light_sources {
        //     let mut connection_1 = Vector3::new(0.0, 0.0, 0.0);
        //     ls.v1.sub_to(&reference_point, &mut connection_1);
        //     let mut connection_2 = Vector3::new(0.0, 0.0, 0.0);
        //     ls.v2.sub_to(&reference_point, &mut connection_2);
        //     let mut connection_3 = Vector3::new(0.0, 0.0, 0.0);
        //     ls.v3.sub_to(&reference_point, &mut connection_3);
        //     if (ls.n1.dot(&connection_1) < 0.0 && connection_1.dot(&triangle.face_normal(&reference_point)) > 0.0)
        //             || (ls.n2.dot(&connection_2) < 0.0 && connection_2.dot(&triangle.face_normal(&reference_point)) > 0.0)
        //             || (ls.n3.dot(&connection_3) < 0.0 && connection_3.dot(&triangle.face_normal(&reference_point)) > 0.0) {
        //         ls_candidates.push(ls.clone());
        //     }
        // }
        let rndm_ls_index = (utilities::random()*self.light_sources.len() as f64) as usize;
        let rndm_ls = &self.light_sources[rndm_ls_index];
        (rndm_ls, 1.0/self.light_sources.len() as f64)
    }
    /// Evaluate intersection according to path tracing with next event estimation
    ///
    /// Terminate rays after max_bounces bounces
    fn eval_intersection(
        &self,
        intersection: objects::Intersection,
        light_emission: bool,
    ) -> Vector3<f64> {
        info!("Evaluating ray-scene-intersection: {intersection:?}, consider light emission: {light_emission}");
        // init intersection point
        let point_of_intersection = intersection.ray.base
            +intersection.ray.direction*intersection.ray_scaling;
        // eval intersection regarding direct lighting
        let randiance_emitted = if light_emission && intersection.triangle.mat.is_light_source() {
            Some(intersection.triangle.mat.light_emission())
        } else {
            None
        };
        // check number of bounces
        let lighting = if intersection.ray.n_bounces > self.max_bounces { // terminate ray
            Vector3::new(0.0, 0.0, 0.0)
        } else { // continue tracing ray
            // get face normal
            let face_normal = intersection.triangle.face_normal(&point_of_intersection);
            trace!("Face normal: {face_normal}");
            // eval primary lighting
            let (ls, probability_ls) = self.get_random_light_source();
                // &point_of_intersection,
                // &intersection.triangle);
            let (rndm_point_on_ls, ray_to_ls, area) = ls
                .random_point_on_surface(&point_of_intersection, intersection.ray.n_bounces);
            let ls_intersection = Intersection {
                ray: &ray_to_ls,
                triangle: ls,
                ray_scaling: (rndm_point_on_ls-point_of_intersection).norm(),
            };
            let ls_face_normal = ls.face_normal(&rndm_point_on_ls);
            let randiance_direct = match self.trace_to_lightsource(&ray_to_ls, &ls_intersection){
                Some(randiance_direct)
                if ray_to_ls.direction.dot(&(face_normal/face_normal.norm())) > 0.0
                && ray_to_ls.direction.dot(&(ls_face_normal/ls_face_normal.norm())) < 0.0 => {
                    trace!("randiance from ls: {randiance_direct}");
                    1.0/(probability_ls*1.0/area) // pdf
                    *intersection.triangle.mat.brdf().component_mul(&randiance_direct) // brdf * radiance
                    *ray_to_ls.direction.dot(&(face_normal/face_normal.norm())) // cos(omega_i, normal_p)
                    *(-ray_to_ls.direction.dot(&(ls_face_normal/ls_face_normal.norm()))) // cos(-omega_i, normal_x)
                    /ls_intersection.ray_scaling.powi(2) // devision by r^2
                },
                _ => {
                    Vector3::new(0.0, 0.0, 0.0)
                }
            };
            trace!("randiance from ls after surface interaction: {randiance_direct}");
            // eval secondary lighting
            let ray_refl =  Ray::new_into_hemisphere(
                point_of_intersection,
                face_normal,
                intersection.ray.n_bounces+1,
                SampleWeighting::Cosine);

            let randiance_indirect = match self.trace(&ray_refl, false){
                Some(randiance_indirect) => {
                    1.0/ray_refl.direction.dot(&(face_normal/face_normal.norm())) // pdf
                    *intersection.triangle.mat.brdf().component_mul(&randiance_indirect) // brdf * radiance
                    *ray_refl.direction.dot(&(face_normal/face_normal.norm())) // cos(omega_i, normal_p)
                },
                None => {
                    Vector3::new(0.0, 0.0, 0.0)
                }
            };
            trace!("randiance from indirect lighting after surface interaction: {randiance_indirect}");
            // combine with randiance_emitted
            randiance_direct+randiance_indirect
        };
        // combine lighting and randiance_emitted
        if let Some(randiance_emitted) = randiance_emitted {
            debug!("Ray-scene-intersection evaluated: radiance = {}", randiance_emitted + lighting);
            randiance_emitted + lighting
        } else {
            debug!("Ray-scene-intersection evaluated: radiance = {}", lighting);
            lighting
        }
    }
    /// Trace ray into the scene
    ///
    /// Intersect ray with the scene, split ray according to stratification rules and evaluate intersections with scene.
    fn trace(&self, ray: &Ray, eval_light_emission: bool) -> Option<Vector3<f64>> {
        info!("Tracing ray: {ray:?}");
        // determine first intersection with an object and evaluate it in case it exists
        if let Some(intersec) = self.bbv_tree_root.intersect(ray) {
            return Some(self.eval_intersection(
                intersec,
                eval_light_emission));
        }
        None
    }
    fn trace_to_lightsource(&self, ray: &Ray, intersection: &objects::Intersection) -> Option<Vector3<f64>> {
        info!("Tracing ray to light source: {ray:?}");
        // check if lightsource is concealed by other object
        if self.bbv_tree_root.is_concealed(ray, intersection.ray_scaling) { // if true, do not do anything
            debug!("Light source is concealed");
            None
        } else { // since is not concealed eval intersection with ls
            debug!("Direct connection between light source and ray origin exists");
            // todo check orientation of ray and ls
            Some(intersection.triangle.mat.light_emission())
        }
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
            info!("Tracing ray number {rays_numbered:?} in new thread");
            let mut color_stack = ColorStack::new(res_x, res_y);
            // total number of pixels
            let n_pixels = (res_x as u64)*(res_y as u64);
            // calc width and height of a pixel
            let pixel_width = scene.view_point.window.size.0/(res_x as f64);
            let pixel_height = scene.view_point.window.size.1/(res_y as f64);
            for i in rays_numbered.clone() {
                // determine associated pixel (start top left in lines to bottom right)
                let pixel_x = (i%res_x as u64) as usize;
                let pixel_y = ((i/res_x as u64)%(res_y as u64)) as usize;

                let ray = scene.view_point
                    .create_rndm_ray_through_pixel(pixel_x, pixel_y, pixel_width, pixel_height);
                debug!("Created ray: {ray:?}");

                // intersect ray with scene objects and eval intersections
                if let Some(ray_data) = scene.trace(&ray, true) {
                    color_stack.accumulated_color[pixel_y][pixel_x] += ray_data;
                }
                // color_stack.inc_height(pixel_x, pixel_y, ray.color);

                if i%n_pixels == n_pixels-1 {
                    bar.inc(1);
                }
            }
            bar.abandon();
            info!("Finished tracing ray number {rays_numbered:?}");
            // send result
            transmitter.send(color_stack).unwrap();
        })
    }
    pub fn look(&self, res_x: u32, res_y: u32, n_rays_per_pixel: u64, n_threads: u8) -> ColorStack {
        info!("Looking into scene with resolution {res_x}x{res_y} and {n_rays_per_pixel} ray per pixel");
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
