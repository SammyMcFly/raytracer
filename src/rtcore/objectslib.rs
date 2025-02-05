//! ObjectsLib
//!
//! Objects, materials and associated functions for rusty ray tracer
//!
//!
// use tracing::{error, debug}; // , error, info, span, trace, warn};

use nalgebra::Vector3;

use super::Ray;


/// Struct representing a triangle in a [[Model]]
struct Triangle<T: num_traits::Num, U: Material> {
    /// Vertex 1 of triangle
    v1: Vector3<T>,
    /// Vertex 2 of triangle
    v2: Vector3<T>,
    /// Vertex 3 of triangle
    v3: Vector3<T>,
    /// Normal vector of vertex 1
    n1: Vector3<T>,
    /// Normal vector of vertex 2
    n2: Vector3<T>,
    /// Normal vector of vertex 3
    n3: Vector3<T>,
    /// Material of vertex 1
    mat1: U,
    /// Material of vertex 2
    mat2: U,
    /// Material of vertex 3
    mat3: U,
}

impl<T: num_traits::Num, U: Material> Triangle<T, U> {
    fn from(
        p1: Vector3<T>,
        p2: Vector3<T>,
        p3: Vector3<T>,
        n1: Vector3<T>,
        n2: Vector3<T>,
        n3: Vector3<T>,
        mat1: U,
        mat2: U,
        mat3: U,
    ) -> Self {
        Self { v1: p1, v2: p2, v3: p3, n1, n2, n3, mat1, mat2, mat3, }
    }
}

/// Struct representing a Model for Ray Tracing
///
/// Container for the necessary model characteristics
struct Model<T: num_traits::Num, U: Material> {
    surface: Vec<Triangle<T, U>>,
    bounding_volume: Vec<Triangle<T, U>>,
}

/// Material
trait Material {
    fn eval_ray_intersection();
}

/// Matt material implementation of a [[Material]]
struct MattMaterial {
    refl: Reflectance,
}

/// Specular material implementation of a [[Material]]
struct SpecMaterial {
    refl: Reflectance,
    e: u8,
}

/// Combination of matt and specular Material implementation of a [[Material]]
struct MattSpecMaterial {
    refl: Reflectance,
    e: u8,
    alpha: f64, // combine beta and alpha?
    beta: f64,
}

/// Absorption
pub struct Reflectance {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Reflectance {
    /// Reflectance constructor
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b, }
    }
}


// Reflection Pattern
// Transmission Pattern (refractive index)


// Scene creation lib
