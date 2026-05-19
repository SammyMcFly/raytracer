//! ObjectsLib
//!
//! Objects, materials and associated functions for rusty ray tracer
//!
//!

use nalgebra::Vector3;
use tracing::{info, debug};// , error, info, span, trace, warn};

use std::sync::Arc;

use std::f64::consts;

use super::Ray;

// use super::objectslib;


/// Struct representing a triangle in the [[Scene]]
#[derive(Debug, Clone)]
pub struct Triangle {
    /// Vertex 1 of triangle
    pub v1: Arc<Vector3<f64>>,
    /// Vertex 2 of triangle
    pub v2: Arc<Vector3<f64>>,
    /// Vertex 3 of triangle
    pub v3: Arc<Vector3<f64>>,
    /// Normal vector of vertex 1
    pub n1: Arc<Vector3<f64>>,
    /// Normal vector of vertex 2
    pub n2: Arc<Vector3<f64>>,
    /// Normal vector of vertex 3
    pub n3: Arc<Vector3<f64>>,
    // Material of triangle
    pub mat: Arc<Material>,
    // /// Material of vertex 1
    // pub mat1: Arc<Material>,
    // /// Material of vertex 2
    // pub mat2: Arc<Material>,
    // /// Material of vertex 3
    // pub mat3: Arc<Material>,
}

impl Triangle {
    pub fn new(
        p1: Arc<Vector3<f64>>,
        p2: Arc<Vector3<f64>>,
        p3: Arc<Vector3<f64>>,
        n1: Arc<Vector3<f64>>,
        n2: Arc<Vector3<f64>>,
        n3: Arc<Vector3<f64>>,
        mat: Arc<Material>,
    ) -> Self {
        Self { v1: p1, v2: p2, v3: p3, n1, n2, n3, mat, }
    }
    pub fn face_normal(&self, at: &Vector3<f64>) -> Vector3<f64> {
        if *self.n1 == *self.n2 && *self.n1 == *self.n3 {
            return *self.n1;
        }
        todo!("Produces Nan sometimes");
        // calculate the parametrization of "at" for a*self.v1+b*self.v2+c*self.v3
        // for that a linear system needs to be solved
        let f_1: f64 = self.v2.y-self.v1.y*self.v2.x/self.v1.x;
        let d_1 = at.x/self.v1.x-self.v2.x/self.v1.x*at.y/f_1+self.v2.x/self.v1.x/self.v1.x*self.v1.y*at.x/f_1;
        let d_2 = self.v2.x/self.v1.x*(self.v3.y-self.v1.y*self.v3.x/self.v1.x)/f_1-self.v3.x/self.v1.x;
        let d_3 = at.y/f_1-self.v1.y/self.v1.x*at.x/f_1;
        let d_4 = 1.0/f_1*(self.v3.y-self.v1.y*self.v3.x/self.v1.x);
        let c = 1.0/(d_2*self.v1.z-d_4*self.v2.z+self.v3.z)*(at.z-d_1*self.v1.z-d_3*self.v2.z);
        let b = d_3-c*d_4;
        let a = d_1+c*d_2;
        // note: a scales self.v1, b scales self.v2, c scales self.v3 to arrive at point "at"
        a*(*self.n1)+b*(*self.n2)+c*(*self.n3)
    }
    pub fn random_point_on_surface(&self, reference_point: &Vector3<f64>, n_bounces: u64) -> (Vector3<f64>, Ray, f64, f64, f64) {
        // calculate area of triangle
        let mut e1 = Vector3::new(0.0, 0.0, 0.0);
        self.v2.sub_to(&*self.v1, &mut e1);
        let mut e2 = Vector3::new(0.0, 0.0, 0.0);
        self.v3.sub_to(&*self.v1, &mut e2);
        let area = 0.5*e1.cross(&e2).norm();
        // determine random point on lightsource
        let p0 = *self.v1;
        let rndm_number_1 = super::utilities::random_incl();
        let rndm_number_2 = super::utilities::random_incl();
        let rndm_point = p0+rndm_number_1*e1+rndm_number_2*e2;
        // apply point symmetry at middle point between self.v2 and self.v3 to achieve a uniform distribution
        // inside the triangle
        let point_for_symmetry = p0+0.5*e1+0.5*e2;
        let rndm_point = if rndm_number_1+rndm_number_2 > 1.0 {
            let vector_for_symmetry = point_for_symmetry-rndm_point;
            rndm_point+2.0*vector_for_symmetry
        } else {
            rndm_point
        };
        // create Intersection
        let ray_to_rndm_p = Ray::new(
            *reference_point,
            rndm_point-reference_point,
            n_bounces+1);
        (rndm_point, ray_to_rndm_p, area, rndm_number_1, rndm_number_2)
    }
}

#[derive(Debug, Clone)]
pub enum BacksideIntersection {
    #[allow(dead_code)]
    Ignore,
    #[allow(dead_code)]
    EndRay,
}

/// Struct that contains every information important to a intersection
#[derive(Debug, Clone)]
pub struct Intersection<'a> {
    pub ray: &'a Ray,
    pub triangle: &'a Triangle,
    pub ray_scaling: f64,
    pub normal: Vector3<f64>,
}

impl<'a> Intersection<'a> {
    pub fn new(
        ray: &'a Ray,
        triangle: &'a Triangle,
        ray_scaling: f64,
        b1: f64,
        b2: f64,
    ) -> Self {
        Intersection {
            ray,
            triangle,
            ray_scaling,
            normal: Self::get_interpolated_normal(b1, b2, &triangle.n1, &triangle.n2, &triangle.n3),
        }
    }

    pub fn get_interpolated_normal(
        b1: f64,
        b2: f64,
        n1: &Vector3<f64>,
        n2: &Vector3<f64>,
        n3: &Vector3<f64>
    ) -> Vector3<f64> {
        let b0 = 1.0 - b1 - b2;
        (b0 * n1 + b1 * n2 + b2 * n3).normalize()
    }

    pub fn from_front(&self) -> bool {
        self.ray.direction.dot(&self.normal) < 0.0
    }
}

/// Axis that devides a [[BBVTLeafPartition]]
enum PartitioningAxis {
    X,
    Y,
    Z,
}

impl PartitioningAxis {
    fn from(extend: &Vector3<f64>) -> Self {
        if extend.x > extend.y && extend.x > extend.z {
            Self::X
        } else if extend.y > extend.x && extend.y > extend.z {
            Self::Y
        } else {
            Self::Z
        }
    }
}

#[derive(Debug, Clone)]
struct BBVTLeafPartition {
    leafs: Vec<BBVTLeaf>,
    bounds_min: Vector3<f64>,
    bounds_max: Vector3<f64>,
    extend: Vector3<f64>,
    surface: f64,
}

impl BBVTLeafPartition {
    fn from(triangles: Vec<Triangle>,) -> Self {
        let mut triangle_bvols: Vec<BBVTLeaf> = Vec::new();
        for triangle in triangles {
            triangle_bvols.push(BBVTLeaf::from(triangle));
        }
        Self::from_leafs(triangle_bvols)
    }
    fn from_leafs(leafs: Vec<BBVTLeaf>) -> Self {
        // determine bounds_min, bounds_max
        let (bounds_min, bounds_max) = Self::bounds(&leafs);
        // Calculate extend
        let extend = bounds_max-bounds_min;
        // Calculate surface
        let surface = 2.0*extend.x*extend.y + 2.0*extend.x*extend.z + 2.0*extend.y*extend.z;
        Self {
            leafs,
            bounds_min,
            bounds_max,
            extend,
            surface,
        }
    }
    fn bounds(partition: &[BBVTLeaf]) -> (Vector3<f64>, Vector3<f64>) {
        let mut bounds_min = partition[0].bounds_min;
        let mut bounds_max = partition[0].bounds_max;
        for bv in partition.iter().skip(1) {
            if bv.bounds_min.x < bounds_min.x {
                bounds_min.x = bv.bounds_min.x;
            }
            if bv.bounds_max.x > bounds_max.x {
                bounds_max.x = bv.bounds_max.x;
            }
            if bv.bounds_min.y < bounds_min.y {
                bounds_min.y = bv.bounds_min.y;
            }
            if bv.bounds_max.y > bounds_max.y {
                bounds_max.y = bv.bounds_max.y;
            }
            if bv.bounds_min.z < bounds_min.z {
                bounds_min.z = bv.bounds_min.z;
            }
            if bv.bounds_max.z > bounds_max.z {
                bounds_max.z = bv.bounds_max.z;
            }
        }
        (bounds_min, bounds_max)
    }
    fn sort_by_axis(&mut self, axis: &PartitioningAxis) {
        self.leafs.sort_by(|a, b| {
            match axis {
                PartitioningAxis::X => {
                    a.centroid.x.partial_cmp(&b.centroid.x).unwrap()
                },
                PartitioningAxis::Y => {
                    a.centroid.y.partial_cmp(&b.centroid.y).unwrap()
                },
                PartitioningAxis::Z => {
                    a.centroid.z.partial_cmp(&b.centroid.z).unwrap()
                },
            }
        });
    }
    fn len(&self) -> usize {
        self.leafs.len()
    }
    /// Split at optimal position according to surface area heuristic
    fn split_by_sah(self, axis: PartitioningAxis) -> (Self, Self) {
        let mut partition = self.clone();
        // Sort triangle_bvols by partitioning axis's value of centroid
        partition.sort_by_axis(&axis);
        // calculate costs of splitting at certain position and store minimum
        let mut first_partition = partition.leafs;
        let second_partition  = first_partition.split_off(1);
        let mut first_partition = Self::from_leafs(first_partition);
        let mut second_partition = Self::from_leafs(second_partition);
        let mut cost = 1.0/8.0+first_partition.surface*first_partition.len() as f64
            +second_partition.surface*second_partition.len() as f64;
        for i in 2..self.leafs.len() {
            let mut p1 = self.leafs.clone();
            let p2  = p1.split_off(i);
            let p1 = Self::from_leafs(p1);
            let p2  = Self::from_leafs(p2);
            let c = 1.0/8.0+p1.surface*p1.len() as f64
                +p2.surface*p2.len() as f64;
            if c < cost {
                cost = c;
                first_partition = p1;
                second_partition = p2;
            }
        }
        (first_partition, second_partition)
    }
}

/// Node in the binary bounding volume tree (aabb-volume)
#[derive(Debug, Clone)]
struct BBVTLeaf {
    triangle: Triangle,
    bounds_min: Vector3<f64>,
    bounds_max: Vector3<f64>,
    centroid: Vector3<f64>,
}

impl BBVTLeaf {
    fn from(triangle: Triangle) -> Self {
        let x_values = vec![
            triangle.v1.x,
            triangle.v2.x,
            triangle.v3.x,
        ];
        let y_values = vec![
            triangle.v1.y,
            triangle.v2.y,
            triangle.v3.y,
        ];
        let z_values = vec![
            triangle.v1.z,
            triangle.v2.z,
            triangle.v3.z,
        ];
        let bounds_min = Vector3::new(
            x_values.clone().into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
            y_values.clone().into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
            z_values.clone().into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
        );
        let bounds_max = Vector3::new(
            x_values.into_iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
            y_values.into_iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
            z_values.into_iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
        );
        Self {
            triangle,
            bounds_min,
            bounds_max,
            centroid: 0.5*(bounds_min+bounds_max),
        }
    }
}

/// Enum describing all varaints, which the child of a [[BBVTNode]] can be
#[derive(Debug, Clone)]
enum BBVTNodeChild {
    Node(Box<BBVTNode>),
    Triangle(Triangle),
}

impl BBVTNodeChild {
    fn intersect<'a>(&'a self, ray: &'a Ray, backside_intersection: &BacksideIntersection) -> Option<Intersection<'a>> {
        match self {
            BBVTNodeChild::Node(node) => {
                node.intersect(ray, backside_intersection)
            },
            BBVTNodeChild::Triangle(tri) => {
                if let Some((t, b1, b2)) = ray.intersect_triangle(
                    &tri.v1,
                    &tri.v2,
                    &tri.v3,
                ) {
                    match backside_intersection {
                        BacksideIntersection::Ignore => {
                            let intersection = Intersection::new(
                                ray,
                                tri,
                                t,
                                b1, b2,
                            );
                            if intersection.from_front() {
                                return Some(Intersection::new(ray, tri, t, b1, b2));
                            }
                        },
                        BacksideIntersection::EndRay => {
                            return Some(Intersection::new(ray, tri, t, b1, b2));
                        },
                    }
                }
                None
            },
        }
    }
    fn is_concealed(&self, ray: &Ray, ray_scaling: f64) -> bool {
        match self {
            BBVTNodeChild::Node(node) => {
                return node.is_concealed(ray, ray_scaling);
            },
            BBVTNodeChild::Triangle(tri) => {
                if let Some((t, b1, b2)) = ray.intersect_triangle(
                    &tri.v1,
                    &tri.v2,
                    &tri.v3,
                ) {
                    let interpolated_n = Intersection::get_interpolated_normal(b1, b2, &tri.n1, &tri.n2, &tri.n3);
                    // Test if hitting from the direction the normal points in
                    let is_front_face = ray.direction.dot(&interpolated_n) < 0.0;
                    if !is_front_face || t < ray_scaling {
                        return true;
                    }
                }
            },
        }
        false
    }
}

/// Node in the binary bounding volume tree (aabb-volume)
#[derive(Debug, Clone)]
pub struct BBVTNode {
    child1: BBVTNodeChild,
    child2: BBVTNodeChild,
    /// Boundaries of the bounding volume
    pub bounds_min: Vector3<f64>,
    pub bounds_max: Vector3<f64>,
}

impl BBVTNode {
    /// Build binary bounding volume tree using the surface area heuristic
    pub fn build_bbv_tree(triangles: Vec<Triangle>) -> Self {
        // Init triangle bounding volumes
        let triangle_bvols = BBVTLeafPartition::from(triangles);

        Self::build_node(triangle_bvols)
    }
    fn build_node(leaf_partition: BBVTLeafPartition) -> Self {
        // Choose partitioning axis (largest extend)
        let axis = PartitioningAxis::from(&leaf_partition.extend);
        // pass on bounds (before move)
        let bounds_min = leaf_partition.bounds_min;
        let bounds_max = leaf_partition.bounds_max;
        // calculate costs of splitting at certain position
        let (partition1, partition2) = leaf_partition.split_by_sah(axis);

        let child1 = if partition1.leafs.len() >= 2 {
            BBVTNodeChild::Node(Box::new(Self::build_node(partition1)))
        } else {
            BBVTNodeChild::Triangle(partition1.leafs.first().unwrap().triangle.clone())
        };
        let child2 = if partition2.leafs.len() >= 2 {
            BBVTNodeChild::Node(Box::new(Self::build_node(partition2)))
        } else {
            BBVTNodeChild::Triangle(partition2.leafs.first().unwrap().triangle.clone())
        };
        Self {
            child1,
            child2,
            bounds_min,
            bounds_max,
        }
    }
    pub fn intersect<'a>(
        &'a self,
        ray: &'a Ray,
        backside_intersection: &BacksideIntersection,
    ) -> Option<Intersection<'a>> {
        info!("Intersecting ray with BV");
        // test if ray enters bounding volume
        if ray.intersect_bv(self.bounds_min, self.bounds_max) {
            debug!("Intersection found");
            // if yes: Test Children
            match (self.child1.intersect(ray, backside_intersection), self.child2.intersect(ray, backside_intersection)) {
                (Some(intersec_1), Some(intersec_2)) => {
                    if intersec_1.ray_scaling < intersec_2.ray_scaling {
                        return Some(intersec_1);
                    } else {
                        return Some(intersec_2);
                    }
                },
                (Some(intersec_1), None) => {
                    return Some(intersec_1);
                },
                (None, Some(intersec_2)) => {
                    return Some(intersec_2);
                },
                (None, None) => {
                    return None;
                },
            }
        }
        debug!("No intersection with BV found");
        // if no: Return None
        None
    }
    pub fn is_concealed(&self, ray: &Ray, ray_scaling: f64) -> bool {
        // test if ray enters bounding volume
        if ray.intersect_bv(self.bounds_min, self.bounds_max) {
            // if yes: Test Children
            if self.child1.is_concealed(ray, ray_scaling) || self.child2.is_concealed(ray, ray_scaling) {
                return true;
            }
        }
        // if no: Return false
        false
    }
}

/// Material
#[derive(Debug, Clone)]
pub enum Material {
    /// Matt material
    Diffuse {
        refl: Reflectance,
    },
    /// Specular material
    Specular {
        refl: Reflectance,
        e: u8,
    },
    /// Combination of matt and specular Material
    SpecularDiffuse {
        refl: Reflectance,
        e: u8,
        alpha: f64, // combine beta and alpha?
        beta: f64,
    },
    DiffuseLightSource {
        refl: Reflectance,
        light_emission: super::Color,
    }
}

impl Material {
    pub fn brdf(&self) -> Vector3<f64> {
        match self {
            Self::Diffuse { refl } => refl/consts::PI,
            Self::Specular { refl, e } => refl/consts::PI, // todo
            Self::SpecularDiffuse { refl, e, alpha, beta } => refl/consts::PI, // todo
            Self::DiffuseLightSource { refl, .. } => refl/consts::PI,
        }
    }
    pub fn is_light_source(&self) -> bool {
        match self {
            Self::DiffuseLightSource { .. } => true,
            _ => false,
        }
    }
    pub fn light_emission(&self) -> super::Color {
        match self {
            Self::DiffuseLightSource { light_emission, .. } => {
                light_emission.clone()
            },
            _ => panic!("Called light emission on non-light-source"),
        }
    }
}

pub type Reflectance = Vector3<f64>;

// /// Absorption/Reflectance/Color
// #[derive(Debug, Clone)]
// pub struct Reflectance {
//     pub r: f64,
//     pub g: f64,
//     pub b: f64,
// }

// impl Reflectance {
//     /// Reflectance constructor
//     pub fn new(r: f64, g: f64, b: f64) -> Self {
//         Self { r, g, b, }
//     }
// }


// Reflection Pattern
// Transmission Pattern (refractive index)


// Scene creation lib
