/// Rusty ray tracer
///
/// Simple ray tracer that supports loading .obj files and rendering into a .png file
///
use clap::Parser;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::FmtSubscriber;
use tracing::info; // , error, info, span, trace, warn, debug}
use nalgebra::Vector3;
use exr::prelude::WritableImage;
use std::io::Error;
use std::fs;

mod rtcore;
mod config;

use rtcore::{PathLength, objects::BacksideIntersection};
use config::{Config, PathLengthConfig};

/// Simple ray tracer that supports loading .obj files and rendering into a .png file
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// File path to input .obj file with paired .mtl file
    ///
    /// Only triangulated object are supported, .mtl file needs to contain a material called: light or Light,
    /// which will be used as light source
    #[arg(short, long, default_value = "./cornell-box.obj")]
    object: String,
    /// Assume graphical coordinates when importing scene from .obj file (y-axis is vertical)
    #[arg(short, long, default_value_t = true)]
    graphical_coordinates: bool,
    /// Rendering configuration
    #[arg(short, long, default_value = "./rendering.toml")]
    config: String,
    /// Log severity level (Options: TRACE, DEBUG, INFO, WARN, ERROR, OFF)
    #[arg(short, long, default_value_t=String::from("OFF"))]
    log: String,
}

fn main() -> Result<(), Error>{
    // parse args
    let args = Args::parse();
    // init logging
    let severity_level = match &args.log[..] {
        "TRACE" => LevelFilter::TRACE,
        "DEBUG" => LevelFilter::DEBUG,
        "INFO" => LevelFilter::INFO,
        "WARN" => LevelFilter::WARN,
        "ERROR" => LevelFilter::ERROR,
        _ => LevelFilter::OFF, // String::from("OFF")
    };
    let subscriber = FmtSubscriber::builder()
        .with_max_level(severity_level)
        .with_writer(std::io::stdout)
        .with_line_number(true)
        // .with_ansi(false)
        // .pretty()
        .finish();
        // .with(debug_log);
    tracing::subscriber::set_global_default(subscriber).expect("Failed to set subscriber");
    info!("Start rusty ray tracer");

    // import .obj file
    let (models, materials) =
        tobj::load_obj(
            &args.object,
            &tobj::LoadOptions::default()
        )
        .expect("Failed to .obj load file");

    for model in &models {
        println!("Model name: {}", model.name);
        println!("Model positions: {:?}", model.mesh.positions);
        println!("Model indices: {:?}", model.mesh.indices);
        println!("Face normals: {:?}", model.mesh.normals);
        println!("Face normal indices: {:?}", model.mesh.normal_indices);
    }

    let materials = match materials {
        Ok(m) => { Some(m) },
        Err(_) => { println!("Failed to load MTL file"); None },
    };

    // Load config
    let config_str = fs::read_to_string(&args.config)
        .expect("Failed to read config file");
    let config: Config = toml::from_str(&config_str)
        .expect("Failed to parse config file");

    // Build scene parameters from config
    let window = rtcore::ViewWindow::from(
        Vector3::from(config.camera.view_direction),
        (config.camera.window_size[0], config.camera.window_size[1]),
    );

    let view_point = rtcore::ViewPoint {
        point: Vector3::from(config.camera.position),
        window,
        distance_to_window_plane: config.camera.focal_distance,
    };

    let path_length: PathLength = match config.render.path_length {
        PathLengthConfig::Fixed { max_bounces } => rtcore::PathLength::Fixed { max_bounces },
        PathLengthConfig::Adaptive { min_bounces, termination_probability } => {
            rtcore::PathLength::Adaptive { min_bounces, termination_probability }
        },
    };

    let backside_intersections = match config.render.backside_intersections.as_str() {
        "end_ray" => BacksideIntersection::EndRay,
        _ => BacksideIntersection::Ignore,
    };

    // Render

    let scene1 = rtcore::Scene::from(
        view_point,
        models,
        materials.unwrap(),
        args.graphical_coordinates,
    );

    // take look into scene
    let colorstack = scene1.look(
        config.render.width,
        config.render.height,
        config.render.rays_per_pixel,
        backside_intersections,
        path_length,
        config.render.channel_bound,
        config.render.tile_size);

    // Write data to .png file
    // let _png: image::ImageBuffer<image::Rgb::<u8>, Vec<u8>> = colorstack.into();
    // _png.save("rendering.png").unwrap();

    // Write data to .exr file
    let _exr = colorstack.into_exr();
    _exr.write()
        .to_file("rendering.exr").unwrap();

    Ok(())
}
