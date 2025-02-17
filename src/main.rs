//! Rusty ray tracer
//!
//! Simple ray tracer that supports loading .obj files and rendering into a .png file
//!
//!
use clap::Parser;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::FmtSubscriber;
use tracing::info; // , error, info, span, trace, warn, debug};

use nalgebra::Vector3;

use exr::prelude::WritableImage;

use std::io::Error;

mod rtcore;



/// Simple ray tracer that supports loading .obj files and rendering into a .png file
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// File path to input .obj file
    #[arg(short, long, default_value = "./test_scene.obj")]
    object: String,
    // /// Horizontal resolution
    // /// Vertical resolution
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

    // create scene
    let window = rtcore::ViewWindow::from(
        Vector3::from([1.0, 1.0, -1.0]),
        (16.0/8.0, 9.0/8.0));

    let view_point = rtcore::ViewPoint {
        point: Vector3::from([-8.0, -8.0, 8.0]),
        window,
        distance_to_window_plane: 2.0,
    };
    let scene1 = rtcore::Scene::from(view_point, models, materials.unwrap(), 1);

    // init img properties
    let imgx = 2560;
    let imgy = 1440;
    // let n_rays = 1E9 as u64;
    let n_rays_per_pixel = 500;
    let n_threads: u8 = 10;

    // take look into scene
    let colorstack = scene1.look(imgx, imgy, n_rays_per_pixel, n_threads);

    // Write data to .png file
    // let _png: image::ImageBuffer<image::Rgb::<u8>, Vec<u8>> = colorstack.into();
    // _png.save("rendering.png").unwrap();

    // Write data to .exr file
    let _exr = colorstack.into_exr();
    _exr.write()
        .to_file("rendering.exr").unwrap();

    Ok(())
}
