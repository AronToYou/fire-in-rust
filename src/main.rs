use std::env;
use indicatif::ProgressIterator;
use fire_in_rust::{new_sim, Params, Field};

fn main() {
    let mut sim = new_sim::<240, 320>(Params {
        h: 1.0,
        dt: 0.5,

        visc: 0.0005,
        k_diff: 0.0001,
        k_buoy: 0.05,
        temp_air: 300.0,
        temp_max: 2200.0,
        k_cool: 0.002,
        vconf: 2.0,
        s: 0.7,
        d_fuel: 1.0,
        d_hgas: 0.1,

        render_scale: 1.0
    });

    // Prints the initialized level set field
    sim.print_field(Field::Phi);
    sim.print_field(Field::U);
    
    let steps: usize = env::args().nth(1).and_then(|s| s.parse::<usize>().ok()).unwrap_or(100);
    println!("Simulating {steps} steps...");
    for _ in (0..steps).progress() {
        sim.step();
    }
    sim.print_field(Field::Phi);
    sim.print_field(Field::U);
    println!("Done!");
}