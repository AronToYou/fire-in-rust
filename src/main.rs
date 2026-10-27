use std::env;
use indicatif::ProgressIterator;
use fire_in_rust::{new_sim, Params, Field};

fn main() {
    let mut sim = new_sim::<120, 160>(Params {
        h: 0.05,    // grid cell length     [m]
        dt: 0.5,    // simulation time step [s]

        visc: 0.0005,      // kinematic viscosity [m^2/s] (for velocity diffusion)
        k_diff: 0.0001,    // diffusion constant  [m^2/s] (for scalar diffusion)
        k_buoy: 0.15,      // buoyancy constant   [N/K]   (for buoyancy force)
        temp_air: 300.0,   // ambient air temperature [K]
        temp_max: 2200.0,  // peak flame temperature  [K]
        k_cool: 3000.0,  // cooling constant        [K/s]
        vconf_f: 16.0,   // vorticity confinement   [1/s] (for fuel)
        vconf_h: 60.0,   // vorticity confinement   [1/s] (for hot gas)
        s: 0.5,         // reaction rate      [m/s] (.44 for propane premixed with oxidizer)
        d_fuel: 1.0,    // denisty of fuel    [kg/m^3]
        d_hgas: 0.01,    // density of hot gas [kg/m^3] (0.1)

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