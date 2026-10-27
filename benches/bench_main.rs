use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_step(c: &mut Criterion) {
	let sim = fire_in_rust::new_sim::<240, 320>(fire_in_rust::Params {
		h: 0.05,    // grid cell length     [m]
        dt: 0.5,    // simulation time step [s]

        visc: 0.0005,      // kinematic viscosity [m^2/s] (for velocity diffusion)
        k_diff: 0.0001,    // diffusion constant  [m^2/s] (for scalar diffusion)
        k_buoy: 0.05,      // buoyancy constant   [N/K]   (for buoyancy force)
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

	c.bench_function("240x320 Simulation Step", |b| {
		
		let mut sim_copy = sim.clone();

		b.iter(|| {
			black_box(&mut sim_copy).step();
		});
	});
}

criterion_group!(benches, bench_step);
criterion_main!(benches);