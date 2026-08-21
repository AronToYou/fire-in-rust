use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_step(c: &mut Criterion) {
	let sim = fire_in_rust::new_sim::<240, 320>(fire_in_rust::Params {
		h: 1.0,
		dt: 0.5,

		visc: 0.0005,
		k_diff: 0.0001,
		k_buoy: 0.05,
		temp_air: 300.0,
		temp_max: 2200.0,
		k_cool: 0.002,
		vconf: 2.0,
		k_react: 0.7,
		d_fuel: 1.0,
		d_hgas: 0.1,

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