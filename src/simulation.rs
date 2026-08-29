#[cfg(debug_assertions)]
use std::time::{SystemTime, UNIX_EPOCH};
use crate::grid::{P, GridDisp, Linterp, IsNan};
use rayon::prelude::*;
const MIN_NORM: f32 = 1e-8;


// --------------------------------------------- Utility Functions ---------------------------------------------
pub fn new_sim<const NX: usize, const NY: usize>(p: Params) -> Sim<NX, NY, impl Fn((f32, f32)) -> (f32, f32) + Clone> {
    let maxx = (NX as f32) - 1.001;
    let maxy = (NY as f32) - 1.001;
    let clamp_xy = move |(x, y): (f32, f32)| (x.clamp(0.0, maxx), y.clamp(0.0, maxy));
    Sim::<NX, NY, _>::new(p, clamp_xy)
}


// ------------------------------------------- Simulation Parameters -------------------------------------------
#[derive(Clone, Copy)]
pub struct Params {
    // Grid Parameters //
    pub h: f32,     // grid cell length     [m]
    pub dt: f32,    // simulation time step [s]

    // Physical Constants //
    pub visc: f32,      // kinematic viscosity [m^2/s] (for velocity diffusion)
    pub k_diff: f32,    // diffusion constant  [m^2/s] (for scalar diffusion)
    pub k_buoy: f32,    // buoyancy constant   [N/K]   (for buoyancy force)
    pub temp_air: f32,  // ambient air temperature [K]
    pub temp_max: f32,  // peak flame temperature  [K]
    pub k_cool: f32,    // cooling constant        [K/s]
    pub vconf: f32,     // vorticity confinement   [1/s]
    pub s: f32,         // reaction rate      [m/s]
    pub d_fuel: f32,    // denisty of fuel    [kg/m^3]
    pub d_hgas: f32,    // density of hot gas [kg/m^3]

    // scaling factor for rendering //
    pub render_scale: f32
}


// --------------------------------------------- Simulation State ---------------------------------------------
#[derive(Clone)]
pub struct Sim<const NX: usize, const NY: usize, C> where C: Fn((f32, f32)) -> (f32, f32) + Clone {
    param: Params,  // Simulation parameters (defined above)
    clamp_xy: C,   // clamping function for coordinates

    u: Box<[[P<f32>; NY]; NX]>,      // velocity field (x, y)
    p: Box<[[f32; NY]; NX]>,        // pressure field
    div_u: Box<[[f32; NY]; NX]>,   // divergence of velocity field (∇·u)

    phi: Box<[[f32; NY]; NX]>,        // level set (+pos in fuel region, -neg outside, 0 at boundary)
    temp_gas: Box<[[f32; NY]; NX]>,  // temperature field (hot gas domain)
    rt: Box<[[f32; NY]; NX]>,       // reaction-time tracker (1 at fuel; decreases after crossing)
    dns: Box<[[f32; NY]; NX]>,     // smoke density (simple)

    tmp: Box<[[f32; NY]; NX]>, tmp2: Box<[[P<f32>; NY]; NX]>,  // temporary intermediate fields for calculations
    
    #[cfg(debug_assertions)]
    ms_counter: u128  // time counter for NaN checking [ms]
}

/// Which field from `Sim` to print
pub enum Field { 
    U,
    P, DivU,
    Phi, Temp, Rt, Dns
}

impl<const NX: usize, const NY: usize, C> Sim<NX, NY, C> where C: Fn((f32, f32)) -> (f32, f32) + Clone {
    fn new(param: Params, clamp_xy: C) -> Self {
        let mut s = Self {
            param, clamp_xy,
            u: Box::new([[P(0.0, 0.0); NY]; NX]),
            p: Box::new([[0.0; NY]; NX]),
            div_u: Box::new([[0.0; NY]; NX]),

            phi: Box::new([[1.0; NY]; NX]),
            temp_gas: Box::new([[param.temp_air; NY]; NX]),
            rt: Box::new([[0.0; NY]; NX]),
            dns: Box::new([[0.0; NY]; NX]),

            tmp: Box::new([[0.0; NY]; NX]), tmp2: Box::new([[P(0.0, 0.0); NY]; NX]),

            #[cfg(debug_assertions)]
            ms_counter: 0
        };
        s.init_fuel_inlet();
        s
    }

    /// Initialize fuel inlet at bottom of domain
    fn init_fuel_inlet(&mut self) {
        for x in 0..NX {
            for y in 0..NY {
                if y < 10 {
                    self.phi[x][y] = 5.0;
                    self.temp_gas[x][y] = self.param.temp_air + 10.0;
                    self.u[x][y] = if y == 0 { P(0.0, 1.5) } else { P(0.0, 0.0) };
                } else {
                    self.phi[x][y] = -5.0;
                }
            }
        }
    }

    /// Apply boundary conditions to velocity field
    pub fn apply_boundary_conditions(&mut self) {
        for x in 0..NX {
            self.u[x][0].1 = 1.5;
            self.u[x][NY-1].1 = 0.0;
        }
        for y in 0..NY {
            self.u[0][y].0 = 0.0;
            self.u[NX-1][y].0 = 0.0;
        }
    }

    /// Perform single full simulation step 
    pub fn step(&mut self) {
        // A) Update level set field //
        self.update_levelset();
        self.check_for_nans("after updating level set");

        // B) Intermediate Velocity calculated //
        self.add_forces();                    // 1. Add Forces {Bouyancy, vorticity confinement}
        self.check_for_nans("after adding forces");
        self.apply_boundary_conditions();    //
        self.semi_lagrangian_advect();      // 2. Semi-Lagrangian advection of velocity fields
        self.check_for_nans("after semi-lagrangian advection");
        self.apply_boundary_conditions();

        // C) Apply Pressure Gradient //
        self.compute_divergence();         // 1. Compute divergence of intermediate velocity field
        self.solve_for_pressure();        // 2. Jacobi iteration to solve Poisson equation for pressure field
        self.check_for_nans("after solving for pressure");
        self.apply_pressure_gradient();  // 3. Apply gradient of pressure field to velocity field
        self.check_for_nans("after applying pressure gradient");
        self.apply_boundary_conditions();
    }


    // ---------------------------------- A) Thin-flame Level Set Propagation ----------------------------------
    /// Updates the level set using upwind one-sided differencing to estimate spatial derivatives
    fn update_levelset(&mut self) {
        let (u, phi) = (&*self.u, &*self.phi);
        let (s, dt, h) = (self.param.s, self.param.dt, self.param.h);
        self.tmp.copy_from_slice(phi);  // @TODO
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // A) 2. (unscaled) Central differencing for normed gradient (∇φ/|∇φ|) //
                let gx = phi[x+1][y] - phi[x-1][y];  // gradient x-component
                let gy = phi[x][y+1] - phi[x][y-1];  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(MIN_NORM);  // gradient norm
                if norm.is_nan() {
                    panic!("NaN detected in level set gradient norm at ({}, {})", x, y);
                }
                
                // A) 3. Velocity of implicit surface (where φ==0) //
                let P(wx, wy) = u[x][y] + P(gx, gy)*(s/norm);
                
                // A) 4.1 (unscaled) Upwind one-sided differencing for gradient (∇φ) //
                let ddx = if wx > 0.0 {
                    phi[x][y] - phi[x-1][y]
                } else {
                    phi[x+1][y] - phi[x][y]
                };
                let ddy = if wy > 0.0 {
                    phi[x][y] - phi[x][y-1]
                } else {
                    phi[x][y+1] - phi[x][y]
                };

                // A) 4.2 (scaled) Application of time derivative //
                self.tmp[x][y] = phi[x][y] - (wx*ddx + wy*ddy)*(dt/h);
            }
        }
        std::mem::swap(&mut *self.phi, &mut *self.tmp);
    }


    // ------------------------------- B) Velocity Update via Stam's 4-step loop -------------------------------
    /// B) 1. Addition of Bouyancy and Vorticity Confinement effects to velocity field
    fn add_forces(&mut self) {
        let (u, force, vorticity) = (&*self.u, &mut *self.tmp2, &mut *self.tmp);
        let (k_buoy, temp_air) = (self.param.k_buoy, self.param.temp_air);
        for x in 0..NX {
            for y in 0..NY {
                // B) 1.1 Buoyancy force α(T - T_air)ŷ //
                force[x][y] = P(0.0, k_buoy*(self.temp_gas[x][y] - temp_air));

                // B) 1.2.1 (unscaled) Vorticity ω //
                let P(_, dv_dx) = match x {
                    0 =>             u[x+1][y],
                    _ if x < NX-1 => u[x+1][y] - u[x-1][y],
                    _ =>                       - u[x-1][y],
                };
                let P(du_dy, _) = match y {
                    0 =>             u[x][y+1],
                    _ if y < NY-1 => u[x][y+1] - u[x][y-1],
                    _ =>                       - u[x][y-1],
                };
                vorticity[x][y] = dv_dx - du_dy;
            }
        }
        let (vconf, h, dt) = (self.param.vconf, self.param.h, self.param.dt);
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // B) 1.2.2 (unscaled) Central differencing for normed gradient N = (∇|ω|/|∇|ω||) //
                let gx = vorticity[x+1][y].abs() - vorticity[x-1][y].abs();  // gradient x-component
                let gy = vorticity[x][y+1].abs() - vorticity[x][y-1].abs();  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(MIN_NORM);  // gradient norm

                // B) 1.2.3 (scaled) Force of vorticity confinement εh(N x ω) //
                force[x][y] += P(-gy, gx)*vorticity[x][y]*(vconf*h/norm);

                // B) 1.3 Add force //
                self.u[x][y] += force[x][y]*dt;
            }
        }
        // Don't forget to add bouyancy force to the left/right boundaries (skip top/bottom due to BCs) //
        for x in [0, NX-1] {
            for y in 1..NY-1 {
                self.u[x][y] += force[x][y]*dt;
            }
        }
    }


    // -------------------------- B) 2. Semi-Lagrangian advection of velocity fields --------------------------
    /// Runge-Kutta 2-stage backtrace, bilinear velocity sampling, clamped at boundaries
    fn semi_lagrangian_advect(&mut self) {
        let (u, phi) = (&*self.u, &*self.phi);
        let (dt, h) = (self.param.dt, self.param.h);
        for x in 0..NX {
            for y in 0..NY {
                let P(x1, y1) = P(x as f32, y as f32);

                // Proceed accordingly to whether implicit surface is crossed //
                if self.sample_bilin(phi, (x1, y1)) > 0.0 {  // already in fuel region initially
                    // B) 2.1 backtrace half-step to midpoint //
                    let P(x0, y0) = P(x1, y1) - u[x][y]*(0.5*dt/h);
                    let u0 = self.sample_bilin(u, (x0, y0));  // B) 2.2 midpoint velocity

                    // B) 2.3 backtrace full step //
                    let P(x0, y0) = P(x1, y1) - u0*(dt/h);
                    self.tmp2[x][y] = self.sample_bilin(u, (x0, y0));  // B) 2.4 final velocity

                } else {
                    // B) 2.1 backtrace half-step to midpoint //
                    let P(x0, y0) = P(x1, y1) - u[x][y]*(0.5*dt/h);
                    let phi_0 = self.sample_bilin(phi, (x0, y0));
                    let u0 = if phi_0 > 0.0 {  // if boundary crossed...
                        self.sample_ghost_velocity(P(x0, y0))  // ...sample ghost velocity
                    } else {
                        self.sample_bilin(u, (x0, y0))
                    };  // B) 2.2 midpoint velocity

                    // B) 2.3 backtrace full step //
                    let P(x0, y0) = P(x1, y1) - u0*(dt/h);
                    let phi_0 = self.sample_bilin(phi, (x0, y0));
                    self.tmp2[x][y] = if phi_0 > 0.0 {  // if boundary crossed...
                        self.sample_ghost_velocity(P(x0, y0))  // ...sample ghost velocity
                    } else {
                        self.sample_bilin(u, (x0, y0))
                    };  // B) 2.4 final velocity
                }
            }
        }
        std::mem::swap(&mut *self.u, &mut *self.tmp2);
    }


    // ---------------------------------- C) Apply Pressure Gradient ----------------------------------
    // C) 1. Compute (scaled) divergence of intermediate velocity field (∇·u)(h/8hΔt)
    fn compute_divergence(&mut self) {
        let (u, phi, div_u) = (&*self.u, &*self.phi, &mut *self.div_u);
        let (d_fuel, d_hgas) = (self.param.d_fuel, self.param.d_hgas);
        let scale = 0.125*self.param.h / (self.param.dt);  // (1/4)(h²/4Δt)/(2h)

        let m = d_fuel*self.param.s;  // Mass flux
        let corr = self.param.dt*m*m*(1.0/d_hgas - 1.0/d_fuel);  // Correction term for pressure when sampling across the implicit surface
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // Calculate negative(positive) correction term when sampling hot-gas(fuel) region from fuel(hot-gas) region
                // @NOTE the sign is backwards since the divergence is subtracted from p
                let (d, count) = if phi[x][y] > 0.0 {
                    (d_fuel,
                    ((phi[x+1][y] < 0.0) as i32) + ((phi[x-1][y] < 0.0) as i32) + ((phi[x][y+1] < 0.0) as i32) + ((phi[x][y-1] < 0.0) as i32))
                } else {
                    (d_hgas,
                    -((phi[x+1][y] > 0.0) as i32) - ((phi[x-1][y] > 0.0) as i32) - ((phi[x][y+1] > 0.0) as i32) - ((phi[x][y-1] > 0.0) as i32))
                };

                // Calculate divergence using central differencing
                let du_dx = u[x+1][y].0 - u[x-1][y].0;
                let dv_dy = u[x][y+1].1 - u[x][y-1].1;
                div_u[x][y] = d * (du_dx + dv_dy) * scale + corr*(count as f32);
            }
        }

        //// Handle boundary conditions for divergence at edges of the grid
        // Corners always remain zero
        // y boundaries (top and bottom)
        for x in 1..NX-1 {
            for y in [0, NY-1] {
                let d = if phi[x][y] > 0.0 { d_fuel } else { d_hgas };
                div_u[x][y] = d * (u[x+1][y].0 - u[x-1][y].0) * scale;
            }
        }
        // x boundaries (left and right)
        for x in [0, NX-1] {
            for y in 1..NY-1 {
                let d = if phi[x][y] > 0.0 { d_fuel } else { d_hgas };
                div_u[x][y] = d * (u[x][y+1].1 - u[x][y-1].1) * scale;
            }
        }
    }

    // C) 2. Jacobi iteration to solve Poisson equation for pressure field (∇²p/ρ = (∇·u)/Δt)
    /// Solves velocity diffusion using the Jacobi method on a 2D grid.
    pub fn solve_for_pressure(&mut self) {
        let (div_u, p, p_next) = (&*self.div_u, &mut *self.p, &mut *self.tmp);
        
        let max_iter = 5000;  // maximum number of Jacobi iterations
        let tolerance = 1e-6;  // convergence tolerance for residual error
        // Main Iterative Loop
        for _iter in 0..max_iter {
            let mut max_residual = 0.0f32;

            // Iterate strictly over interior nodes to preserve boundary conditions
            for x in 1..NX-1 {
                for y in 1..NY-1 {

                    // Jacobi Update Formula: p = 0.25 * (Sum(Neighbors) - h^2)
                    p_next[x][y] = 0.25*(p[x+1][y] + p[x-1][y] + p[x][y+1] + p[x][y-1]) - div_u[x][y];

                    // Calculate the local residual error to track convergence
                    let residual = (p_next[x][y] - p[x][y]).abs();
                    if residual > max_residual {
                        max_residual = residual;
                    }
                }
            }

            // Handle boundaries separately to maintain Neumann boundary conditions (∂p/∂n=0)
            for x in 1..NX-1 {
                p_next[x][0] = 0.25*(p[x+1][0] + p[x-1][0]) - div_u[x][0];
                p_next[x][NY-1] = 0.25*(p[x+1][NY-1] + p[x-1][NY-1]) - div_u[x][NY-1];
            }
            for y in 1..NY-1 {
                p_next[0][y] = 0.25*(p[0][y+1] + p[0][y-1]) - div_u[0][y];
                p_next[NX-1][y] = 0.25*(p[NX-1][y+1] + p[NX-1][y-1]) - div_u[NX-1][y];
            }
            // Handle corners separately to approximate Neumann boundary conditions
            p_next[0][0] =       (p[1][0] +       p[0][1])/2.0;
            p_next[0][NY-1] =    (p[1][NY-1] +    p[0][NY-2])/2.0;
            p_next[NX-1][0] =    (p[NX-2][0] +    p[NX-1][1])/2.0;
            p_next[NX-1][NY-1] = (p[NX-2][NY-1] + p[NX-1][NY-2])/2.0;

            // Efficiently swap memory pointers without reallocating arrays
            std::mem::swap(p, p_next);

            // Check for early convergence termination
            if max_residual < tolerance {
                break;
            }
        }
    }

    // C) 3. Apply gradient of pressure field to velocity field (u = u* - Δt∇p/ρ)
    pub fn apply_pressure_gradient(&mut self) {
        let (u, p) = (&mut *self.u, &*self.p);
        let (dt, h, d_fuel, d_hgas) = (self.param.dt, self.param.h, self.param.d_fuel, self.param.d_hgas);
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                let rho = if self.phi[x][y] > 0.0 { d_fuel } else { d_hgas };
                u[x][y].0 -= dt*(p[x+1][y] - p[x-1][y])/(2.0*h*rho);
                u[x][y].1 -= dt*(p[x][y+1] - p[x][y-1])/(2.0*h*rho);
            }
        }
    }
    
    // ------------------------------------------- Utility Functions -------------------------------------------
    /// Biinear sampling of hot gas 'ghost velocity' within the fuel region
    fn sample_ghost_velocity(&self, p: P<f32>) -> P<f32> {
        let P(x, y) = p;
        let (phi, u) = (&*self.phi, &*self.u);
        let (d_h, d_f, s) = (self.param.d_hgas, self.param.d_fuel, self.param.s);

            // (unscaled) Central differencing for normed gradient at non-integer coordinate
            let nx = self.sample_bilin(phi, (x+1.0, y)) - self.sample_bilin(phi, (x-1.0, y));
            let ny = self.sample_bilin(phi, (x, y+1.0)) - self.sample_bilin(phi, (x, y-1.0));
            let norm = (nx*nx + ny*ny).sqrt().max(MIN_NORM);
            let n = P(nx, ny)*(1.0/norm);  // implicit surface unit normal vector
            
            let u_ghost = self.sample_bilin(u, (x, y)) + n*(d_f/d_h - 1.0)*s;
            
            u_ghost
    }

    /// Bilinear sampling of scalar field, clamped at boundaries
    fn sample_bilin<T: Linterp + IsNan + std::fmt::Display>(&self, f: &[[T; NY]; NX], p: (f32, f32)) -> T {
        let (x, y) = (self.clamp_xy)(p);
        let P(x0, y0) = P(x, y).floor();
        let P(tx, ty) = P(x, y) - P(x0 as f32, y0 as f32);

        let f00 = f[x0][y0];           let f01 = f[x0][y0+1];
        let f10 = f[x0+1][y0];         let f11 = f[x0+1][y0+1];
        let a = f00*(1.0-tx) + f10*tx; let b = f01*(1.0-tx) + f11*tx;
        let samp = a*(1.0-ty) + b*ty;
        if cfg!(debug_assertions) && samp.is_nan() {
            panic!("NaN detected in bilin p:({},{}) clamp(p):({},{}) p0:({},{}) t:({},{}) f:[{},{},{},{}] a:{} b:{} samp:{}", 
                    p.0, p.1, x, y, x0, y0, tx, ty, f00, f01, f10, f11, a, b, samp);
        }
        samp
    }

    /// Check all fields for NaN values and panic if any are found
    pub fn check_for_nans(&mut self, stage: &str) {
        if cfg!(debug_assertions) {
            let start = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_millis();
            let fields = [&*self.p, &*self.div_u, &*self.phi];
            let u = &*self.u;
            rayon::join(|| {
                fields.par_iter().enumerate().for_each(|(i, field)| {
                    field.iter().for_each(|row| {
                            row.iter().any(|val| val.is_nan()).then(|| panic!("NaN detected in field {}", i));
                    });
                });
            }, || {
                u.par_iter().for_each(|row| {
                    row.iter().any(|val| val.is_nan()).then(|| panic!("NaN detected in velocity field"));
                });
            });
            let end = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_millis();
            self.ms_counter += end - start;
            println!("NaN check completed in {} ms", self.ms_counter);
        }
    }

    /// Print a downsampled version of a field to the console
    pub fn print_field(&self, which: Field) {
        let (label, field): (&str, &dyn GridDisp) = match which {
            Field::U => ("Velocity (ux, uy)", &*self.u),
            Field::P => ("Hot Gas Pressure", &*self.p),
            Field::DivU => ("Divergence of Hot Gas Pressure", &*self.div_u),
            Field::Phi => ("Level Set", &*self.phi),
            Field::Temp => ("Temperature", &*self.temp_gas),
            Field::Rt => ("Reaction Parameter", &*self.rt),
            Field::Dns => ("Smoke Density", &*self.dns)
        };
        let (stride_x, stride_y) = ((NX/30).max(1), (NY/30).max(1));
        println!("{label} field ({NX}x{NY}), sampled with stride: ({stride_x}, {stride_y})");
        for y in (0..NY).rev().step_by(stride_y) {
            for x in (0..NX).step_by(stride_x) {
                print!("{:4.1} ", field.at(x, y));
            }
            println!();
        }
    }
}
