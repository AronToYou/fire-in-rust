- only do first 2 steps of Stam's 4-step loop
	- instead solve Poisson for pressure like in `fire.pdf`

- [ ] Fuel should be injected from solids using normal velocity:
$$V_f = V_{solid} + (\rho_{solid}/\rho_f - 1)S$$

This implementation utilizes flat 1D vectors instead of nested vectors (`Vec<Vec<f64>>`) to guarantee memory locality and cache efficiency, which is critical for high-performance scientific computing in Rust.

I don't understand why you say that $p_{i,j}$ is an unknown in the five-point-difference approximation of the Poisson $\nabla^2p$.


Page 2 last piece of text of ns.pdf explains why the pressure is not required.
- also why Poisson on velocity
	- followed by projection
	- results in the time derivative of u



## Physically Based Fire Simulation solves
$$\vec{u}_t = -(\vec{u}\cdot\nabla)\vec{u} - \frac{1}{\rho}\nabla p + \vec{f}$$
- Solves **Pressure Poisson** explicitly

## Jos Stam solves
$$\vec{u}_t = \mathbf{P}(-(\vec{u}\cdot\nabla)\vec{u} + \nu\nabla^2\vec{u} + \vec{f})$$
- includes **diffusion**
- avoids pressure term, instead applies **Projection** step
$$\nabla^2q = \nabla\cdot\vec{w}_3$$
$$\vec{w}_4 = \vec{w}_3 - \nabla q$$