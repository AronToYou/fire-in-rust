- $\phi$ as a signed distance function
	- distance to closest point on the surface (the surface is a list of points)
	- $k$ width of narrow band
		- $\pm6$ from front is good balance
		- proportional to $F$?
		- based on values of $\phi$, add new points when necessary           
	- once zero-band reaches Narrow Band boundary
		- stop and rebuild
- only do first 2 steps of Stam's 4-step loop
	- instead solve Poisson for pressure like in `fire.pdf`
- Page `140` of *Level set methods and fast marching methods* points towards
	- extension velocity field to maintain signed distance function $|\nabla\phi|=1$

Sussman, M., Smereka, P. and Osher, S.J., A Level Set Method for
Computing Solutions to Incompressible Two-Phase Flow, J. Comp. Phys.
114, pp. 146-159, 1994.

Chang, Y.C., Hou, T.Y., Merriman, B., and Osher, S.J., A Level Set
Formulation of Eulerian Interface Capturing Methods for Incompressible
Fluid Flows, Jour. Comp. Phys., 124, pp. 449-464, 1996.
## CFL( Courant-Friedrichs-Lewy) Condition
$$\pdv{w}{t} + u\pdv{w}{x} = 0$$
$$C = \frac{\Delta t}{h}\leq C_{max} \approx  1$$
## 18 Combustion (page 240/264)
- $S_f$ depends on local curvature! (page 244/269)
	- perturbation in closer contact with hot gas -> faster heating -> faster reaction
	- $S_f^o$ - speed of planar flame
	- $\kappa$- local curvature
	- $L$ constant
$$S_f = S_f^o(1 - L\kappa)$$
- hydrodynamic advection
- reaction
- [ ] Fuel should be injected from solids using normal velocity (accounts for expansion):
$$V_f = V_{solid} + (\rho_{solid}/\rho_f - 1)S$$

This implementation utilizes flat 1D vectors instead of nested vectors (`Vec<Vec<f64>>`) to guarantee memory locality and cache efficiency, which is critical for high-performance scientific computing in Rust.

I don't understand why you say that $p_{i,j}$ is an unknown in the five-point-difference approximation of the Poisson $\nabla^2p$.


Page 2 last piece of text of ns.pdf explains why the pressure is not required.
- also why Poisson on velocity
	- followed by projection
	- results in the time derivative of u

# Parameters
---
$c_T \approx 3000K/s$ - cooling constant
$\alpha \approx 0.15 m/(Ks^2)$


# Fast Marching Level Set Method
---
### A Fast Marching Level Set Method for Monotonically Advancing Fronts
- Due to upwind differencing, information propagates from smaller to larger values of $T$
- hence choosing smallest value of $T$ ensures the most accurate value of $T$
1. Initialize
	1)  (Alive points: shaded points): Let $A$ be the set of all grid points $\{i,j=1\}$; 
		- set $T_{i,1} = 0.0$ for all points in $A$.
	2) (Narrow band points: circles): Let *Narrow Band* be the set of all grid points $\{i,j=2\}$;
		- set $T_{i,2} = dy/F_{ij}$ for all points in *Narrow Band*.
	3) (Far away points: rectangles): Let *Far Away* be the set of all grid points $\{i,j>2\}$;
		- set $T_{ij} = \infty$ for all points in *Far Away*
2. Marching Forward
	1) Begin Loop: Let $(i_{min},j_{min})$ be the point in *Narrow Band* with the smallest value for $T$.
	2) Add the point $(i_{min},j_{min})$ to $A$;
		- remove it from *Narrow Band*.
	3) Tag as neighbors points $(i_{min}\pm1,j_{min})$ / $(i_{min},j_{min}\pm1)$ in *Narrow Band* or *Far Away*.
		- If a neighbor is in *Far Away*, remove it from that list and add it to *Narrow Band*.
	4) Recompute $T$ at all neighbors according to Eq. 8, selecting the largest possible solution to the quadratic equation.
	5) Return to top of Loop.

$$\phi_t + F(x,y,z)|\nabla\phi| = 0$$
$$\phi(\vec{x},t=0) = \Gamma$$

In this framework, start with
- front at $y=0$
- $F(x,y) > 0$ and periodic in $x$
Upwind gradient is approximated as (**Eq. 8**):
$$\big[\max(\max(D_{ij}^{-x}T,\:0),\;-\min(D_{ij}^{+x}T,\:0))^2 + \max(\max(D_{ij}^{-y}T,\:0),\;-\min(D_{ij}^{+y}T,\:0))^2\big] = 1/F_{ij}^2$$
where
$$D_{ij}^{+x}\phi = (\phi_{i+1,\ j} - \phi_{i,j})/(\Delta x)$$
$$D_{ij}^{-x}\phi = (\phi_{i,j} - \phi_{i-1,\ j})/(\Delta x)$$
### Level set methods and fast marching methods (book)
1) Tag *Alive* points in narrow band.
2) Build *Land Mines* to indicate near edge.
3) Initialize *Far Away* points outside the narrow band with large positive (negative) values if outside (inside) the front.
4) Solve level set equation until *Land Mine* hit
5) Rebuild and loop
# Differences
---
## Physically Based Fire Simulation solves

$$\vec{u}_t = -(\vec{u}\cdot\nabla)\vec{u} - \frac{1}{\rho}\nabla p + \vec{f}$$
- Solves **Pressure Poisson** explicitly

## Jos Stam solves
$$\vec{u}_t = \mathbf{P}(-(\vec{u}\cdot\nabla)\vec{u} + \nu\nabla^2\vec{u} + \vec{f})$$
- includes **diffusion**
- avoids pressure term, instead applies **Projection** step
$$\nabla^2q = \nabla\cdot\vec{w}_3$$
$$\vec{w}_4 = \vec{w}_3 - \nabla q$$