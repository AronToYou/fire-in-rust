- [x] replace `Vec` with `[]` ✅ 2025-09-14
	- [x] `Box<[[f32; NY]; NX]>` ✅ 2025-09-14
- [x] P convertible between `usize` & `f32` ✅ 2025-09-15
- [x] Rust Obsidian Vault within repo ✅ 2026-01-03
	- [x] README linked to main page in book ✅ 2026-01-03
- [x] Write down equations for ghost ✅ 2026-08-22
- [ ] print rows
- [ ] Update README to match Extended README
- [ ] write benchmarks
- [ ] write unit tests
- [ ] determine why `NaN` values occur
	- [ ] parallelize the `NaN` check
- [ ] [documentation](https://doc.rust-lang.org/rust-by-example/meta/doc.html)
- [ ] [enums](https://stackoverflow.com/questions/28028854/how-do-i-match-enum-values-with-an-integer)
- [ ] switch to face centered velocity fields
	- [ ] fix boundary conditions

## Directions
- Field trait
	- apply trait to all field objects or special Field type
	- bilinear sample implemented for trait/type
- Where to define min
	- in clamp
	- in `P<T>`
- Where to define clamp
	- in `Sim` method
	- in closure field

## Other

- [ ] trait of indexing
- [ ] trait of printing
- [x] const clamp impl ✅ 2025-09-14
- [ ] stack too much Sim

- [ ] generalize to 3D
- [ ] generic Field type
- [ ] somehow include physics package
- [ ] add smoke density
	- [ ] density dependent [gravity force](https://www.researchgate.net/publication/2390581_Visual_Simulation_of_Smoke)