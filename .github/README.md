<a href="https://codeberg.org/ryarza/common-envelope-wind-tunnel">
    <img alt="Get it on Codeberg" src="https://get-it-on.codeberg.org/get-it-on-white-on-black.svg" height="60">
</a>

This repository mirrors the <a href=https://codeberg.org/ryarza/common-envelope-wind-tunnel>original repository</a> hosted on <a href=https://codeberg.org/>Codeberg</a>, a <a href="https://www.gnu.org/philosophy/free-sw.html">free</a> and nonprofit platform.

---

Setup for the FLASH hydrodynamics code to simulate flow around objects embedded in stellar envelopes. Written by Morgan MacLeod ([@morganemacleod](https://github.com/morganemacleod)) and Andrea Antoni ([@andreaantoni](https://github.com/andreaantoni)).

## Installation
This repository contains two folders:
- `cewt`: The common envelope wind tunnel setup. This folder should go inside `source/Simulation/SimulationMain/`.
- `cewindgrav`: A physics module dependency of the main setup. This folder should go inside `source/physics/Gravity/GravityMain/`.

## Example FLASH setup line
```
./setup cewt -3d -auto -maxblocks=1000 -objdir=obj_cewt +pm4dev +uhd
    cd obj_cewt; make
```

## Important flash.par parameters
Simulations are mostly characterized by (see publications below for definitions)
- `sim_q`: mass ratio between the embedded object and the mass enclosed by its orbit.
- `sim_epsilon_grad`: dimensionless density gradient (number of density scale heights within an accretion radius).
- `sim_radius`: ratio between geometrical and gravitational radii for the embedded object.
- `sim_bctype`: whether the surface of the object behaves as a reflective (`sim_bctype=1`) or sink (`sim_bctype=2`) boundary. The reflective boundary requires the unsplit solver (`+uhd` in the setup line).

## Citation
If you use this setup please cite the following publications:
- [MacLeod et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...803...41M)
- [MacLeod et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...838...56M)
- [Yarza et al. 2022](https://arxiv.org/abs/2203.11227)
