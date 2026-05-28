# 3D AA6061-T6 Torsion Beam — Linear Elasticity with AMR

<p align="center">
  <img src="https://img.shields.io/badge/FreeFEM++-Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/AA6061--T6-Aluminium%20Alloy-silver?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/3D-Linear%20Elasticity-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Torque-500%20Nm-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/AMR-Von%20Mises%20Driven-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/ParaView-VTK%20Export-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-MIT-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A 3D finite element simulation of <b>linear elastic torsion</b> of an
  AA6061-T6 aluminium alloy beam using FreeFEM++.
  A torque of <b>500 Nm</b> is applied to the free end while the other end is fully clamped.
  <b>Adaptive mesh refinement (AMR)</b> driven by the von Mises stress gradient
  automatically concentrates elements at the clamp corners and loaded face.
  All fields are exported as <i>VTU/PVD files</i> for smooth 30-frame ParaView animation.
</p>

<img width="1008" height="772" alt="torque" src="https://github.com/user-attachments/assets/ac1354f1-43b4-4fa6-ac2d-f0e82cebf4a4" />

---

## Physics

The simulation solves the 3D static linear elasticity equations (Lamé system), capturing
the full torsional deformation and stress distribution in the beam. Key features include:

- 3D static linear elasticity solved with a single fully coupled matrix solve
- Taylor-Hood P1 × P1 × P1 element triple — stable for 3D solid mechanics
- Traction boundary condition derived from correct polar moment formula: `C = T / Ip`
- Adaptive mesh refinement (AMR) every 8 steps driven by von Mises stress gradient
- Fine cells concentrate at clamp corners (stress concentration); coarse cells fill beam interior
- UMFPACK direct sparse solver via `varf` + `matrix^-1` — official FreeFEM module pattern
- VTU/PVD export for ParaView animation with 5 output fields per frame
- Deformation baked into saved mesh at 200× scale — no Warp filter needed in ParaView

---

## Geometry

```
        z
        |
        |      ← torsion face (label 2, x=Lx)
        |   /
        |  /
        | /_________ y
       /
      /  ← clamp face (label 4, x=0)
     x

  Beam: 200 mm × 40 mm × 40 mm (Lx × Ly × Lz)
  Clamp:   left face  x = 0      (all DOF fixed)
  Torsion: right face x = 0.200  (distributed traction)
  Free:    top, bottom, front, back (natural Neumann)
```

- Boundary labels (from `buildlayers`):
  - Label 0: Top + bottom z-faces — free (natural Neumann)
  - Label 1: Front face (y=0) — free
  - Label 2: Right face (x=Lx) — torsion traction applied here
  - Label 3: Back face (y=Ly) — free
  - Label 4: Left face (x=0) — fully clamped (ux=uy=uz=0)

---

## Material Properties — AA6061-T6

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Young's modulus | E | 69.0 | GPa |
| Poisson ratio | ν | 0.33 | — |
| Shear modulus | G = E/2(1+ν) | 25.94 | GPa |
| Density | ρ | 2700 | kg/m³ |
| Yield stress | σ_y | 276 | MPa |
| Lame λ | λ | 51.13 | GPa |
| Lame μ | μ | 25.94 | GPa |

---

## Process Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Max torque | T_max | 500 Nm | Near-yield torque for AA6061 |
| Load steps | nStep | 30 | Ramped 0 → T_max |
| Traction coefficient | C_max = T/Ip | 1.172×10⁹ N/m² | Correct polar moment formula |
| Polar moment | Ip = Ly·Lz·(Ly²+Lz²)/12 | 4.267×10⁻⁷ m⁴ | Square cross-section |
| Warp scale | — | 200× | Visual scale baked into saved mesh |
| AMR frequency | adaptEvery | 8 steps | Adapt every 8 load steps |

---

## Mesh

| Parameter | Value | Description |
|---|---|---|
| 2D base mesh | B1(15)+B2(6)+B3(15)+B4(6) | Unstructured triangular in x-y plane |
| Layers | 6 | Extruded in z-direction |
| Initial vertices | ~1,100 | Before first AMR pass |
| Max 2D vertices | 8,000 | Hard cap (nbvxAMR) |
| Method | buildlayers | Prismatic extrusion of 2D mesh |

### AMR Settings

| Setting | Value | Effect |
|---|---|---|
| errAMR | 0.05 | Target interpolation error |
| hminAMR | 0.003 m | Minimum element size (3mm) |
| hmaxAMR | 0.020 m | Maximum element size (20mm) |
| nbvxAMR | 8,000 | Maximum 2D vertex count |
| Refine on | VonMises at z=Lz/2 | Stress-driven isotropic refinement |

Fine cells concentrate at the clamp corners (x=0) where stress concentration is highest
and at the torsion face (x=Lx). Coarse cells fill the beam interior and far from boundaries.

---

## Finite Element Spaces

| Space | Element | Fields |
|---|---|---|
| Uh | P1 × P1 × P1 | Displacement (ux, uy, uz) |
| Sh | P1 | All derived nodal scalar fields |
| Ph | P0 | Element-centred von Mises (intermediate) |

P1 elements used throughout — avoids the vector-component assignment errors
common with product spaces in FreeFEM 3D. Official FreeFEM elasticity module pattern.

---

## Governing Equations

### Linear Elasticity — Lamé System

```
−∇·σ = f            in Ω

σ_ij = λ δ_ij ∇·u + 2μ ε_ij

ε_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i)
```

### Variational Form

```
∫_Ω [λ (∇·u)(∇·v) + 2μ ε(u):ε(v)] dV = ∫_Γ_torsion t·v dA
```

### Torsion Traction (correct formula)

```
t = (0,  −C·(z − Lz/2),  +C·(y − Ly/2))

Torque: T = ∫∫ [(y−Ly/2)·tz − (z−Lz/2)·ty] dA
          = C · Ip

Polar moment: Ip = Ly·Lz·(Ly² + Lz²)/12

Therefore: C = T / Ip
```

### Boundary Conditions

```
Clamp   (label 4, x=0):   ux = uy = uz = 0     (Dirichlet, fully fixed)
Torsion (label 2, x=Lx):  t = (0, −C(z−Lz/2), +C(y−Ly/2))  (Neumann)
Free    (labels 0,1,3):    natural Neumann (zero traction)
```

---

## Numerical Method

| Aspect | Choice |
|---|---|
| Spatial discretisation | Finite Element Method (FEM) |
| Element type | P1 × P1 × P1 (linear tetrahedral) |
| Assembly | `varf` + `matrix^-1` (official FreeFEM module pattern) |
| Linear solver | UMFPACK direct sparse solver |
| Mesh construction | `buildlayers` (prismatic extrusion) |
| Mesh adaptation | `adaptmesh()` — error-driven isotropic refinement |
| Loading | Incremental: 30 steps, C ramped 0 → C_max |
| Deformation output | `movemesh3` with 200× scale baked in |

---

## Output Fields (ParaView)

Each `.vtu` frame contains the following fields:

| Field | Type | Description | Units |
|---|---|---|---|
| DispMag | P1 nodal | Displacement magnitude √(ux²+uy²+uz²) | m |
| VonMises | P1 nodal | Von Mises stress (smooth nodal projection) | Pa |
| Ux | P1 nodal | X-displacement component | m |
| Uy | P1 nodal | Y-displacement component (shows twist) | m |
| Uz | P1 nodal | Z-displacement component (shows twist) | m |

---

## Expected Physics

```
Step  1/30 :  T =  16.7 Nm   umax = 0.010 mm   VM_max =  9.2 MPa   twist = 0.020°
Step 10/30 :  T = 166.7 Nm   umax = 0.101 mm   VM_max = 91.9 MPa   twist = 0.204°
Step 20/30 :  T = 333.3 Nm   umax = 0.202 mm   VM_max = 183.8 MPa  twist = 0.407°
Step 30/30 :  T = 500.0 Nm   umax = 0.303 mm   VM_max = 275.7 MPa  twist = 0.611°

Analytical twist angle (Roark): θ = T·L / GJ = 0.611°  ✓
GJ = G × 0.1406 × a⁴ = 25.94e9 × 3.599e-7 = 9337 N·m²
```

---

## Repository Structure

```
aa6061_torsion.edp                     # Main FreeFEM++ simulation script
README.md                              # This file

D:\freefem++\aa6061_v4\
├── beam.pvd                           # Master animation (open this in ParaView)
├── frame000.vtu                       # T = 16.7 Nm  (step 1)
├── frame001.vtu                       # T = 33.3 Nm  (step 2)
├── ...
└── frame029.vtu                       # T = 500.0 Nm (step 30 — full torque)
```

---

## How to Run

### Requirements

- FreeFEM++ v4.10 or later: https://freefem.org
- ParaView v5.x or later: https://www.paraview.org

### Step 1 — Run the simulation

```bash
FreeFem++ aa6061_torsion.edp
```

The script will:
1. Build initial 3D prismatic mesh (~1,100 vertices, 6 layers)
2. Run a quick pre-solve at full load to determine global colour range
3. For each of 30 load steps: assemble stiffness matrix, solve K·u = F, compute post-processing
4. Apply AMR every 8 steps driven by von Mises stress at beam mid-plane
5. Save deformed mesh (200× warp) + all 5 fields as VTU per step
6. Write complete PVD collection file at the end

Console output example:
```
Step  1/30  T=16.7Nm   umax=0.010mm  twist=0.020deg  VM=9.2MPa   nv=1161
Step  8/30  T=133.3Nm  umax=0.081mm  twist=0.163deg  VM=73.6MPa  nv=1161
  AMR: nv=4800
Step  9/30  T=150.0Nm  umax=0.091mm  twist=0.183deg  VM=82.8MPa  nv=4800
...
Step 30/30  T=500.0Nm  umax=0.303mm  twist=0.611deg  VM=275.7MPa nv=4800
```

### Step 2 — Open in ParaView

1. `File > Open` → navigate to `D:\freefem++\aa6061_v4\`
2. Select `beam.pvd` → OK
3. Choose PVD Reader when prompted → OK
4. Click `Apply`
5. Set colour field to `VonMises`
6. Click `Rescale to Data Range Over All Timesteps` for consistent colour scale

### Step 3 — Visualise the deformation fields

**Option A — Von Mises stress (recommended first view)**
```
Color by VonMises
Colormap: Rainbow, range 0 to 276e6 Pa (yield stress)
Press Play
→ Stress grows from zero at clamp (fixed end)
→ Peak stress at clamp corners — classic torsion stress concentration
→ Uniform shear stress in beam interior at full load
```

**Option B — Twist displacement**
```
Color by Uz  (or Uy)
Colormap: Cool-Warm (diverging), centred at zero
→ Top face moves in −z (blue), bottom face moves in +z (red)
→ Left face (clamp) stays zero
→ Right face shows maximum rotation
→ Progressive twist builds from frame 0 to frame 29
```

**Option C — Displacement magnitude**
```
Color by DispMag
Colormap: Blue → Red
→ Zero at clamp (blue), maximum at free end corners (red)
→ Radial distribution: corners move more than face centres
```

**Option D — AMR mesh evolution**
```
Representation: Wireframe
→ Coarse uniform mesh in early frames
→ After frame 8: fine cells appear at clamp corners
→ After frame 16: further refinement tracks stress concentration
```

Press `Play` and set animation speed to **Slowest** for smooth colour gradient transition.

---

## What to Look for in Results

### Torsion Stress Distribution
At full load (step 30), the von Mises stress field shows the classical Saint-Venant
torsion pattern: zero stress at the beam centroid axis, maximum shear at the face
centres, and slightly reduced stress at the corners. The clamp face (x=0) shows the
highest overall stress due to the warping constraint.

### Progressive Twist
Watching `Uz` or `Uy` across 30 frames shows the twist angle building linearly
with load — expected for linear elasticity. The free end rotates while the clamped end
remains perfectly fixed. The total twist at full load is 0.611° as predicted analytically.

### AMR Stress Concentration
The `Wireframe` view after step 8 shows the mesh refining at the clamp corners
exactly where the stress gradient is steepest. This improves accuracy of the peak
stress prediction without increasing global mesh density.

### Volume Conservation
The beam cross-section at the free end should remain approximately square (small
Poisson distortion at ν=0.33). Any significant distortion of the cross-section shape
indicates numerical issues with the mesh.

---

## Common Errors and Fixes

| Error | Cause | Fix |
|---|---|---|
| `umax = 1e+12` | Clamp BC not applied (wrong label) | Run label scan script to confirm label 4 is x=0 face |
| `out_of_memory` in UMFPACK | AMR mesh too large | Reduce `nbvxAMR` to 5000 or lower |
| No visible deformation | Warp scale too small | Increase `warpScale` to 500 or 1000 |
| `Error interpolation vectorial` | Self-assign after AMR | Remove re-interpolation — static problem needs none |
| PVD XML parse error | File written with append | Write entire PVD in one `ofstream` block at end |
| Stress = 0 everywhere | Wrong traction formula | Use `C = T / Ip` not `C = T / (0.5·Ly·Area)` |
| Twist < 0.01 degrees | `C` too small (wrong Ip) | Verify `Ip = Ly·Lz·(Ly²+Lz²)/12` |

---

## Extending the Model

| Extension | What to change |
|---|---|
| Non-linear material | Replace linear `σ = C:ε` with Neo-Hookean or Johnson-Cook |
| Higher torque | Increase `Tmax` beyond 500 Nm (check VM < yield before using linear) |
| Rectangular cross-section | Change `Ly ≠ Lz`; update `Ip` formula accordingly |
| Combined loading | Add bending traction to inlet face alongside torsion |
| Finer mesh | Increase `nbvxAMR` to 20,000 and `nz layers` to 10 |
| Dynamic analysis | Add inertia term `ρ/dt²·u` to variational form (Newmark) |
| Thermal loading | Add `(3λ+2μ)·αT·ΔT·∇·v` thermal strain term |
| Different alloy | Change E, nu, rho; update yield stress reference in README |
| Fatigue analysis | Export stress tensor components; post-process with Goodman criterion |

---

## Citation

If you use this simulation, please cite:

```bibtex
@software{mishra_2026_aa6061_torsion,
  author    = {Mishra, A.},
  title     = {3D AA6061-T6 Torsion Beam — Linear Elasticity with AMR},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20428327},
  url       = {https://doi.org/10.5281/zenodo.20428327}
}
```

Plain text citation:

> Mishra, A. (2026). *3D AA6061-T6 Torsion Beam — Linear Elasticity with AMR*. Zenodo. https://doi.org/10.5281/zenodo.20428327

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20428327.svg)](https://doi.org/10.5281/zenodo.20428327)

---

## Author

**akshansh11**
GitHub: https://github.com/akshansh11

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
</p>

You are free to:

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

Under the following terms:

- **Attribution** — You must give appropriate credit and provide a link to this repository
- **NonCommercial** — You may not use the material for commercial purposes

Copyright 2026 akshansh11. All rights reserved for commercial use.
