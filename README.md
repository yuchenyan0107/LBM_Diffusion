# LBM_Diffusion
A Lattice-Boltzmann Method tool to simulate diffusion-advection of multi-species systems.

The algorithm is largely inspired by Musubi, an open-source LBM project.
https://geb.inf.tu-dresden.de/doxy/musubi/page/features/multispecies.html

Zudrop, J., Masilamani, K., Roller, S., & Asinari, P. (2017). A robust lattice Boltzmann method for parallel simulations of multicomponent flows in complex geometries. Computers & Fluids, 153, 20–33. https://doi.org/10.1016/j.compfluid.2017.04.021

The implementation of the Dirichlet and Robin boundary conditions are based on:
Zhang, T., Shi, B., Guo, Z., Chai, Z., & Lu, J. (2012). General bounce-back scheme for concentration boundary condition in the lattice-Boltzmann method. Physical Review E, 85(1), 016701. https://doi.org/10.1103/physreve.85.016701

## 2D_LBM_MS file structure

The `2D_LBM_MS` folder contains the 2D Maxwell-Stefan LBM implementation, example notebooks, and sample outputs/results:

```text
2D_LBM_MS/
  sim_functions/               # Core simulation package
    MS_LBM_functions.py        # Main simulation loop and workflow functions
    eq_and_ms.py               # Equilibrium distribution and Maxwell-Stefan terms
    boundary_conditions.py     # Boundary-condition handlers
    BC_2.py                    # Additional/alternative BC implementations
    parameters.py              # Parameter definitions and setup helpers
    plotting_functions.py      # Plotting and visualization utilities
    common.py                  # Shared utility functions
    __init__.py                # Package entry

  frame_output/                # Generated frame images from simulations

  *.ipynb                      # Case-study and demo notebooks
                               # (e.g., shear decay, stripe diffusion, thick absorption)
  *.npy                        # Saved simulation states/results for reuse and analysis
  MS_LBM_Algorithm_Notes.md    # Notes on model/algorithm details
```

### Notebook quick guide (`2D_LBM_MS`)

- `parameters.ipynb`: Converts physical inputs to LBM-scale parameters and reports key dimensionless numbers (e.g., `Re`, `Pe`).
- `shear_decay.ipynb`: Runs a laminar shear-wave decay test and fits decay to verify effective viscosity.
- `stripe_diffusion.ipynb`: Simulates diffusion of a narrow/Gaussian stripe and checks diffusion scaling over time.
- `thick_absorption_demo.ipynb`: Small 2-species demo of flow + bottom absorption boundary behavior.
- `thick_absorption_4_species.ipynb`: Main 4-species thick-absorption case with flow, concentration frames, and saved end state.
- `PS_thick_4_species.ipynb`: Parameter-sweep version of reaction chamber scale flow
- `thick_concentration_analysis.ipynb`: Post-processes saved simulation arrays (e.g., `f_simulation.npy`) to analyze concentration profiles.
- `sample_mask1.ipynb`: 4-species masked-geometry simulation example and result export (`f_sample_simulation.npy`).
- `PS_mask.ipynb`: Parameter-sweep version of the masked 4-species case (sweeps absorption coefficient and saves result lists).


Here we give some examples in jupyter notebooks to demonstrate the current implementation:

1. (Almost) pure diffusion in LBM: 
we compared the result from the pure diffusion model ($\nabla^2c=0$) which people used to study selective area growth (SAG) in MOCVD, to the diffusion in LBM.
The shape and the thickness of these two models matches in general, but with some differences.
With this demo, we showed that the Dirichlet and Robin boundary conditions are correctly modeled

2. The laminar shear wave: we used the decay in flow speed to calculate the kinematic viscosity, and confirmed that it matches the theoretical value.
$\nu = c_{s}^{2}(\frac{1}{rel}-\frac{1}{2})$, where the effective relaxation $rel = \frac{\lambda}{1+\theta\lambda}$.

3. Flow around cylinder: we tried to demonstrate the turbulence phenomena with the flow around cylinder. 
The result (average flow field and Reynolds stresses) was compared with the literature. 
Nguyen, Q. D., & Lei, C. (2022). A PIV study of blockage ratio effects on flow over a confined circular cylinder at low Reynolds numbers. Experiments in Fluids, 64(1). https://doi.org/10.1007/s00348-022-03548-w
(Note: the results didn't exactly match, we are still looking into it)

4. A simple case of the diffusion of a narrow stripe. We verified the diffusion follows the square root of time.

5. Tri-species diffusion demonstrating the Maxwell-Stefan effect. 
A chamber is initialized with uniform 50% nitrogen, 50% carbon dioxide in the left half, and 50% hydrogen in the right half.
The result shows that despite the nitrogen concentration is driven by the diffusion of other species, because of differences in the Maxwell-Stefan coefficients.
(We will try to find method to validate this result quantitatively)
