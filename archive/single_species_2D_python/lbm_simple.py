"""
MRT runner for the Maxwell–Stefan LBM example.

What this script does
---------------------
- Builds the three-species 1D/2D slab setup used in the Fortran reference.
- Initialises partial pressures with a tanh profile, converts them to
  equilibria (zero velocity), and runs the θ-splitting + MRT Maxwell–Stefan
  step in a loop to demonstrate diffusion and coupling.

Key routines
------------
- `initialise_pressures`: reproduces MIXLBM’s hyperbolic-tangent initial data.
- `initialise_from_partial_pressures`: loads feq(ρσ, 0) into the distributions.
- `maxwell_stefan_step`: performs one θ-split step with MRT collisions; see
  comments in `functions/maxwell_stefan.py` for detailed per-step notes.
"""

import numpy as np

from functions.maxwell_stefan import (
    MaxwellStefanState,
    initialise_from_partial_pressures,
    maxwell_stefan_step,
)


def initialise_pressures(nx: int, ny: int, num_species: int, molecular_weights: np.ndarray) -> np.ndarray:
    """
    Reproduce the hyperbolic tangent initialisation from MIXLBM.f90.

    This prescribes species partial pressures pσ(x) which the solver converts
    to densities via ρσ = 3 pσ / φσ (see `initialise_from_partial_pressures`).
    """
    psigma = np.zeros((num_species, nx, ny), dtype=np.float64)
    for i in range(nx):
        profile = np.tanh(i - nx / 2.0)
        for j in range(ny):
            psigma[0, i, j] = 0.319 * (1.0 - profile) / 2.0 + 1e-4
            psigma[1, i, j] = 0.528 * (1.0 - profile) / 2.0 + 1e-4
            psigma[2, i, j] = 0.847 * (1.0 + profile) / 2.0 + 0.153
    return psigma


def main() -> None:
    species = 3
    nx = 60
    ny = 3
    iterations = 10_000
    output_stride = 1_000

    molecular_weights = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    phi = 1.0 / molecular_weights

    # Initial partial pressures pσ for each species.
    psigma_init = initialise_pressures(nx, ny, species, molecular_weights)

    # Load feq(ρσ, 0) into populations and precompute macroscopic fields.
    state = MaxwellStefanState(species, nx, ny)
    initialise_from_partial_pressures(state, phi, psigma_init)

    for step in range(iterations + 1):
        if step % output_stride == 0:
            rsigma = state.rsigma
            rho_means = np.mean(rsigma, axis=(1, 2))
            print(f"step {step:6d} mean densities: " + ", ".join(f"{value:.6f}" for value in rho_means))

        if step == iterations:
            break

        # One θ-splitting step with MRT collision and MS coupling.
        maxwell_stefan_step(
            state,
            phi=phi,
            molecular_weights=molecular_weights,
            nN=20,
            nB=15,
            theta=0.5,
        )


if __name__ == "__main__":
    main()
