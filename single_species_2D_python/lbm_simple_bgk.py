"""
BGK runner for the Maxwell–Stefan LBM example.

This is identical to `lbm_simple.py` in setup but calls the BGK-based step
(`maxwell_stefan_bgk_step`) instead of the MRT one, illustrating that the
Maxwell–Stefan θ-splitting works without MRT.
"""

import numpy as np

from functions.maxwell_stefan import (
    MaxwellStefanState,
    initialise_from_partial_pressures,
)
from functions.maxwell_stefan_bgk import maxwell_stefan_bgk_step


def initialise_pressures(nx: int, ny: int, num_species: int, molecular_weights: np.ndarray) -> np.ndarray:
    """
    Reproduce the hyperbolic tangent initialisation from MIXLBM.f90.

    Same as the MRT runner; the solver converts these to densities via
    ρσ = 3 pσ / φσ at initialisation.
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
            print(f"[BGK] step {step:6d} mean densities: " + ", ".join(f"{value:.6f}" for value in rho_means))

        if step == iterations:
            break

        # One θ-splitting step with BGK collision and MS coupling.
        maxwell_stefan_bgk_step(
            state,
            phi=phi,
            molecular_weights=molecular_weights,
            nN=20,
            nB=15,
            theta=0.5,
        )


if __name__ == "__main__":
    main()
