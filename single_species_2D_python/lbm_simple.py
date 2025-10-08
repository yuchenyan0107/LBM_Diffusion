import numpy as np

from functions.maxwell_stefan import (
    MaxwellStefanState,
    initialise_from_partial_pressures,
    maxwell_stefan_step,
)


def initialise_pressures(nx: int, ny: int, num_species: int, molecular_weights: np.ndarray) -> np.ndarray:
    """
    Reproduce the hyperbolic tangent initialisation from MIXLBM.f90.
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

    psigma_init = initialise_pressures(nx, ny, species, molecular_weights)

    state = MaxwellStefanState(species, nx, ny)
    initialise_from_partial_pressures(state, phi, psigma_init)

    for step in range(iterations + 1):
        if step % output_stride == 0:
            rsigma = state.rsigma
            rho_means = np.mean(rsigma, axis=(1, 2))
            print(f"step {step:6d} mean densities: " + ", ".join(f"{value:.6f}" for value in rho_means))

        if step == iterations:
            break

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

