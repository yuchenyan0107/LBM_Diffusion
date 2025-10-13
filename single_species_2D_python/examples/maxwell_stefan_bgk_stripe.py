"""Simple Maxwell–Stefan BGK example with a drifting binary mixture.

The setup matches the user's request:

- Domain: rectangular lattice with periodic boundaries along ``x`` and
  bounce-back (no-slip) boundaries along ``y`` (enforced inside the BGK step).
- Components: two species labelled ``A`` and ``B`` with different molecular
  weights.
- Initial state: the lattice is filled with ``A`` except for a narrow vertical
  stripe rich in ``B`` in the centre of the domain.
- Flow: a weak bulk velocity is imposed along ``x`` to illustrate the
  competition between advection and Maxwell–Stefan diffusion.

Running the script produces a PNG image with the final concentrations for both
species using different colour maps.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

from functions.maxwell_stefan import (
    MaxwellStefanState,
    equilibrium_distribution,
    initialise_from_partial_pressures,
)
from functions.maxwell_stefan_bgk import maxwell_stefan_bgk_step


@dataclass
class SimulationConfig:
    """Container for the physical and numerical parameters of the run."""

    nx: int = 60
    ny: int = 24
    timesteps: int = 160
    output_stride: int = 20
    stripe_width: int = 6
    base_pressure_a: float = 0.95
    base_pressure_b: float = 0.05
    stripe_pressure_a: float = 0.05
    stripe_pressure_b: float = 0.95
    molecular_weights: tuple[float, float] = (28.0, 32.0)
    initial_ux: float = 0.02
    initial_uy: float = 0.0
    theta: float = 0.5
    nN: int = 20
    nB: int = 15


def build_partial_pressures(config: SimulationConfig) -> np.ndarray:
    """Create the initial partial pressure field for species A and B."""

    psigma = np.empty((2, config.nx, config.ny), dtype=np.float64)
    psigma[0, :, :] = config.base_pressure_a
    psigma[1, :, :] = config.base_pressure_b

    stripe_start = config.nx // 2 - config.stripe_width // 2
    stripe_end = stripe_start + config.stripe_width
    psigma[0, stripe_start:stripe_end, :] = config.stripe_pressure_a
    psigma[1, stripe_start:stripe_end, :] = config.stripe_pressure_b
    return psigma


def impose_initial_velocity(
    state: MaxwellStefanState,
    phi: np.ndarray,
    ux: float,
    uy: float,
) -> None:
    """Overwrite the populations with an equilibrium carrying the desired flow."""

    num_species, nx, ny = state.rsigma.shape
    for s in range(num_species):
        for i in range(nx):
            for j in range(ny):
                rho = state.rsigma[s, i, j]
                feq = equilibrium_distribution(rho, phi[s], ux, uy)
                state.fd[s, :, i, j] = feq
                state.fd_new[s, :, i, j] = feq
                state.uxsigma[s, i, j] = ux
                state.uysigma[s, i, j] = uy


def run_simulation(config: SimulationConfig) -> tuple[np.ndarray, np.ndarray]:
    """Run the Maxwell–Stefan BGK solver and return the final densities."""

    molecular_weights = np.array(config.molecular_weights, dtype=np.float64)
    phi = 1.0 / molecular_weights

    psigma_init = build_partial_pressures(config)
    state = MaxwellStefanState(num_species=2, nx=config.nx, ny=config.ny)
    initialise_from_partial_pressures(state, phi, psigma_init)
    impose_initial_velocity(state, phi, config.initial_ux, config.initial_uy)

    for step in range(config.timesteps + 1):
        if step % config.output_stride == 0:
            mean_a = np.mean(state.rsigma[0])
            mean_b = np.mean(state.rsigma[1])
            print(
                f"step {step:5d}: mean density A={mean_a:.6f}, B={mean_b:.6f}"
            )

        if step == config.timesteps:
            break

        maxwell_stefan_bgk_step(
            state,
            phi=phi,
            molecular_weights=molecular_weights,
            nN=config.nN,
            nB=config.nB,
            theta=config.theta,
        )

    return state.rsigma.copy(), phi


def plot_concentrations(rsigma: np.ndarray, output_path: Path) -> None:
    """Save a two-panel figure showing the densities for both species."""

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True, sharey=True)
    species_labels = ["Species A", "Species B"]
    cmaps = ["Blues", "Reds"]

    for idx, ax in enumerate(axes):
        im = ax.imshow(
            rsigma[idx].T,
            origin="lower",
            cmap=cmaps[idx],
            aspect="auto",
        )
        ax.set_title(f"{species_labels[idx]} density")
        ax.set_xlabel("x")
        if idx == 0:
            ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Maxwell–Stefan BGK: advective stripe mixing")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved concentration map to {output_path}")


def main() -> None:
    config = SimulationConfig()
    rsigma, _ = run_simulation(config)

    output_dir = Path(__file__).resolve().parent
    plot_concentrations(rsigma, output_dir / "stripe_concentrations.png")


if __name__ == "__main__":
    main()
