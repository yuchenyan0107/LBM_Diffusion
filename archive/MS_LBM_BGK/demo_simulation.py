import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from functions.MS_LBM_functions import (
    D2Q9_CX,
    D2Q9_CY,
    OPPOSITE,
    theta as THETA,
    calculate_moment,
    calculate_m_mix,
    calculate_CHI,
    calculate_lambda,
    calculate_u_star,
    equilibrium,
    calculate_g_dagger,
    post_stream_Chi_S,
    distribution_semi_implicit,
    solve_ms_fluxes,
)


# --- Streaming with bounce-back on top/bottom, periodic in x ---------------
def lattice_stream_bounce_y(f: np.ndarray) -> np.ndarray:
    """Stream populations with periodic x and no-slip walls at y boundaries.

    This mirrors the behavior in the reference BGK implementation where
    populations attempting to stream past the horizontal walls are reflected
    to their opposite direction at the same node. It is based on the existing
    lattice_stream but adjusts only the y-boundary handling.
    """
    num_species, q, nx, ny = f.shape
    f_streamed = np.zeros_like(f)

    for k in range(q):
        dx = int(D2Q9_CX[k])
        dy = int(D2Q9_CY[k])

        # Periodic shift in x using roll for simplicity
        fk = np.roll(f[:, k, :, :], shift=dx, axis=1)  # (N, nx, ny)

        if dy == 0:
            # No y movement; copy directly into channel k
            f_streamed[:, k, :, :] = fk
            continue

        if dy == 1:
            # Push interior northwards
            f_streamed[:, k, :, 1:ny] += fk[:, :, 0 : ny - 1]
            # Bounce-back at top boundary: reflect into opposite channel at same node
            f_streamed[:, OPPOSITE[k], :, ny - 1] += fk[:, :, ny - 1]
        elif dy == -1:
            # Push interior southwards
            f_streamed[:, k, :, 0 : ny - 1] += fk[:, :, 1:ny]
            # Bounce-back at bottom boundary
            f_streamed[:, OPPOSITE[k], :, 0] += fk[:, :, 0]
        else:
            # D2Q9 dy is in {-1,0,1}; this branch is defensive
            f_streamed[:, k, :, :] += fk

    return f_streamed


# --- Configuration ----------------------------------------------------------
@dataclass
class Config:
    nx: int = 60
    ny: int = 24
    steps: int = 240
    output_stride: int = 12  # "each dozens of frames"
    stripe_width: int = 6
    # Set A as center stripe, B elsewhere
    base_pressure_A: float = 0.05
    base_pressure_B: float = 0.95
    stripe_pressure_A: float = 0.95
    stripe_pressure_B: float = 0.05
    molecular_weights: tuple[float, float] = (28.0, 32.0)
    initial_ux: float = 0.02
    initial_uy: float = 0.0
    nB: int = 15
    theta: float = THETA
    frames_dir: str = "demo_frames/array_stripe"


def initialise_stripe(config: Config):
    """Build initial distributions for a vertical A stripe and B elsewhere.

    Returns
    -------
    f : np.ndarray
        Populations, shape (2, 9, nx, ny)
    phi : np.ndarray
        Equation-of-state factors, shape (2,)
    """
    nx, ny = config.nx, config.ny
    N = 2
    phi = 1.0 / np.array(config.molecular_weights, dtype=np.float64)  # EOS factors

    # Partial pressures per species
    pA = np.full((nx, ny), config.base_pressure_A, dtype=np.float64)
    pB = np.full((nx, ny), config.base_pressure_B, dtype=np.float64)
    s0 = nx // 2 - config.stripe_width // 2
    s1 = s0 + config.stripe_width
    pA[s0:s1, :] = config.stripe_pressure_A
    pB[s0:s1, :] = config.stripe_pressure_B

    # Convert partial pressures to densities: p_sigma = phi_sigma * rho_sigma / 3
    rhoA = 3.0 * pA / phi[0]
    rhoB = 3.0 * pB / phi[1]
    rho_s = np.stack([rhoA, rhoB], axis=0)  # (2, nx, ny)

    # Build initial equilibrium with a small rightward velocity
    f = np.zeros((N, 9, nx, ny), dtype=np.float64)
    ux_s = np.full((N, nx, ny), config.initial_ux, dtype=np.float64)
    uy_s = np.full((N, nx, ny), config.initial_uy, dtype=np.float64)
    feq = equilibrium(f, rho_s, phi, ux_s, uy_s)
    f[...] = feq
    return f, phi


def bgk_step(
    f: np.ndarray,
    phi: np.ndarray,
    molecular_weights: np.ndarray,
    nB: int,
    theta: float,
    stream_fn=lattice_stream_bounce_y,
) -> np.ndarray:
    """One BGK time step using the array-based MS_LBM functions.

    This follows the notebook pipeline and keeps boundary handling in the
    streaming stage with periodic x and bounce-back on horizontal walls.
    """
    # Moments and mixture fields
    rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(f, phi)
    m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weights)
    CHI_sc = calculate_CHI(m_mix, molecular_weights, nB)
    lambda_s = calculate_lambda(rho_mix, p_mix, molecular_weights, nB)

    # First half-step: velocities to u*
    ux_star_s, uy_star_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s)
    feq = equilibrium(f, rho_s, phi, ux_star_s, uy_star_s)
    g_dagger_s = calculate_g_dagger(f, feq, lambda_s)

    # Streaming with walls in y
    f_streamed = stream_fn(f)

    rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(f_streamed, phi)
    m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weights)
    CHI_sc = calculate_CHI(m_mix, molecular_weights, nB)
    lambda_s = calculate_lambda(rho_mix, p_mix, molecular_weights, nB)

    # Maxwell–Stefan coupling after streaming
    Chi_S = post_stream_Chi_S(CHI_sc, rho_s, rho_mix)
    ux_dagger, uy_dagger, _, _ = solve_ms_fluxes(
        lambda_s, Chi_S, CHI_sc, rho_s, rho_mix, ux_s, uy_s, theta=theta
    )
    ux_star_dagger_s, uy_star_dagger_s = calculate_u_star(
        CHI_sc, rho_s, rho_mix, ux_dagger, uy_dagger
    )
    feq_dagger = equilibrium(f_streamed, rho_s, phi, ux_star_dagger_s, uy_star_dagger_s)

    # Second half-step
    f_new = distribution_semi_implicit(feq_dagger, g_dagger_s, lambda_s)
    return f_new


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_frame(concA: np.ndarray, frame_idx: int, out_dir: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    im = ax.imshow(concA.T, origin="lower", cmap="viridis", aspect="auto")
    ax.set_title(f"Species A concentration (frame {frame_idx:04d})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig( f"demo_frames/frame_{frame_idx:04d}.png", dpi=160)
    plt.close(fig)


def main():
    config = Config()

    # Output directory for frames
    #out_dir = Path(__file__).resolve().parent / config.frames_dir
    #ensure_dir(out_dir)

    out_dir = 'demo_frames'

    # Setup
    f, phi = initialise_stripe(config)
    molecular_weights = np.array(config.molecular_weights, dtype=np.float64)

    # Run
    for step in tqdm(range(config.steps + 1)):
        # Diagnostics and frame writing
        if step % config.output_stride == 1:
            rho_s, _, _, rho_mix, _ = calculate_moment(f, phi)
            rho_mix_safe = np.where(rho_mix > 0.0, rho_mix, 1.0)
            concA = rho_s[0] / rho_mix_safe
            save_frame(concA, step, out_dir)

        if step == config.steps:
            break

        f = bgk_step(f, phi, molecular_weights, nB=config.nB, theta=config.theta)

    print(f"Saved frames to {out_dir}")


if __name__ == "__main__":
    main()
