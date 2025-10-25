from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

#from MS_LBM_BGK.MS_LBM import molecular_weight
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
    lattice_stream,
)
from tqdm import tqdm

def lattice_stream_bounce_xy(f: np.ndarray) -> np.ndarray:
    """Stream with bounce-back on all four boundaries (no-slip walls).

    For each direction k, if the target (i+cx[k], j+cy[k]) is outside the
    domain, reflect into the opposite direction at (i, j). Otherwise, stream
    into the interior target.
    """
    num_species, q, nx, ny = f.shape
    out = np.zeros_like(f)
    for k in range(q):
        dx = int(D2Q9_CX[k])
        dy = int(D2Q9_CY[k])
        if dx == 0 and dy == 0:
            out[:, k, :, :] += f[:, k, :, :]
            continue
        for i in range(nx):
            ii = i + dx
            for j in range(ny):
                jj = j + dy
                if (0 <= ii < nx) and (0 <= jj < ny):
                    out[:, k, ii, jj] += f[:, k, i, j]
                else:
                    out[:, OPPOSITE[k], i, j] += f[:, k, i, j]
    return out

@dataclass
class Config:
    nx: int = 200
    ny: int = 100
    steps: int = 1000
    output_stride: int = 50
    molecular_weights: tuple[float, float, float] = (28.0, 2.0, 44.0)  # N2, H2, CO2
    # Left half: 50% N2 + 50% H2, Right half: 50% N2 + 50% CO2
    #left_frac: tuple[float, float, float] = (0.5, 0.5, 0.0)
    #right_frac: tuple[float, float, float] = (0.5, 0.0, 0.5)
    left_frac: tuple[float, float, float] = (0.5, 0.4, 0.1)
    right_frac: tuple[float, float, float] = (0.5, 0.1, 0.4)
    total_pressure: float = 1.0
    theta: float = THETA
    nB: int = 10
    frames_dir: str = "demo_frames/three_species_mixing"

def initialise_chamber(config: Config):
    nx, ny = config.nx, config.ny
    species = 3
    phi = 1.0 / np.array(config.molecular_weights, dtype=np.float64)

    # Partial pressures for each species
    psigma = np.zeros((species, nx, ny), dtype=np.float64)
    mid = nx // 2
    left = np.array(config.left_frac, dtype=np.float64) * config.total_pressure
    right = np.array(config.right_frac, dtype=np.float64) * config.total_pressure
    for s in range(species):
        psigma[s, :mid, :] = left[s]
        psigma[s, mid:, :] = right[s]

    # Convert partial pressures to densities: p_s = phi_s * rho_s / 3
    rho_s = np.zeros_like(psigma)
    for s in range(species):
        rho_s[s] = 3.0 * psigma[s] / phi[s]

    # Initial populations from equilibrium with zero velocity
    f = np.zeros((species, 9, nx, ny), dtype=np.float64)
    ux_s = np.zeros((species, nx, ny), dtype=np.float64)
    uy_s = np.zeros((species, nx, ny), dtype=np.float64)
    feq = equilibrium(f, rho_s, phi, ux_s, uy_s)
    f[...] = feq
    return f, phi

def bgk_step(
    f: np.ndarray,
    phi: np.ndarray,
    molecular_weights: np.ndarray,
    nB: int,
    theta: float,
    stream_fn=lattice_stream_bounce_xy,
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

def save_concentration_frames(f: np.ndarray, phi: np.ndarray, frame_idx: int, out_dir: Path, molecular_weights) -> None:
    rho_s, _, _, rho_mix, _ = calculate_moment(f, phi)
    rho_mix_safe = np.where(rho_mix > 0.0, rho_mix, 1.0)
    conc = rho_s #/ rho_mix_safe[None, :, :]

    labels = ["N2", "H2", "CO2"]
    cmaps = ["Blues", "Greens", "Reds"]

    fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True, sharey=True)
    for s in range(3):
        im = axes[s].imshow((conc[s].T)/molecular_weights[s], origin="lower", cmap=cmaps[s], aspect="auto")
        axes[s].set_title(f"{labels[s]} concentration")
        axes[s].set_xlabel("x")
        axes[s].set_ylabel("y")
        fig.colorbar(im, ax=axes[s], fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(f"demo_frames_triple/frame_{frame_idx:04d}.png", dpi=160)
    plt.close(fig)

def main():
    cfg = Config()
    out_dir = 'demo_frames'

    f, phi = initialise_chamber(cfg)
    molecular_weights = np.array(cfg.molecular_weights, dtype=np.float64)

    for step in tqdm(range(cfg.steps + 1)):
        if step % cfg.output_stride == 0:
            save_concentration_frames(f, phi, step, out_dir, molecular_weights)
        if step == cfg.steps:
            break
        f = bgk_step(f, phi, molecular_weights, nB=cfg.nB, theta=cfg.theta)

    print(f"Saved frames to {out_dir}")