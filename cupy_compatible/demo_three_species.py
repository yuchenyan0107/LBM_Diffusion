from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from functions import (
    xp,
    cp,
    equilibrium,
    calculate_moment,
    bgk_step,
    lattice_stream_bounce_xy,
)

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
    left_frac: tuple[float, float, float] = (0.5, 0.25, 0.35)
    right_frac: tuple[float, float, float] = (0.5, 0.35, 0.25)
    total_pressure: float = 1.0
    theta: float = 0.5
    nB: int = 2

def initialise_chamber(config: Config):
    nx, ny = config.nx, config.ny
    species = 3
    phi = 1.0 / xp.array(config.molecular_weights, dtype=xp.float64)

    # Partial pressures for each species
    psigma = xp.zeros((species, nx, ny), dtype=xp.float64)
    mid = nx // 2
    left = xp.array(config.left_frac, dtype=xp.float64) * config.total_pressure
    right = xp.array(config.right_frac, dtype=xp.float64) * config.total_pressure
    for s in range(species):
        psigma[s, :mid, :] = left[s]
        psigma[s, mid:, :] = right[s]

    # Convert partial pressures to densities: p_s = phi_s * rho_s / 3
    rho_s = xp.zeros_like(psigma)
    for s in range(species):
        rho_s[s] = 3.0 * psigma[s] / phi[s]

    # Initial populations from equilibrium with zero velocity
    f = xp.zeros((species, 9, nx, ny), dtype=xp.float64)
    ux_s = xp.zeros((species, nx, ny), dtype=xp.float64)
    uy_s = xp.zeros((species, nx, ny), dtype=xp.float64)
    feq = equilibrium(f, rho_s, phi, ux_s, uy_s)
    f[...] = feq
    return f, phi

def _to_numpy(arr):
    if cp is not None and hasattr(cp, "asnumpy") and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return np.asarray(arr)

def save_concentration_frames(f: np.ndarray, phi: np.ndarray, frame_idx: int, out_dir: Path, molecular_weights) -> None:
    rho_s, _, _, rho_mix, _ = calculate_moment(f, phi)
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    conc = rho_s #/ rho_mix_safe[None, :, :]

    conc_cpu = _to_numpy(conc)
    mw_cpu = _to_numpy(molecular_weights)
    labels = ["N2", "H2", "CO2"]
    cmaps = ["Blues", "Greens", "Reds"]

    fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True, sharey=True)
    for s in range(3):
        im = axes[s].imshow((conc_cpu[s].T)/mw_cpu[s], origin="lower", cmap=cmaps[s], aspect="auto")
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
    molecular_weights = xp.array(cfg.molecular_weights, dtype=xp.float64)

    for step in tqdm(range(cfg.steps + 1)):
        if step % cfg.output_stride == 0:
            save_concentration_frames(f, phi, step, out_dir, molecular_weights)
        if step == cfg.steps:
            break
        #f = bgk_step(f, phi, molecular_weights, nB=cfg.nB, theta=cfg.theta)
        f = bgk_step(f, molecular_weights, phi, cfg.nB, lattice_stream_bounce_xy)

    print(f"Saved frames to {out_dir}")
