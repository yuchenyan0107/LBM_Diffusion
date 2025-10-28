from .MS_LBM_functions import *

def save_concentration_frames(f: np.ndarray, phi: np.ndarray, frame_idx: int, molecular_weights) -> None:
    rho_s, _, _, rho_mix, _ = calculate_moment(f, phi)
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    conc = rho_s #/ rho_mix_safe[None, :, :]

    conc_cpu = to_numpy(conc)
    mw_cpu = to_numpy(molecular_weights)
    labels = ["N2", "H2", "CO2"]
    cmaps = ["Blues", "Greens", "Reds"]

    fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True, sharey=True)
    for s in range(3):
        #im = axes[s].imshow((conc_cpu[s].T)/mw_cpu[s], origin="lower", cmap=cmaps[s], aspect="auto")
        im = axes[s].imshow((conc_cpu[s].T), origin="lower", cmap=cmaps[s], aspect="auto")
        axes[s].set_title(f"{labels[s]} concentration")
        axes[s].set_xlabel("x")
        axes[s].set_ylabel("y")
        fig.colorbar(im, ax=axes[s], fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(f"demo_frames_triple/frame_{frame_idx:04d}.png", dpi=160)
    plt.close(fig)

