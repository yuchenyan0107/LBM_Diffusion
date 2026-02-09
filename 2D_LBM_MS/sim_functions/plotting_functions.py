from .MS_LBM_functions import *

def save_concentration_frames(f: np.ndarray, frame_idx: int, lbm_config) -> None:
    rho_s, _, _, rho_mix, _ = calculate_moment(f, lbm_config)
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    conc = rho_s #/ rho_mix_safe[None, :, :]

    conc_cpu = to_numpy(conc)
    labels = ["H2", "TMIn", "PH3", "CH4"]
    cmaps = ["Blues", "Greens", "Reds", "Greys"]

    fig, axes = plt.subplots(4, 1, figsize=(6, 9), sharex=True, sharey=True)
    for s in range(f.shape[0]):
        #im = axes[s].imshow((conc_cpu[s].T)/mw_cpu[s], origin="lower", cmap=cmaps[s], aspect="auto")
        im = axes[s].imshow((conc_cpu[s].T), origin="lower", cmap=cmaps[s], aspect="auto")
        axes[s].set_title(f"{labels[s]} concentration")
        axes[s].set_xlabel("x")
        axes[s].set_ylabel("y")
        fig.colorbar(im, ax=axes[s], fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(f"frame_output/frame_{frame_idx:04d}.png", dpi=160)
    plt.close(fig)

def plot_vector(ux, uy, frame_idx, skip=10, scale=1, cmap='viridis', show_bg=False, dx=1.0, dy=1.0, zoom=None, shapes = None, save = True):
    """
    Plot a 2D vector field where ux, uy have shape (nx, ny) = (x, y).
    - x increases left→right, y increases bottom→top.
    - Arrows colored by speed |u|.
    """
    ux = np.asarray(ux)
    uy = np.asarray(uy)
    if ux.shape != uy.shape or ux.ndim != 2:
        raise ValueError("ux and uy must be 2D arrays with the same shape (nx, ny).")

    nx, ny = ux.shape
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy

    # Build grid in standard plotting orientation (rows = y, cols = x)
    X, Y = np.meshgrid(x, y, indexing='xy')  # shapes (ny, nx)

    # Transpose velocity components to match (ny, nx)
    U = ux.T  # (ny, nx)
    V = uy.T  # (ny, nx)

    # Thinning
    s = max(1, int(skip))
    Xs, Ys = X[::s, ::s], Y[::s, ::s]
    Us, Vs = U[::s, ::s] * zoom, V[::s, ::s] * zoom
    speed = np.hypot(Us, Vs)

    fig, ax = plt.subplots(figsize=(7, 5), dpi=150)

    if show_bg:
        bg = np.hypot(ux, uy).T  # (ny, nx)
        ax.imshow(
            bg, origin='lower',
            extent=[x.min(), x.max(), y.min(), y.max()],
            alpha=0.35, cmap=cmap
        )

    q = ax.quiver(
        Xs, Ys, Us, Vs, speed,
        angles='xy', scale_units='width', scale=scale,
        width=0.003, headwidth=3, headlength=4, cmap=cmap
    )

    if shapes is not None:
        ax.imshow((1 - shapes).T, cmap="gray", alpha=0.3)

    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Velocity field (color = |u|)')
    fig.colorbar(q, ax=ax, label='speed')

    plt.tight_layout()
    if save == True:
        fig.savefig(f"frame_output/vector_{frame_idx:04d}.png", dpi=160)
    else:
        plt.show()
    plt.close(fig)
    return fig, ax

