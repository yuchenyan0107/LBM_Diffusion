from .common import *

def equilibrium(f, rho_s, phi, ux_star_s, uy_star_s):

    #N_species = len(phi)
    feq = xp.zeros_like(f, dtype=xp.float32)

    u_sq = ux_star_s ** 2 + uy_star_s ** 2 # (N_species, nx, ny)
    cu = D2Q9_CX[None, :, None, None] * ux_star_s[:, None, :, :] + D2Q9_CY[None, :, None, None] * uy_star_s[:, None, :, :]

    feq = w[None, :, None, None] * rho_s[:, None, :, :] * (
        phi[:, None, None, None]
        + 3 * cu
        + 4.5 * cu**2
        -1.5 * u_sq[:, None, :, :]
    )

    feq[:, 0, :, :] = w[0] * rho_s * ((9 - 5 * phi[:, None, None]) / 4 - 1.5 * u_sq)

    return feq

def post_stream_Chi_S(CHI_sc, rho_s, rho_mix):
    Chi_S = xp.zeros_like(rho_s, dtype=xp.float32)

    for s in range(rho_s.shape[0]):
        Chi_S[s] = xp.sum(CHI_sc[s, :, :, :] * rho_s / rho_mix[None, :, :], axis=0)

    return Chi_S

def calculate_CHI(m_mix, molecular_weight, nB):

    N_species = len(molecular_weight)
    nx, ny = m_mix.shape
    CHI_sc = xp.zeros((N_species, N_species, nx, ny), dtype=xp.float32)

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)

        for c in range(N_species):
            B_sc = B_binary_resistivity(molecular_weight[s], molecular_weight[c], nB)

            CHI_sc[s, c, :, :] = m_mix **2 * B_sc / (molecular_weight[s] * molecular_weight[c] * B_ss)

    return CHI_sc

def calculate_lambda(rho_mix, p_mix, molecular_weight, nB):
    N_species = len(molecular_weight)
    nx, ny = rho_mix.shape
    lambda_s = xp.zeros((N_species, nx, ny), dtype=xp.float32)

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)
        lambda_s[s, :, :] = p_mix * B_ss / rho_mix

    return lambda_s

def solve_ms_fluxes(lambda_s, Chi_S, CHI_sc, rho_s, rho_mix, ux_s, uy_s, theta=theta):
    """
    Solve the Maxwell-Stefan flux system A * j = g across the full lattice.

    Parameters
    ----------
    lambda_s : np.ndarray
        Relaxation factors per species (N_species, nx, ny).
    Chi_S : np.ndarray
        Species-projected resistivity sums matching CHI_sigma (N_species, nx, ny).
    CHI_sc : np.ndarray
        Pairwise resistivity ratios (N_species, N_species, nx, ny).
    rho_s : np.ndarray
        Species densities (N_species, nx, ny).
    rho_mix : np.ndarray
        Mixture density (nx, ny).
    ux_s, uy_s : np.ndarray
        Species velocities before the MS coupling step (N_species, nx, ny).
    theta : float, optional
        Splitting parameter (defaults to global theta).

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Updated species velocities (ux, uy) and the flux vectors (jx, jy).
    """
    num_species, nx, ny = rho_s.shape
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = safe_divide(
        rho_s,
        rho_mix_safe[None, :, :],
        mask=(rho_mix[None, :, :] > 0.0),
    )

    theta_lambda = xp.nan_to_num(theta * lambda_s, nan=0.0, posinf=0.0, neginf=0.0)
    theta_lambda_exp = xp.moveaxis(theta_lambda, 0, -1)  # (nx, ny, N)
    rho_ratio_exp = xp.moveaxis(rho_ratio, 0, -1)
    chi_s = xp.nan_to_num(xp.moveaxis(Chi_S, 0, -1), nan=0.0, posinf=0.0, neginf=0.0)
    CHI = xp.nan_to_num(CHI_sc.transpose(2, 3, 0, 1), nan=0.0, posinf=0.0, neginf=0.0)

    valid = rho_mix > 0.0
    theta_lambda_exp = xp.where(valid[:, :, None], theta_lambda_exp, 0.0)
    rho_ratio_exp = xp.where(valid[:, :, None], rho_ratio_exp, 0.0)
    chi_s = xp.where(valid[:, :, None], chi_s, 0.0)

    A = -theta_lambda_exp[:, :, :, None] * rho_ratio_exp[:, :, :, None] * CHI
    for idx in range(num_species):
        A[:, :, idx, idx] += 1.0 + theta_lambda_exp[:, :, idx] * chi_s[:, :, idx]

    gjx = xp.moveaxis(rho_s * ux_s, 0, -1)
    gjy = xp.moveaxis(rho_s * uy_s, 0, -1)
    gjx = xp.where(valid[:, :, None], gjx, 0.0)
    gjy = xp.where(valid[:, :, None], gjy, 0.0)

    jx = xp.linalg.solve(A, gjx)
    jy = xp.linalg.solve(A, gjy)

    jx_species = xp.moveaxis(jx, -1, 0)
    jy_species = xp.moveaxis(jy, -1, 0)

    ux_updated = safe_divide(jx_species, rho_s, mask=(rho_s > 0.0))
    uy_updated = safe_divide(jy_species, rho_s, mask=(rho_s > 0.0))

    ux_updated = xp.nan_to_num(ux_updated)
    uy_updated = xp.nan_to_num(uy_updated)
    return ux_updated, uy_updated, jx_species, jy_species

def B_binary_resistivity(m1: float, m2: float, nB: int) -> float:
    """
    This list is copied from the musubi example
    Maxwell–Stefan binary resistivity lookup B(σ, ν).

    Parameters
    ----------
    m1, m2:
        Molecular weights of the interacting species.
    nB:
        Index (1-based) into the empirical resistivity table.

    Returns
    -------
    float:
        Resistivity value consistent with MIXLBM.f90::B.
    """
    table = xp.array(
        [
            0.500000000000000,
            0.601274801898418,
            0.723062774795962,
            0.869518853351124,
            1.045639552591273,
            1.257433429682934,
            1.512126072666108,
            1.818406609575491,
            2.186724147886553,
            2.629644257653945,
            3.162277660168372,
            3.802795747331057,
            4.573050519273251,
            5.499320090094956,
            6.613205195495662,
            7.952707287670479,
            9.563524997900334,
            11.500613197126169,
            13.830057843624722,
            16.631330580338215,
        ],
        dtype=xp.float32,
    )
    factor = table[nB]
    return factor * 10.0 * (1.0 / m1 + 1.0 / m2) ** (-0.5)