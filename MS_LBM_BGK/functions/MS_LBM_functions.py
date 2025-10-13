import numpy as np

D2Q9_CX = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int64)
D2Q9_CY = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)

w = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
#w = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPPOSITE = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64)
theta = 0.5


def lattice_stream(f):

    f_streamed = np.zeros_like(f)
    for i in range(w.shape[0]):
        f_streamed[:, i, :, :] = np.roll(f[:, i, :, :], (D2Q9_CX[i], D2Q9_CY[i]), axis = (1,2))

    return f_streamed

def calculate_moment(f, phi):

    rho_s = np.sum(f, axis = 1)
    ux_s = np.sum(f * D2Q9_CX.reshape(1,9,1,1), axis = 1) / rho_s
    uy_s = np.sum(f * D2Q9_CY.reshape(1,9,1,1), axis = 1) / rho_s

    rho_mix = np.clip(np.sum(rho_s, axis = 0), a_min = 0, a_max= np.inf)
    p_mix = np.clip(np.sum(rho_s * phi[:, None, None]/3, axis = 0), a_min = 0, a_max = np.inf)


    return rho_s, ux_s, uy_s, rho_mix, p_mix

def calculate_m_mix(rho_s, rho_mix, molecular_weight):

    inv_M_sum = np.sum(rho_s / ( molecular_weight[:, None, None] * rho_mix[None, :, :] ), axis = 0)

    return np.clip(1/inv_M_sum, a_min = 0, a_max= np.inf)

def B_binary_resistivity(m1: float, m2: float, nB: int) -> float:
    """
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
    table = np.array(
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
        dtype=np.float64,
    )
    factor = table[nB - 1]
    return factor * 10.0 * (1.0 / m1 + 1.0 / m2) ** (-0.5)

def calculate_CHI(m_mix, molecular_weight, nB):

    N_species = len(molecular_weight)
    nx, ny = m_mix.shape
    CHI_sc = np.zeros((N_species, N_species, nx, ny))

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)

        for c in range(N_species):
            B_sc = B_binary_resistivity(molecular_weight[s], molecular_weight[c], nB)

            CHI_sc[s, c, :, :] = m_mix **2 * B_sc / (molecular_weight[s] * molecular_weight[c] * B_ss)

    return CHI_sc

def calculate_lambda(rho_mix, p_mix, molecular_weight, nB):
    N_species = len(molecular_weight)
    nx, ny = rho_mix.shape
    lambda_s = np.zeros((N_species, nx, ny))

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)
        lambda_s[s, :, :] = p_mix * B_ss / rho_mix

    return lambda_s

'''
def calculate_u_star_original(CHI_sc, rho_s, rho_mix, ux_s, uy_s):
    N_species = rho_s.shape[0]
    nx, ny = rho_mix.shape
    ux_star_s = np.zeros((N_species, nx, ny))
    uy_star_s = np.zeros((N_species, nx, ny))
    for s in range(N_species):
        delta_ux = ux_s[s][None,:] - ux_s[:]
        delta_uy = uy_s[s][None,:] - uy_s[:]
        ux_star_s[s] = ux_s[s] + np.sum(CHI_sc[s, :, :, :] * rho_s * delta_ux, axis = 0)
        uy_star_s[s] = uy_s[s] + np.sum(CHI_sc[s, :, :, :] * rho_s * delta_uy, axis = 0)
    return ux_star_s, uy_star_s
'''
def calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s):
    """
    Compute the modified species velocities u* used in the BGK half-step.

    The expression matches the reference implementation:
        u*_s = u_s + sum_v CHI_{s,v} (rho_v / rho_mix) (u_v - u_s)
    where CHI_sc carries the per-node resistivity ratios and rho_mix is the
    mixture density.
    """
    rho_mix_safe = np.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = np.divide(
        rho_s,
        rho_mix_safe[None, :, :],
        out=np.zeros_like(rho_s),
        where=rho_mix_safe[None, :, :] > 0.0,
    )

    delta_ux = ux_s[None, :, :, :] - ux_s[:, None, :, :]
    delta_uy = uy_s[None, :, :, :] - uy_s[:, None, :, :]

    chi_weight = CHI_sc * rho_ratio[None, :, :, :]
    ux_star = ux_s + np.sum(chi_weight * delta_ux, axis=1)
    uy_star = uy_s + np.sum(chi_weight * delta_uy, axis=1)
    return ux_star, uy_star

def equilibrium(f, rho_s, phi, ux_star_s, uy_star_s):

    #N_species = len(phi)
    feq = np.zeros_like(f)

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

def calculate_g_dagger(f, feq, lambda_s):

    g_s = f - theta * lambda_s[:, None, None, None] * (feq - f)
    g_dagger = g_s + (lambda_s/(1 + theta * lambda_s))[:, None, None, None] * (feq - g_s)

    return g_dagger


def post_stream_Chi_S(CHI_sc, rho_s, rho_mix):
    Chi_S = np.zeros_like(rho_s)

    for s in range(rho_s.shape[0]):
        Chi_S[s] = np.sum(CHI_sc[s, :, :, :] * rho_s / rho_mix[None, :, :], axis=0)

    return Chi_S

def distribution_semi_implicit(feq_dagger, g_dagger_s, lambda_s):
    f_new = (
        (g_dagger_s + theta * lambda_s[:, None, :, :] * feq_dagger)
        / ((1 + theta * lambda_s)[:, None, :, :])
    )

    return f_new

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
    rho_mix_safe = np.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = np.divide(
        rho_s,
        rho_mix_safe[None, :, :],
        out=np.zeros_like(rho_s),
        where=rho_mix_safe[None, :, :] > 0.0,
    )

    theta_lambda = np.nan_to_num(theta * lambda_s, nan=0.0, posinf=0.0, neginf=0.0)
    theta_lambda_exp = np.moveaxis(theta_lambda, 0, -1)  # (nx, ny, N)
    rho_ratio_exp = np.moveaxis(rho_ratio, 0, -1)
    chi_s = np.nan_to_num(np.moveaxis(Chi_S, 0, -1), nan=0.0, posinf=0.0, neginf=0.0)
    CHI = np.nan_to_num(CHI_sc.transpose(2, 3, 0, 1), nan=0.0, posinf=0.0, neginf=0.0)

    valid = rho_mix > 0.0
    theta_lambda_exp = np.where(valid[:, :, None], theta_lambda_exp, 0.0)
    rho_ratio_exp = np.where(valid[:, :, None], rho_ratio_exp, 0.0)
    chi_s = np.where(valid[:, :, None], chi_s, 0.0)

    A = -theta_lambda_exp[:, :, :, None] * rho_ratio_exp[:, :, :, None] * CHI
    for idx in range(num_species):
        A[:, :, idx, idx] += 1.0 + theta_lambda_exp[:, :, idx] * chi_s[:, :, idx]

    gjx = np.moveaxis(rho_s * ux_s, 0, -1)
    gjy = np.moveaxis(rho_s * uy_s, 0, -1)
    gjx = np.where(valid[:, :, None], gjx, 0.0)
    gjy = np.where(valid[:, :, None], gjy, 0.0)

    jx = np.linalg.solve(A, gjx)
    jy = np.linalg.solve(A, gjy)

    jx_species = np.moveaxis(jx, -1, 0)
    jy_species = np.moveaxis(jy, -1, 0)

    ux_updated = np.divide(
        jx_species,
        rho_s,
        out=np.zeros_like(ux_s),
        where=rho_s > 0.0,
    )
    uy_updated = np.divide(
        jy_species,
        rho_s,
        out=np.zeros_like(uy_s),
        where=rho_s > 0.0,
    )

    ux_updated = np.nan_to_num(ux_updated)
    uy_updated = np.nan_to_num(uy_updated)
    return ux_updated, uy_updated, jx_species, jy_species

