import numpy as np
import matplotlib.pyplot as plt
use_cupy = True

if use_cupy:

    try:
        import cupy as cp
        xp = cp
        print("Using cupy")
    except ImportError:
        cp = None
        xp = np
else:
    cp = None
    xp = np

def to_numpy(arr):
    if cp is not None and hasattr(cp, "asnumpy") and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return np.asarray(arr)


D2Q9_CX = xp.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=xp.int64) # CW: up, ur, r, dr, d.....
D2Q9_CY = xp.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=xp.int64)

w = xp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
OPPOSITE = xp.array([0,5,6,7,8,1,2,3,4], dtype=xp.int64)
theta = 0.5


def _safe_divide(numerator, denominator, mask=None):
    """
    Elementwise division that avoids unsupported CuPy `where`/`out` kwargs.
    Returns zero where the mask is False.
    """
    if mask is None:
        mask = denominator != 0
    denom_safe = xp.where(mask, denominator, xp.ones_like(denominator))
    result = numerator / denom_safe
    return xp.where(mask, result, xp.zeros_like(result))


def lattice_stream(f, phi, step, non_absorb_mask, bc_top, bc_bottom):

    f_streamed = xp.zeros_like(f)
    for i in range(w.shape[0]): # for each species
        f_streamed[:, i, :, :] = xp.roll(f[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    return f_streamed

def lattice_stream_BC_full(g_dagger_s, phi, step, non_absorb_mask, bc_top, bc_bottom):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]):
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    #bottom absorption:
    b1, b2, b3, reflection_boundary = bc_bottom

    g_streamed = top_bottom_boundary('bottom', g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask)

    # top:
    b1, b2, b3, reflection_boundary = bc_top

    g_streamed = top_bottom_boundary('top', g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask)

    '''
    if step % 200 == 0:
        molecular_weight = xp.array([28.0, 2.0, 44.0])
        nB = 2
        rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(g_streamed, phi)
        m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weight)

        CHI_sc = calculate_CHI(m_mix, molecular_weight, nB)

        ux_star_s, uy_star_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s)

        #plot_vector(cp.asnumpy(ux_star_s[0, :, :]), cp.asnumpy(uy_star_s[0, :, :]), zoom = 2)
        #plt.plot(cp.asnumpy(ux_star_s[0, :, :])[:, 200])
    '''
    return g_streamed

def top_bottom_boundary(direction, g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask):

    if direction == 'top':
        boundary_indexes = xp.array([4,5,6])
        row_index = -1
    elif direction == 'bottom':
        boundary_indexes = xp.array([1,2,8])
        row_index = 0

    dx = 1
    ne = 1

    for s in range(g_dagger_s.shape[0]):

        if reflection_boundary[s] == 0:  # if boundary is not reflection:

            if b1[s] == 0:  # Dirichlet boundary
                C_w = xp.ones(g_dagger_s.shape[2]) * b3[s] / b2[s]

            else: # absorption
                C_f = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)  # last component, bottom row

                C_w = (C_f + ne * dx * b3[s] / b1[s] / 2) / (1 + ne * dx * b2[s] / b1[s] / 2)

            # equilibrium based on concentration
            f_w_bottom = eq_single_boundary(C_w, phi[s], xp.zeros(g_dagger_s.shape[2]),
                                            xp.zeros(g_dagger_s.shape[2]))

            for i in boundary_indexes:
                g_streamed[s, i, :, row_index] = -g_dagger_s[s, OPPOSITE[i], :, row_index] + 2 * f_w_bottom[i]

            if b1[s] != 0:  # for absorption condition, overwrite masked area to no-slip wall
                #nx, a = g_dagger_s.shape[2], g_dagger_s.shape[2] // 5
                #non_absorbing = (xp.arange(nx) // a) % 2

                for i in boundary_indexes:
                    g_streamed[s, i, (non_absorb_mask == 1), row_index] = g_dagger_s[s, OPPOSITE[i], (non_absorb_mask == 1), row_index]

        else:  # no-slip wall:
            for i in boundary_indexes:
                g_streamed[s, i, :, row_index] = g_dagger_s[s, OPPOSITE[i], :, row_index]

    return g_streamed


def eq_single_boundary(C_w, phi_s, ux_star, uy_star):
    u_sq = ux_star ** 2 + uy_star ** 2  # (nx or...y...)
    cu = D2Q9_CX[:, None] * ux_star[None, :] + D2Q9_CY[:, None] * uy_star[None, :]

    f_w_eq = w[:, None] * C_w[None, :] * (
            phi_s
            + 3 * cu
            + 4.5 * cu ** 2
            - 1.5 * u_sq
    )

    f_w_eq[0, :] = w[0] * C_w * ((9 - 5 * phi_s) / 4 - 1.5 * u_sq)

    return f_w_eq

def calculate_moment(f, phi):

    rho_s = xp.sum(f, axis = 1)
    ux_s = xp.sum(f * D2Q9_CX.reshape(1,9,1,1), axis = 1) / rho_s
    uy_s = xp.sum(f * D2Q9_CY.reshape(1,9,1,1), axis = 1) / rho_s

    rho_mix = xp.clip(xp.sum(rho_s, axis = 0), a_min = 0, a_max= xp.inf)
    p_mix = xp.clip(xp.sum(rho_s * phi[:, None, None]/3, axis = 0), a_min = 0, a_max = xp.inf)


    return rho_s, ux_s, uy_s, rho_mix, p_mix

def calculate_m_mix(rho_s, rho_mix, molecular_weight):

    denom = molecular_weight[:, None, None] * rho_mix[None, :, :]
    mask = denom > 0
    divider_inv_M = _safe_divide(rho_s, denom, mask=mask)

    inv_M_sum = xp.sum(divider_inv_M, axis = 0)

    mask_mix = inv_M_sum > 0
    m_mix = _safe_divide(xp.ones_like(inv_M_sum), inv_M_sum, mask=mask_mix)

    return m_mix

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
        dtype=xp.float64,
    )
    factor = table[nB - 1]
    return factor * 10.0 * (1.0 / m1 + 1.0 / m2) ** (-0.5)

def calculate_CHI(m_mix, molecular_weight, nB):

    N_species = len(molecular_weight)
    nx, ny = m_mix.shape
    CHI_sc = xp.zeros((N_species, N_species, nx, ny))

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)

        for c in range(N_species):
            B_sc = B_binary_resistivity(molecular_weight[s], molecular_weight[c], nB)

            CHI_sc[s, c, :, :] = m_mix **2 * B_sc / (molecular_weight[s] * molecular_weight[c] * B_ss)

    return CHI_sc

def calculate_lambda(rho_mix, p_mix, molecular_weight, nB):
    N_species = len(molecular_weight)
    nx, ny = rho_mix.shape
    lambda_s = xp.zeros((N_species, nx, ny))

    for s in range(N_species):
        B_ss = B_binary_resistivity(molecular_weight[s], molecular_weight[s], nB)
        lambda_s[s, :, :] = p_mix * B_ss / rho_mix

    return lambda_s

def calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s):
    """
    Compute the modified species velocities u* used in the BGK half-step.

    The expression matches the reference implementation:
        u*_s = u_s + sum_v CHI_{s,v} (rho_v / rho_mix) (u_v - u_s)
    where CHI_sc carries the per-node resistivity ratios and rho_mix is the
    mixture density.
    """
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = _safe_divide(
        rho_s,
        rho_mix_safe[None, :, :],
        mask=(rho_mix[None, :, :] > 0.0),
    )

    delta_ux = ux_s[None, :, :, :] - ux_s[:, None, :, :]
    delta_uy = uy_s[None, :, :, :] - uy_s[:, None, :, :]

    chi_weight = CHI_sc * rho_ratio[None, :, :, :]
    ux_star = ux_s + xp.sum(chi_weight * delta_ux, axis=1)
    uy_star = uy_s + xp.sum(chi_weight * delta_uy, axis=1)
    return ux_star, uy_star

def equilibrium(f, rho_s, phi, ux_star_s, uy_star_s):

    #N_species = len(phi)
    feq = xp.zeros_like(f)

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

    g_s = f - theta * lambda_s[:, None, :, :] * (feq - f)
    g_dagger = g_s + (lambda_s/(1 + theta * lambda_s))[:, None, :, :] * (feq - g_s)

    return g_dagger


def post_stream_Chi_S(CHI_sc, rho_s, rho_mix):
    Chi_S = xp.zeros_like(rho_s)

    for s in range(rho_s.shape[0]):
        Chi_S[s] = xp.sum(CHI_sc[s, :, :, :] * rho_s / rho_mix[None, :, :], axis=0)

    return Chi_S

def distribution_semi_implicit(feq_dagger, g_dagger_s, lambda_s):
    f_new = (
        (g_dagger_s + theta * lambda_s[:, None, :, :] * feq_dagger)
        / ((1 + theta * lambda_s)[:, None, :, :])
    )

    f_new = xp.clip(f_new, a_min = 0, a_max = xp.inf)

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
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = _safe_divide(
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

    ux_updated = _safe_divide(jx_species, rho_s, mask=(rho_s > 0.0))
    uy_updated = _safe_divide(jy_species, rho_s, mask=(rho_s > 0.0))

    ux_updated = xp.nan_to_num(ux_updated)
    uy_updated = xp.nan_to_num(uy_updated)
    return ux_updated, uy_updated, jx_species, jy_species

def bgk_step(f, molecular_weight, phi, nB, stream_fn, step, non_absorb_mask, bc_top, bc_bottom):

    #################### Before streaming ####################

    rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(f, phi)
    m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weight)
    CHI_sc = calculate_CHI(m_mix, molecular_weight, nB)
    lambda_s = calculate_lambda(rho_mix, p_mix, molecular_weight, nB)

    ux_star_s, uy_star_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s)
    feq = equilibrium(f, rho_s, phi, ux_star_s, uy_star_s)

    g_dagger_s =calculate_g_dagger(f, feq, lambda_s)

    #################### After streaming ####################

    g_streamed = stream_fn(g_dagger_s, phi, step,
                           non_absorb_mask, bc_top, bc_bottom)

    rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(g_streamed, phi)
    m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weight)

    lambda_s = calculate_lambda(rho_mix, p_mix, molecular_weight, nB)
    CHI_sc = calculate_CHI(m_mix, molecular_weight, nB)

    Chi_S = post_stream_Chi_S(CHI_sc, rho_s, rho_mix)
    ux_dagger, uy_dagger, jx_s, jy_s = solve_ms_fluxes(lambda_s, Chi_S, CHI_sc, rho_s, rho_mix, ux_s, uy_s, theta=theta)

    ux_star_dagger_s, uy_star_dagger_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_dagger, uy_dagger)

    feq_dagger = equilibrium(g_streamed, rho_s, phi, ux_star_dagger_s, uy_star_dagger_s)
    f_new = distribution_semi_implicit(feq_dagger, g_streamed, lambda_s)

    return f_new
