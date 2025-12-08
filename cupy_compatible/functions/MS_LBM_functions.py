from .use_cupy import *


D2Q9_CX = xp.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=xp.float32) # CW: up, ur, r, dr, d.....
D2Q9_CY = xp.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=xp.float32)

w = xp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36], dtype=xp.float32)
OPPOSITE = xp.array([0,5,6,7,8,1,2,3,4], dtype=xp.float32)
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

def calculate_moment(f, phi):

    rho_s = xp.sum(f, axis = 1, dtype = xp.float32)
    ux_s = xp.sum(f * D2Q9_CX.reshape(1,9,1,1), axis = 1, dtype=xp.float32) / rho_s
    uy_s = xp.sum(f * D2Q9_CY.reshape(1,9,1,1), axis = 1, dtype=xp.float32) / rho_s

    rho_mix = xp.clip(xp.sum(rho_s, axis = 0, dtype=xp.float32), a_min = 0, a_max= xp.inf)
    p_mix = xp.clip(xp.sum(rho_s * phi[:, None, None]/3, axis = 0, dtype=xp.float32), a_min = 0, a_max = xp.inf)

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
        dtype=xp.float32,
    )
    factor = table[nB]
    return factor * 10.0 * (1.0 / m1 + 1.0 / m2) ** (-0.5)

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

def calculate_g_dagger(f, feq, lambda_s):

    g_s = f - theta * lambda_s[:, None, :, :] * (feq - f)
    g_dagger = g_s + (lambda_s/(1 + theta * lambda_s))[:, None, :, :] * (feq - g_s)

    return g_dagger


def post_stream_Chi_S(CHI_sc, rho_s, rho_mix):
    Chi_S = xp.zeros_like(rho_s, dtype=xp.float32)

    for s in range(rho_s.shape[0]):
        Chi_S[s] = xp.sum(CHI_sc[s, :, :, :] * rho_s / rho_mix[None, :, :], axis=0)

    return Chi_S

def distribution_semi_implicit(feq_dagger, g_dagger_s, lambda_s):
    f_new = (
        (g_dagger_s + theta * lambda_s[:, None, :, :] * feq_dagger)
        / ((1 + theta * lambda_s)[:, None, :, :])
    )

    if xp.any(f_new < 0):
        print("clipped")

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

    #print(f_new.dtype)

    #print("lambda ", xp.max(lambda_s), xp.min(lambda_s))

    #print("vx ", xp.mean(ux_s, axis = (-2, -1)))
    #print(xp.mean(uy_s, axis = (-2, -1)))

    return f_new

'''
def post_stream_Chi_S(CHI_sc, rho_s, rho_mix):
    """
    Compute Chi_S[s, x, y] = sum_c CHI_sc[s, c, x, y] * rho_c(x,y)/rho_mix(x,y)
    (MSGas uses these sums in the diagonal entries of A.)
    """
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)

    # rho_ratio[c,x,y] = rho_c / rho_mix
    rho_ratio = _safe_divide(
        rho_s,
        rho_mix_safe[None, :, :],
        mask=(rho_mix[None, :, :] > 0.0),
    )  # shape (N, nx, ny)

    # Broadcast over the second species index of CHI_sc (c index)
    rho_ratio_exp = rho_ratio[None, :, :, :]          # shape (1, N, nx, ny)

    # Sum over c: Chi_S[s] = sum_c CHI_sc[s,c] * rho_c/rho_mix
    Chi_S = xp.sum(CHI_sc * rho_ratio_exp, axis=1)    # -> (N, nx, ny)

    # Zero out nodes with rho_mix == 0
    Chi_S = xp.where(rho_mix[None, :, :] > 0.0, Chi_S, 0.0)

    return Chi_S
'''
def bgk_step_alt(f, molecular_weight, phi, nB, stream_fn, step,
             non_absorb_mask, bc_top, bc_bottom):
    """
    Single BGK+MS step in the same structure as MSGas:

    1) Compute moments from pre-collision PDFs f.
    2) Build Maxwell–Stefan coefficients (m_mix, CHI_sc, lambda_s, Chi_S).
    3) Solve A * j = g for species momenta j (MS-corrected fluxes).
    4) Convert to velocities u_sigma and build drift velocities u*_sigma.
    5) Build equilibrium feq with u*_sigma.
    6) Do one θ-scheme BGK update using lam_fac = lambda/(1 + theta*lambda).
    7) Stream the post-collision PDFs with user-provided stream_fn.

    This mirrors mus_compute_MSGas_module.fpp::bgk_advRel_d3q19f3_MSGas
    up to:
      - 2D instead of 3D
      - generic N species instead of 3
      - collide-then-stream instead of fused pull-stream+collide.
    """

    # ------------------ 1) Macroscopic fields from f ------------------ #
    rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(f, phi)

    # ------------------ 2) MS coefficients (same as MSGas) ------------ #
    m_mix   = calculate_m_mix(rho_s, rho_mix, molecular_weight)
    CHI_sc  = calculate_CHI(m_mix, molecular_weight, nB)
    lambda_s = calculate_lambda(rho_mix, p_mix, molecular_weight, nB)
    # Avoid NaNs/Infs from empty nodes
    lambda_s = xp.nan_to_num(lambda_s, nan=0.0, posinf=0.0, neginf=0.0)

    Chi_S = post_stream_Chi_S(CHI_sc, rho_s, rho_mix)

    # ------------------ 3) Solve MS flux LSE A * j = g ---------------- #
    # This matches the MSGas matrix:
    #   A_ss = 1 + θ λ_s Σ_{c≠s} χ_{s,c} (ρ_c/ρ)
    #   A_sc = -θ λ_s (ρ_s/ρ) χ_{s,c}     (c ≠ s)
    # with RHS gσ = ρσ uσ.
    ux_ms, uy_ms, jx_s, jy_s = solve_ms_fluxes(
        lambda_s, Chi_S, CHI_sc, rho_s, rho_mix, ux_s, uy_s, theta=theta
    )

    # ------------------ 4) Drift velocities u*_σ (uxstar, uystar) ----- #
    # MSGas: uxstar(σ) = ux_ms(σ) + Σ_c χ_{σc} (ρ_c/ρ) [ux_ms(c)–ux_ms(σ)]
    ux_star_s, uy_star_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_ms, uy_ms)

    # ------------------ 5) Equilibrium with MS-corrected u* ----------- #
    feq = equilibrium(f, rho_s, phi, ux_star_s, uy_star_s)

    # ------------------ 6) θ-scheme BGK update (lam_fac) --------------- #
    # MSGas uses lam_fac = λ/(1+θλ), lamfac_o = 1 - lam_fac
    lam_fac  = lambda_s / (1.0 + theta * lambda_s)
    lamfac_o = 1.0 - lam_fac

    # Broadcast over D2Q9 directions
    lam_fac_b  = lam_fac[:, None, :, :]   # (Ns, 1, nx, ny)
    lamfac_o_b = lamfac_o[:, None, :, :]

    # f_post = (1 - lam_fac) f + lam_fac feq
    f_post_collision = lamfac_o_b * f + lam_fac_b * feq
    f_post_collision = xp.clip(f_post_collision, a_min=0.0, a_max=xp.inf)

    # ------------------ 7) Streaming + BCs ---------------------------- #
    # MSGas does a pull-stream inside the same kernel; here we do
    # collide-then-stream, which is equivalent for this case.
    f_streamed = stream_fn(
        f_post_collision,
        phi,
        step,
        non_absorb_mask,
        bc_top,
        bc_bottom,
    )

    return f_streamed

