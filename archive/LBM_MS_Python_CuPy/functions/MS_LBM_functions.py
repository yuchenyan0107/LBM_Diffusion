from .common import *
from .eq_and_ms import *


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
    divider_inv_M = safe_divide(rho_s, denom, mask=mask)

    inv_M_sum = xp.sum(divider_inv_M, axis = 0)

    mask_mix = inv_M_sum > 0
    m_mix = safe_divide(xp.ones_like(inv_M_sum), inv_M_sum, mask=mask_mix)

    return m_mix

def calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s):
    """
    Compute the modified species velocities u* used in the BGK half-step.

    The expression matches the reference implementation:
        u*_s = u_s + sum_v CHI_{s,v} (rho_v / rho_mix) (u_v - u_s)
    where CHI_sc carries the per-node resistivity ratios and rho_mix is the
    mixture density.
    """
    rho_mix_safe = xp.where(rho_mix > 0.0, rho_mix, 1.0)
    rho_ratio = safe_divide(
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


def calculate_g_dagger(f, feq, lambda_s):

    g_s = f - theta * lambda_s[:, None, :, :] * (feq - f)
    g_dagger = g_s + (lambda_s/(1 + theta * lambda_s))[:, None, :, :] * (feq - g_s)

    return g_dagger

def distribution_semi_implicit(feq_dagger, g_dagger_s, lambda_s):
    f_new = (
        (g_dagger_s + theta * lambda_s[:, None, :, :] * feq_dagger)
        / ((1 + theta * lambda_s)[:, None, :, :])
    )

    if xp.any(f_new < 0):
        print("clipped")

    f_new = xp.clip(f_new, a_min = 0, a_max = xp.inf)

    return f_new

def bgk_step(f, molecular_weight, phi, nB, stream_fn, step, non_absorb_mask, bc_top, bc_bottom, inlet_vx = 0):

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
                           non_absorb_mask, bc_top, bc_bottom, inlet_vx)

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