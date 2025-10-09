"""
Maxwell–Stefan lattice Boltzmann operators using a BGK collision model.

What this implements
--------------------
- Same Maxwell–Stefan coupling and θ-splitting structure as the MRT version
  in :mod:`.maxwell_stefan`, but with a scalar BGK collision: ``A = λ I``.
- The θ-transform becomes scalar: ``G = (I + θ λ I)^{-1}`` and the two half
  steps reduce to compact per-direction formulas. The MS linear solve (χ-based
  flux coupling) stays identical.

Conceptual equation references (PDF)
-----------------------------------
- θ-splitting and transform: ``G = (I + θ A)^{-1}``.
- BGK collision operator: ``Ω = −λ (f − f_eq)`` (here embedded in θ-splitting).
- Maxwell–Stefan coupling: χ_{σν} projection and linear system on the fluxes.
"""

from __future__ import annotations

import numpy as np

from .maxwell_stefan import (
    BB_OPPOSITE,
    D2Q9_CX,
    D2Q9_CY,
    THETA,
    B_binary_resistivity,
    MaxwellStefanState,
    equilibrium_distribution,
    hydrodynamic_moments,
    solve_linear_system,
)


def _theta_split_bgk(populations: np.ndarray, feq: np.ndarray, lam: float, theta: float) -> np.ndarray:
    """Apply the θ-splitting BGK collision to ``populations`` with equilibrium ``feq``.

    The closed-form expression is obtained by setting the collision matrix ``A``
    to ``λ I`` in the θ-splitting scheme used in MIXLBM. Introducing the
    auxiliary state ``g = (I - θA) f + θA f_eq`` and solving ``(I + θA) f⁺ = g``
    yields the compact per-direction update used here (no matrix multiplications
    are needed in BGK).
    """
    if lam <= 0.0:
        return populations.copy()

    f_temp = populations + theta * lam * (populations - feq)
    correction = lam / (1.0 + theta * lam)
    return f_temp + correction * (feq - f_temp)


def maxwell_stefan_bgk_step(
    state: MaxwellStefanState,
    phi: np.ndarray,
    molecular_weights: np.ndarray,
    nN: int = 20,
    nB: int = 15,
    theta: float = THETA,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Advance the Maxwell–Stefan diffusion model by one time step using BGK collisions.

    Parameters
    ----------
    state:
        Lattice Boltzmann state (populations are read from ``state.fd``).
    phi:
        Equation-of-state factors φ_σ relating pressure and density.
    molecular_weights:
        Molecular weights M_σ for each species.
    nN, nB:
        Lookup indices for the kinematic viscosity and resistivity tables,
        matching the original MIXLBM implementation. ``nN`` is accepted for
        interface parity but is not used in the BGK update.
    theta:
        Splitting parameter; θ=0.5 delivers second-order temporal accuracy.

    Returns
    -------
    tuple(np.ndarray, np.ndarray, np.ndarray):
        Updated species densities together with the mixture velocity
        components (u_mix,x, u_mix,y).
    """
    num_species, _, nx, ny = state.fd.shape

    fd = state.fd
    fd_new = np.zeros_like(fd)

    # Scratch arrays reused during neighbour gathering and MS projections.
    rsigma_local = np.zeros(num_species, dtype=np.float64)
    psigma_local = np.zeros(num_species, dtype=np.float64)
    uxsigma_local = np.zeros(num_species, dtype=np.float64)
    uysigma_local = np.zeros(num_species, dtype=np.float64)
    uxstar = np.zeros(num_species, dtype=np.float64)
    uystar = np.zeros(num_species, dtype=np.float64)
    lambda_sigma = np.zeros(num_species, dtype=np.float64)
    CHI = np.zeros((num_species, num_species), dtype=np.float64)

    # --- First θ-half-step (f -> g) + streaming ----------------------------
    for i in range(nx):
        for j in range(ny):
            for k in range(9):
                iI = (i + D2Q9_CX[k]) % nx
                jI = (j + D2Q9_CY[k]) % ny

                # Recover neighbour moments per species.
                for s in range(num_species):
                    f_neighbor = fd[s, :, iI, jI]
                    rho, ux, uy = hydrodynamic_moments(f_neighbor)
                    rsigma_local[s] = rho
                    psigma_local[s] = phi[s] * rho / 3.0
                    uxsigma_local[s] = ux
                    uysigma_local[s] = uy

                # Mixture moments and molecular weight for MS χ build.
                r_mix = np.sum(rsigma_local)
                p_mix = np.sum(psigma_local)
                if r_mix <= 0.0:
                    r_mix = 1.0
                if p_mix <= 0.0:
                    p_mix = 1.0

                inv_MM_sum = 0.0
                # Build χ_{σν} from the tabulated resistivities.
                for s in range(num_species):
                    inv_MM_sum += (rsigma_local[s] / r_mix) / molecular_weights[s]
                mixture_mass = 1.0 / inv_MM_sum if inv_MM_sum != 0.0 else 0.0

                # MS drift-corrected velocities u*_σ.
                for s in range(num_species):
                    for vs in range(num_species):
                        resistor = B_binary_resistivity(
                            molecular_weights[s],
                            molecular_weights[vs],
                            nB,
                        )
                        ref = B_binary_resistivity(
                            molecular_weights[s],
                            molecular_weights[s],
                            nB,
                        )
                        CHI[s, vs] = (
                            (mixture_mass * mixture_mass)
                            / (molecular_weights[s] * molecular_weights[vs])
                            * resistor
                            / ref
                        )

                for s in range(num_species):
                    uxstar[s] = uxsigma_local[s]
                    uystar[s] = uysigma_local[s]
                    for vs in range(num_species):
                        factor = CHI[s, vs] * (rsigma_local[vs] / r_mix)
                        uxstar[s] += factor * (uxsigma_local[vs] - uxsigma_local[s])
                        uystar[s] += factor * (uysigma_local[vs] - uysigma_local[s])

                    # BGK relaxation frequency for species σ.
                    lambda_sigma[s] = p_mix * B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[s],
                        nB,
                    ) / r_mix

                    feq = equilibrium_distribution(
                        rsigma_local[s],
                        phi[s],
                        uxstar[s],
                        uystar[s],
                    )
                    # Apply θ-split BGK and stream only the opposite component.
                    f_coll = _theta_split_bgk(fd[s, :, iI, jI], feq, lambda_sigma[s], theta)
                    fd_new[s, BB_OPPOSITE[k], i, j] = f_coll[BB_OPPOSITE[k]]

    rsigma_post = np.zeros((num_species, nx, ny), dtype=np.float64)
    uxsigma_post = np.zeros_like(rsigma_post)
    uysigma_post = np.zeros_like(rsigma_post)
    psigma_post = np.zeros_like(rsigma_post)

    # --- Reconstruct moments after streaming (g⁺) -------------------------
    for i in range(nx):
        for j in range(ny):
            for s in range(num_species):
                rho, ux, uy = hydrodynamic_moments(fd_new[s, :, i, j])
                rsigma_post[s, i, j] = rho
                uxsigma_post[s, i, j] = ux
                uysigma_post[s, i, j] = uy
                psigma_post[s, i, j] = phi[s] * rho / 3.0

    # --- Maxwell–Stefan coupling and final θ-half-step --------------------
    for i in range(nx):
        for j in range(ny):
            r_mix = np.sum(rsigma_post[:, i, j])
            p_mix = np.sum(psigma_post[:, i, j])
            if r_mix <= 0.0:
                continue
            if p_mix <= 0.0:
                p_mix = 1.0

            inv_MM_sum = 0.0
            for s in range(num_species):
                inv_MM_sum += (rsigma_post[s, i, j] / r_mix) / molecular_weights[s]
            mixture_mass = 1.0 / inv_MM_sum if inv_MM_sum != 0.0 else 0.0

            for s in range(num_species):
                lambda_sigma[s] = p_mix * B_binary_resistivity(
                    molecular_weights[s],
                    molecular_weights[s],
                    nB,
                ) / r_mix

            for s in range(num_species):
                for vs in range(num_species):
                    resistor = B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[vs],
                        nB,
                    )
                    ref = B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[s],
                        nB,
                    )
                    CHI[s, vs] = (
                        (mixture_mass * mixture_mass)
                        / (molecular_weights[s] * molecular_weights[vs])
                        * resistor
                        / ref
                    )

            CHIsigma = np.zeros(num_species, dtype=np.float64)
            for s in range(num_species):
                temp = 0.0
                for vs in range(num_species):
                    temp += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix)
                CHIsigma[s] = temp

            A_matrix = np.zeros((num_species, num_species), dtype=np.float64)
            for s in range(num_species):
                for vs in range(num_species):
                    if s == vs:
                        A_matrix[s, vs] = 1.0 + theta * lambda_sigma[s] * CHIsigma[s]
                    else:
                        A_matrix[s, vs] = 0.0
                    A_matrix[s, vs] -= (
                        theta
                        * lambda_sigma[s]
                        * (rsigma_post[s, i, j] / r_mix)
                        * CHI[s, vs]
                    )

            # Flux moments entering the linear solve (per axis).
            gjx = rsigma_post[:, i, j] * uxsigma_post[:, i, j]
            gjy = rsigma_post[:, i, j] * uysigma_post[:, i, j]

            jx = solve_linear_system(A_matrix, gjx)
            jy = solve_linear_system(A_matrix, gjy)

            for s in range(num_species):
                if rsigma_post[s, i, j] > 0.0:
                    uxsigma_post[s, i, j] = jx[s] / rsigma_post[s, i, j]
                    uysigma_post[s, i, j] = jy[s] / rsigma_post[s, i, j]
                else:
                    uxsigma_post[s, i, j] = 0.0
                    uysigma_post[s, i, j] = 0.0

            for s in range(num_species):
                uxstar[s] = uxsigma_post[s, i, j]
                uystar[s] = uysigma_post[s, i, j]
                for vs in range(num_species):
                    uxstar[s] += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix) * (
                        uxsigma_post[vs, i, j] - uxsigma_post[s, i, j]
                    )
                    uystar[s] += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix) * (
                        uysigma_post[vs, i, j] - uysigma_post[s, i, j]
                    )

                # Second θ-half-step in BGK form at the node (no streaming here).
                feq = equilibrium_distribution(
                    rsigma_post[s, i, j],
                    phi[s],
                    uxstar[s],
                    uystar[s],
                )
                relaxation = lambda_sigma[s]
                denom = 1.0 + theta * relaxation
                if denom <= 0.0:
                    denom = 1.0
                fd_new[s, :, i, j] = (fd_new[s, :, i, j] + theta * relaxation * feq) / denom

    state.fd = fd_new
    state.fd_new = fd_new.copy()
    state.rsigma = rsigma_post
    state.psigma = psigma_post
    state.uxsigma = uxsigma_post
    state.uysigma = uysigma_post

    rho_mix = np.sum(rsigma_post, axis=0)
    ux_mix = np.sum(rsigma_post * uxsigma_post, axis=0) / np.where(rho_mix > 0.0, rho_mix, 1.0)
    uy_mix = np.sum(rsigma_post * uysigma_post, axis=0) / np.where(rho_mix > 0.0, rho_mix, 1.0)

    return rsigma_post, ux_mix, uy_mix


__all__ = [
    "maxwell_stefan_bgk_step",
    "MaxwellStefanState",
]
