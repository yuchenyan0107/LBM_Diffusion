"""
Maxwell–Stefan lattice Boltzmann operators for multispecies diffusion (MRT).

Overview and references
-----------------------
- This module ports the MIXLBM Fortran reference (see
  `matlab_example/MIXLBM_rome08/MIXLBM.f90`) using the same θ-splitting
  strategy and notation. The intended theoretical reference is the
  accompanying PDF "Multi-species Lattice Boltzmann Models and Practical
  Examples" (Asinari et al.).
- The Maxwell–Stefan coupling is enforced by a linear system on the species
  mass fluxes; the LBM collision operator is handled with an MRT matrix ``A``
  and the θ-transform ``G = (I + θ A)^{-1}`` (called TRANS in the Fortran).
- Notation kept from MIXLBM:
    σ -> species index
    k -> lattice direction (D2Q9)
    rsigma, psigma -> densities and partial pressures per species
    uxsigma, uysigma -> velocity components per species
    λ (lamd), ν (lamn), ζ (lamb) -> MRT relaxation frequencies
    χ -> Maxwell–Stefan resistivities between species
    θ -> splitting parameter (0.5 for second order)

Where equations map conceptually (PDF)
-------------------------------------
- Hydrodynamic moments (ρ, ρu): standard D2Q9 zeroth/first moments of f.
- Equilibrium feq: low-Mach equilibrium with pressure factor φσ.
- MRT collision: linear operator A (see MIXLBM subroutine `MRT`).
- θ-splitting: construct G = (I + θ A)^{-1} and apply half-steps in f-space.
- Maxwell–Stefan: χ tensor and flux-coupling linear system in species space.

Equation numbers may differ between PDF versions; where possible the comments
point to the named relations (θ-transform, MRT, Maxwell–Stefan χ-coupling)
rather than a fixed numbering.
"""

import numpy as np


# Discrete velocities for the D2Q9 lattice stencil.
D2Q9_CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int64)
D2Q9_CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int64)
# Bounce-back mapping: index of the direction opposite to k.
BB_OPPOSITE = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64)
# θ parameter for the splitting scheme (0.5 -> second order accuracy).
THETA = 0.5


def hydrodynamic_moments(distribution: np.ndarray) -> tuple[float, float, float]:
    '''Compute density and velocity components for a single species.

    Parameters
    ----------
    distribution:
        Array of the nine populations f_k at one lattice node.

    Returns
    -------
    rho, ux, uy:
        Zeroth and first-order moments of f_k (density and velocities).'''
    rho = np.sum(distribution)
    if rho <= 0.0:
        return 0.0, 0.0, 0.0

    ux = np.dot(D2Q9_CX, distribution) / rho
    uy = np.dot(D2Q9_CY, distribution) / rho
    return rho, ux, uy


def equilibrium_distribution(rho: float, phi: float, ux: float, uy: float) -> np.ndarray:
    '''Construct the multispecies equilibrium distribution on D2Q9.

    Parameters
    ----------
    rho:
        Species density ρ_σ.
    phi:
        Equation-of-state factor φ_σ (relates pressure and density).
    ux, uy:
        Velocity components for the equilibrium.

    Returns
    -------
    np.ndarray:
        Equilibrium populations f_eq,k.'''
    u_sq = ux * ux + uy * uy

    feq = np.empty(9, dtype=np.float64)
    feq[0] = (4.0 / 9.0) * rho * ((9.0 - 5.0 * phi) / 4.0 - 1.5 * u_sq)

    cu = ux
    feq[1] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = uy
    feq[2] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux
    feq[3] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -uy
    feq[4] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)

    cu = ux + uy
    feq[5] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux + uy
    feq[6] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux - uy
    feq[7] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = ux - uy
    feq[8] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)

    return feq


def mrt_matrix(lamd: float, lamn: float, lamb: float) -> np.ndarray:
    '''Build the velocity-space collision matrix A for the MRT operator.

    Parameters
    ----------
    lamd:
        Relaxation frequency for density-like moments.
    lamn:
        Relaxation frequency for shear moments.
    lamb:
        Relaxation frequency for bulk moments.

    Returns
    -------
    np.ndarray:
        9x9 collision matrix matching MIXLBM.f90::MRT.'''
    lam0 = 0.0
    lam34 = 1.0

    A = np.zeros((9, 9), dtype=np.float64)

    A[0, 0] = lam0
    A[0, 1] = lam0 - lamb
    A[0, 2] = lam0 - lamb
    A[0, 3] = lam0 - lamb
    A[0, 4] = lam0 - lamb
    A[0, 5] = lam0 - 2.0 * lamb + lam34
    A[0, 6] = lam0 - 2.0 * lamb + lam34
    A[0, 7] = lam0 - 2.0 * lamb + lam34
    A[0, 8] = lam0 - 2.0 * lamb + lam34

    A[1, 1] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[1, 2] = lamb / 4.0 - lamn / 4.0
    A[1, 3] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[1, 4] = lamb / 4.0 - lamn / 4.0
    A[1, 5] = lamd / 2.0 + lamb / 2.0 - lam34
    A[1, 6] = -lamd / 2.0 + lamb / 2.0
    A[1, 7] = -lamd / 2.0 + lamb / 2.0
    A[1, 8] = lamd / 2.0 + lamb / 2.0 - lam34

    A[2, 1] = lamb / 4.0 - lamn / 4.0
    A[2, 2] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[2, 3] = lamb / 4.0 - lamn / 4.0
    A[2, 4] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[2, 5] = lamd / 2.0 + lamb / 2.0 - lam34
    A[2, 6] = lamd / 2.0 + lamb / 2.0 - lam34
    A[2, 7] = -lamd / 2.0 + lamb / 2.0
    A[2, 8] = -lamd / 2.0 + lamb / 2.0

    A[3, 1] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[3, 2] = lamb / 4.0 - lamn / 4.0
    A[3, 3] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[3, 4] = lamb / 4.0 - lamn / 4.0
    A[3, 5] = -lamd / 2.0 + lamb / 2.0
    A[3, 6] = lamd / 2.0 + lamb / 2.0 - lam34
    A[3, 7] = lamd / 2.0 + lamb / 2.0 - lam34
    A[3, 8] = -lamd / 2.0 + lamb / 2.0

    A[4, 1] = lamb / 4.0 - lamn / 4.0
    A[4, 2] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[4, 3] = lamb / 4.0 - lamn / 4.0
    A[4, 4] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[4, 5] = -lamd / 2.0 + lamb / 2.0
    A[4, 6] = -lamd / 2.0 + lamb / 2.0
    A[4, 7] = lamd / 2.0 + lamb / 2.0 - lam34
    A[4, 8] = lamd / 2.0 + lamb / 2.0 - lam34

    A[5, 5] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[5, 6] = -lamn / 4.0 + lam34 / 4.0
    A[5, 7] = lamn / 4.0 - lam34 / 4.0
    A[5, 8] = -lamn / 4.0 + lam34 / 4.0

    A[6, 5] = -lamn / 4.0 + lam34 / 4.0
    A[6, 6] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[6, 7] = -lamn / 4.0 + lam34 / 4.0
    A[6, 8] = lamn / 4.0 - lam34 / 4.0

    A[7, 5] = lamn / 4.0 - lam34 / 4.0
    A[7, 6] = -lamn / 4.0 + lam34 / 4.0
    A[7, 7] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[7, 8] = -lamn / 4.0 + lam34 / 4.0

    A[8, 5] = -lamn / 4.0 + lam34 / 4.0
    A[8, 6] = lamn / 4.0 - lam34 / 4.0
    A[8, 7] = -lamn / 4.0 + lam34 / 4.0
    A[8, 8] = lamn / 4.0 + 3.0 / 4.0 * lam34

    return A


def theta_transformation(lamd: float, lamn: float, lamb: float, theta: float = THETA) -> np.ndarray:
    '''Construct (I + θA)^-1, the matrix used in the θ-splitting solver.

    Parameters
    ----------
    lamd, lamn, lamb:
        Relaxation rates used when creating the MRT collision matrix.
    theta:
        Splitting parameter; default keeps parity with MIXLBM (0.5).

    Returns
    -------
    np.ndarray:
        Inverse transformation matrix applied to moments.'''
    
    lam0 = 0.0
    lam34 = 1.0

    # Matrix G is the discrete analog of the θ-transform used to advance
    # the linear collision operator with second-order accuracy.
    G = np.zeros((9, 9), dtype=np.float64)

    lamd_theta = lamd * theta
    lamn_theta = lamn * theta
    lamb_theta = lamb * theta
    lam34_theta = lam34 * theta

    def denom(*vals: float) -> float:
        prod = 1.0
        for value in vals:
            prod *= (value + 1.0)
        return prod

    G[0, 0] = 1.0 / (1.0 + lam0 * theta)
    coeff_01 = -theta * (lam0 - lamb) / ((theta * lamb) + 1.0) / (1.0 + lam0 * theta)
    for idx in (1, 2, 3, 4):
        G[0, idx] = coeff_01
    coeff_05 = -(
        -theta * lamb * lam34
        + 2.0 * theta * lam0 * lam34
        - theta * lam0 * lamb
        + lam0
        - 2.0 * lamb
        + lam34
    ) * theta / ((theta * lamb) + 1.0) / ((theta * lam34) + 1.0) / (1.0 + lam0 * theta)
    for idx in (5, 6, 7, 8):
        G[0, idx] = coeff_05

    common1 = denom(lamd_theta, lamb_theta, lamn_theta) * 4.0
    G[1, 1] = (
        theta * theta * lamd * lamn
        + theta * theta * lamb * lamd
        + 2.0 * theta * theta * lamn * lamb
        + 2.0 * theta * lamd
        + 3.0 * theta * lamn
        + 3.0 * theta * lamb
        + 4.0
    ) / common1
    G[1, 2] = -(lamb - lamn) * theta / ((theta * lamn) + 1.0) / ((theta * lamb) + 1.0) / 4.0
    G[1, 3] = (
        theta * lamd * lamn
        + theta * lamb * lamd
        - 2.0 * theta * lamn * lamb
        - lamn
        - lamb
        + 2.0 * lamd
    ) * theta / ((theta * lamn) + 1.0) / ((theta * lamb) + 1.0) / (theta * lamd + 1.0) / 4.0
    G[1, 4] = G[1, 2]
    G[1, 5] = (
        theta * lamd * lam34
        - 2.0 * theta * lamb * lamd
        + theta * lamb * lam34
        - lamd
        + 2.0 * lam34
        - lamb
    ) * theta / ((theta * lamb) + 1.0) / ((theta * lam34) + 1.0) / (theta * lamd + 1.0) / 2.0
    G[1, 6] = (lamd - lamb) * theta / ((theta * lamb) + 1.0) / (theta * lamd + 1.0) / 2.0
    G[1, 7] = G[1, 6]
    G[1, 8] = G[1, 5]

    G[2, 1] = G[1, 2]
    G[2, 2] = G[1, 1]
    G[2, 3] = G[1, 2]
    G[2, 4] = G[1, 3]
    G[2, 5] = G[1, 5]
    G[2, 6] = G[1, 5]
    G[2, 7] = G[1, 6]
    G[2, 8] = G[1, 6]

    G[3, 1] = G[1, 3]
    G[3, 2] = G[1, 2]
    G[3, 3] = G[1, 1]
    G[3, 4] = G[1, 2]
    G[3, 5] = G[1, 6]
    G[3, 6] = G[1, 5]
    G[3, 7] = G[1, 5]
    G[3, 8] = G[1, 6]

    G[4, 1] = G[1, 2]
    G[4, 2] = G[1, 3]
    G[4, 3] = G[1, 2]
    G[4, 4] = G[1, 1]
    G[4, 5] = G[1, 6]
    G[4, 6] = G[1, 6]
    G[4, 7] = G[1, 5]
    G[4, 8] = G[1, 5]

    base = ((theta * lamn) + 1.0) * ((theta * lam34) + 1.0) * 4.0
    G[5, 5] = (theta * lam34 + 3.0 * theta * lamn + 4.0) / base
    G[5, 6] = -(lam34 - lamn) * theta / base
    G[5, 7] = (lam34 - lamn) * theta / base
    G[5, 8] = G[5, 6]

    G[6, 5] = -G[5, 6]
    G[6, 6] = G[5, 5]
    G[6, 7] = -G[5, 6]
    G[6, 8] = G[5, 7]

    G[7, 5] = (lam34 - lamn) * theta / base
    G[7, 6] = -G[5, 6]
    G[7, 7] = G[5, 5]
    G[7, 8] = -G[5, 6]

    G[8, 5] = -G[5, 6]
    G[8, 6] = G[5, 7]
    G[8, 7] = -G[5, 6]
    G[8, 8] = G[5, 5]

    return G


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


def D_fick_diffusivity(m: float) -> float:
    """
    Fick diffusivity lookup used by alternative diffusion models.

    Returns
    -------
    float:
        Diffusivity value for species with molecular weight `m`.
    """
    table = np.array(
        [
            0.010000000000000,
            0.013645831365889,
            0.018620871366629,
            0.025409727055493,
            0.034673685045253,
            0.047315125896148,
            0.064565422903466,
            0.088104887300801,
            0.120226443461741,
            0.164058977319954,
            0.223872113856834,
            0.305492111321552,
            0.416869383470336,
            0.568852930843842,
            0.776247116628692,
            1.059253725177290,
            1.445439770745929,
            1.972422736114855,
            2.691534803926918,
            3.672823004980850,
        ],
        dtype=np.float64,
    )
    return table[4] * 0.2 / m


def NU_kinematic_viscosity(nN: int) -> float:
    """
    Return ν from the kinematic viscosity lookup table.

    Returns
    -------
    float:
        Kinematic viscosity corresponding to table index `nN`.
    """
    table = np.array(
        [
            1.0000,
            1.1708,
            1.3707,
            1.6048,
            1.8789,
            2.1998,
            2.5754,
            3.0153,
            3.5302,
            4.1331,
            4.8390,
            5.6654,
            6.6329,
            7.7657,
            9.0919,
            10.6446,
            12.4625,
            14.5908,
            17.0826,
            20.0000,
        ],
        dtype=np.float64,
    )
    return table[nN - 1]


def XI_bulk_viscosity() -> float:
    """
    Constant bulk viscosity ξ used in the MRT operator.

    Returns
    -------
    float:
        Constant bulk viscosity (0.4 in lattice units).
    """
    return 0.4


def solve_linear_system(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """
    Solve the dense linear system `matrix @ x = rhs`.

    Used for the θ-splitting correction that couples species fluxes.

    Returns
    -------
    np.ndarray:
        Solution vector x.
    """
    return np.linalg.solve(matrix, rhs)


class MaxwellStefanState:
    """
    Container for the lattice populations and macroscopic fields.

    Parameters
    ----------
    num_species:
        Number of chemical species tracked by the solver.
    nx, ny:
        Lattice dimensions.

    Attributes
    ----------
    fd, fd_new:
        Distribution functions f_σ,k before and after collision/streaming.
    rsigma, psigma:
        Densities and partial pressures per species.
    uxsigma, uysigma:
        Species velocity components stored for diagnostics / reuse.
    """

    def __init__(self, num_species: int, nx: int, ny: int):
        self.fd = np.zeros((num_species, 9, nx, ny), dtype=np.float64)
        self.fd_new = np.zeros_like(self.fd)
        self.rsigma = np.zeros((num_species, nx, ny), dtype=np.float64)
        self.psigma = np.zeros_like(self.rsigma)
        self.uxsigma = np.zeros_like(self.rsigma)
        self.uysigma = np.zeros_like(self.rsigma)


def initialise_from_partial_pressures(
    state: MaxwellStefanState,
    phi: np.ndarray,
    psigma_init: np.ndarray,
) -> None:
    """
    Populate a `MaxwellStefanState` from prescribed partial pressures.

    Parameters
    ----------
    state:
        Simulation state to initialise in-place.
    phi:
        Equation-of-state coefficients φ_σ.
    psigma_init:
        Partial pressure field with shape (species, nx, ny).

    Notes
    -----
    The routine mirrors MIXLBM's "simple initialization" by loading an equilibrium
    with zero velocities. It immediately computes the corresponding macroscopic
    fields for later reuse.
    """
    num_species, nx, ny = psigma_init.shape
    for s in range(num_species):
        # Convert partial pressures to densities using ρ_σ = 3 p_σ / φ_σ.
        rho_sigma = 3.0 * psigma_init[s] / phi[s]
        feq = equilibrium_distribution(1.0, phi[s], 0.0, 0.0)
        for i in range(nx):
            for j in range(ny):
                # Initialise populations with the local equilibrium and cache moments.
                feq_local = equilibrium_distribution(rho_sigma[i, j], phi[s], 0.0, 0.0)
                state.fd[s, :, i, j] = feq_local
                state.fd_new[s, :, i, j] = feq_local
                state.rsigma[s, i, j], state.uxsigma[s, i, j], state.uysigma[s, i, j] = hydrodynamic_moments(
                    feq_local
                )
                state.psigma[s, i, j] = phi[s] * state.rsigma[s, i, j] / 3.0


def maxwell_stefan_step(
    state: MaxwellStefanState,
    phi: np.ndarray,
    molecular_weights: np.ndarray,
    nN: int = 20,
    nB: int = 15,
    theta: float = THETA,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform one θ-splitting Maxwell–Stefan collision/streaming step.

    Parameters
    ----------
    state:
        Current lattice Boltzmann state (populations read from `state.fd`).
    phi:
        Equation-of-state factors φ_σ.
    molecular_weights:
        Species molecular weights M_σ.
    nN, nB:
        Indices into the viscosity and resistivity lookup tables (1-based,
        matching the Fortran code).
    theta:
        Splitting parameter; defaults to THETA (0.5).

    Returns
    -------
    tuple:
        rsigma_post, ux_mix, uy_mix where rsigma_post are the updated species
        densities and ux_mix/uy_mix are barycentric velocity components.
        Boundary handling is periodic; other boundaries need additional logic.

    Algorithmic outline (with equation cues)
    ---------------------------------------
    1) Neighbour sampling and moments: gather f at (i+ckx, j+cky) for all k,
       compute per-species (ρσ, uσ) and mixture fields (ρ, p). [standard LBM
       moments]
    2) Maxwell–Stefan drift pre-correction: build χ_{σν} and compute u*_σ =
       uσ + Σν χ_{σν} x_ν (uν − uσ), where x_ν are species fractions. [MS model]
    3) First θ half-step (collision in f-space):
         f* = f − θ A (feq(ρσ, u*_σ) − f)
       Transform A → Ã = A (I + θ A)^{-1} (via G = (I + θ A)^{-1}) and stream
       the BB(k) component. [θ-transform, MRT]
    4) Post-streaming moments (g⁺): compute ρσ, uσ from f_new.
    5) Maxwell–Stefan coupling: solve the linear system in species space for
       fluxes jσ using θ-corrected matrix (diagonal term 1 + θ λσ CHIσ − off-diag
       θ λσ xσ χ_{σν}). Update uσ = jσ / ρσ. [MS linear solve]
    6) Second θ half-step: map corrected moments back to populations using
       feq(ρσ, u*_σ) and the same transform G. [θ-transform]
    """
    num_species, _, nx, ny = state.fd.shape

    fd = state.fd
    fd_new = np.zeros_like(fd)

    # Temporary storage reused across lattice links to avoid reallocations.
    rsigma_local = np.zeros(num_species, dtype=np.float64)  # species densities ρ_σ
    psigma_local = np.zeros(num_species, dtype=np.float64)  # partial pressures p_σ
    uxsigma_local = np.zeros(num_species, dtype=np.float64)  # velocity components u_σx
    uysigma_local = np.zeros(num_species, dtype=np.float64)  # velocity components u_σy
    uxstar = np.zeros(num_species, dtype=np.float64)  # χ-corrected drift velocity x-component
    uystar = np.zeros(num_species, dtype=np.float64)  # χ-corrected drift velocity y-component
    lambda_sigma = np.zeros(num_species, dtype=np.float64)  # Maxwell–Stefan relaxation λ_σ
    lamd_sigma = np.zeros(num_species, dtype=np.float64)  # MRT λ_d (density mode)
    lamn_sigma = np.zeros(num_species, dtype=np.float64)  # MRT λ_n (shear mode)
    lamb_sigma = np.zeros(num_species, dtype=np.float64)  # MRT λ_b (bulk mode)
    CHI = np.zeros((num_species, num_species), dtype=np.float64)  # resistivity matrix χ_{σν}

    # Transport coefficients used for every lattice node in this call.
    nu_value = NU_kinematic_viscosity(nN)
    xi_value = XI_bulk_viscosity()

    for i in range(nx):
        for j in range(ny):
            for k in range(9):
                # Periodic neighbour reached by lattice direction k.
                iI = (i + D2Q9_CX[k]) % nx
                jI = (j + D2Q9_CY[k]) % ny

                # Recover neighbour macroscopic fields for each species.
                for s in range(num_species):
                    f_neighbor = fd[s, :, iI, jI]
                    rho, ux, uy = hydrodynamic_moments(f_neighbor)
                    rsigma_local[s] = rho
                    psigma_local[s] = phi[s] * rho / 3.0
                    uxsigma_local[s] = ux
                    uysigma_local[s] = uy

                # Guard against zero density/pressure when forming mixture values.
                r_mix = np.sum(rsigma_local)
                p_mix = np.sum(psigma_local)
                if r_mix <= 0.0:
                    r_mix = 1.0
                if p_mix <= 0.0:
                    p_mix = 1.0

                # Barycentric velocity from species contributions: u = Σσ(ρσ uσ)/ρ.
                ux_mix = np.dot(rsigma_local, uxsigma_local) / r_mix
                uy_mix = np.dot(rsigma_local, uysigma_local) / r_mix

                # Mixture molecular weight via harmonic mean of species masses.
                inv_MM_sum = 0.0
                for s in range(num_species):
                    inv_MM_sum += (rsigma_local[s] / r_mix) / molecular_weights[s]
                mixture_mass = 1.0 / inv_MM_sum if inv_MM_sum != 0.0 else 0.0

                # Maxwell–Stefan resistivity coefficients χ_{σν}
                # (see MIXLBM’s CHI build and PDF MS relations).
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

                # Drift velocities corrected by χ coupling: u*_σ.
                for s in range(num_species):
                    uxstar[s] = uxsigma_local[s]
                    uystar[s] = uysigma_local[s]
                    for vs in range(num_species):
                        factor = CHI[s, vs] * (rsigma_local[vs] / r_mix)
                        uxstar[s] += factor * (uxsigma_local[vs] - uxsigma_local[s])
                        uystar[s] += factor * (uysigma_local[vs] - uysigma_local[s])

                # Relaxation frequencies for MRT collision (lamd, lamn, lamb)
                # used in the matrix A and transform G.
                for s in range(num_species):
                    lambda_sigma[s] = p_mix * B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[s],
                        nB,
                    ) / r_mix
                    lamd_sigma[s] = lambda_sigma[s]
                    lamn_sigma[s] = 1.0 / (3.0 * nu_value)
                    lamb_sigma[s] = (2.0 - phi[s]) / (3.0 * xi_value)

                # θ-splitting first half-step (f → g) followed by streaming:
                # compute f*, then f_coll = f* + Ã (feq − f*), then stream BB(k).
                for s in range(num_species):
                    f_neighbor = fd[s, :, iI, jI]
                    feq = equilibrium_distribution(
                        rsigma_local[s],
                        phi[s],
                        uxstar[s],
                        uystar[s],
                    )
                    A = mrt_matrix(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s])
                    G = theta_transformation(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s], theta=theta)
                    delta = feq - f_neighbor
                    f_temp = f_neighbor - theta * (A @ delta)
                    A_tilde = A @ G
                    delta_tilde = feq - f_temp
                    f_coll = f_temp + A_tilde @ delta_tilde
                    fd_new[s, BB_OPPOSITE[k], i, j] = f_coll[BB_OPPOSITE[k]]

    rsigma_post = np.zeros((num_species, nx, ny), dtype=np.float64)
    uxsigma_post = np.zeros_like(rsigma_post)
    uysigma_post = np.zeros_like(rsigma_post)
    psigma_post = np.zeros_like(rsigma_post)

    # Reconstruct moments from post-streaming populations (g⁺ at node (i,j)).
    for i in range(nx):
        for j in range(ny):
            # Gather post-streaming moments for each species at node (i, j).
            for s in range(num_species):
                # Post-streaming moments per species (g^+ -> f^+).
                rho, ux, uy = hydrodynamic_moments(fd_new[s, :, i, j])
                rsigma_post[s, i, j] = rho
                uxsigma_post[s, i, j] = ux
                uysigma_post[s, i, j] = uy
                psigma_post[s, i, j] = phi[s] * rho / 3.0

            # Rebuild mixture properties at this node.
            r_mix = np.sum(rsigma_post[:, i, j])
            p_mix = np.sum(psigma_post[:, i, j])
            if r_mix <= 0.0:
                continue
            if p_mix <= 0.0:
                p_mix = 1.0

            # Mixture molecular weight and MRT frequencies get refreshed.
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
                lamd_sigma[s] = lambda_sigma[s]
                lamn_sigma[s] = 1.0 / (3.0 * nu_value)
                lamb_sigma[s] = (2.0 - phi[s]) / (3.0 * xi_value)

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

            # χ projections required by the MS linear solver (CHIsigma).
            CHIsigma = np.zeros(num_species, dtype=np.float64)
            for s in range(num_species):
                temp = 0.0
                for vs in range(num_species):
                    temp += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix)
                CHIsigma[s] = temp

            # Build the θ-corrected linear system for species fluxes (MS coupling):
            # diagonal: 1 + θ λσ CHIsigma; off-diagonal: − θ λσ xσ χ_{σν}.
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

            # Flux moments entering the linear solve: g_j = ρσ uσ (per axis).
            gjx = rsigma_post[:, i, j] * uxsigma_post[:, i, j]
            gjy = rsigma_post[:, i, j] * uysigma_post[:, i, j]

            jx = solve_linear_system(A_matrix, gjx)
            jy = solve_linear_system(A_matrix, gjy)

            for s in range(num_species):
                uxsigma_post[s, i, j] = jx[s] / rsigma_post[s, i, j] if rsigma_post[s, i, j] > 0 else 0.0
                uysigma_post[s, i, j] = jy[s] / rsigma_post[s, i, j] if rsigma_post[s, i, j] > 0 else 0.0

            # Re-evaluate drift velocities with the corrected fluxes (u*_σ again).
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

            # Second θ-half-step: map corrected moments back to populations using G.
            for s in range(num_species):
                feq = equilibrium_distribution(rsigma_post[s, i, j], phi[s], uxstar[s], uystar[s])
                A = mrt_matrix(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s])
                G = theta_transformation(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s], theta=theta)
                tempV = theta * (A @ feq)
                tempV1 = fd_new[s, :, i, j] + tempV
                fd_new[s, :, i, j] = G @ tempV1

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
