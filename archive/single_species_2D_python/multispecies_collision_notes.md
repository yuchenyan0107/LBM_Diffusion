# Multispecies BGK Collision in the Simple LBM Code

These lecture notes unpack how the multi-component helpers in
`functions/lbm_functions.py` work together:

- `compute_multispecies_moments`
- `equilibrium_multispecies`
- `collision_multispecies_BGK`

The structure mirrors the classic single-species D2Q9 derivation while
highlighting the extra bookkeeping needed once several species share the
same lattice.

## 1. Notation

- species index: $ \sigma = 1, \dots, N_s $
- velocity index: $ q = 0, \dots, 8 $
- distributions: $ F_{\sigma q}(x) = F[\sigma, q, ix, iy] $
- lattice velocities: $ \mathbf{c}_q = (c_{qx}, c_{qy}) $
- weights: $ w_q $
- species density: $ \rho_\sigma = \sum_q F_{\sigma q} $
- total density: $ \rho = \sum_\sigma \rho_\sigma $
- species velocity: $ \mathbf{u}_\sigma = (u_{\sigma x}, u_{\sigma y}) $
- barycentric velocity: $ \mathbf{u} = (u_x, u_y) $
- compressibility factor: $ \phi_\sigma $ (from `LocalEquilibrium.m`)

Throughout, $\sum_q$ means $\sum_{q=0}^8$ and $\sum_\sigma$ means
$\sum_{\sigma=1}^{N_s}$.

## 2. `compute_multispecies_moments`

**Goal:** recover species densities and velocities plus the mixture velocity.

### 2.1 Species densities

$$
\rho_\sigma(x) = \sum_q F_{\sigma q}(x)
$$

### 2.2 Species momenta and velocities

$$
\rho_\sigma(x)\,u_{\sigma x}(x) = \sum_q c_{qx} F_{\sigma q}(x), \qquad
\rho_\sigma(x)\,u_{\sigma y}(x) = \sum_q c_{qy} F_{\sigma q}(x)
$$

Velocities are obtained by dividing by $\rho_\sigma$ (with a small
regulariser for near-zero densities).

### 2.3 Mixture density and velocity

$$
\rho(x) = \sum_\sigma \rho_\sigma(x), \qquad
\mathbf{u}(x) =
\frac{\sum_\sigma \rho_\sigma(x)\,\mathbf{u}_\sigma(x)}
     {\sum_\sigma \rho_\sigma(x)}
$$

The function returns these quantities so downstream code does not need to
recompute them.

## 3. `equilibrium_multispecies`

**Goal:** build the discrete equilibrium $f^{\mathrm{eq}}_{\sigma q}$ for every
species.

The low-Mach Hermite expansion derived in `LocalEquilibrium.m` gives

$$
f^{\mathrm{eq}}_{\sigma 0} =
\frac{4}{9}\,\rho_\sigma
\left(\frac{9 - 5\phi_\sigma}{4} - \frac{3}{2}\lVert \mathbf{u} \rVert^{2}\right),
$$

$$
f^{\mathrm{eq}}_{\sigma q} =
w_q \rho_\sigma
\left(
  \phi_\sigma
  + 3\,\mathbf{c}_q \cdot \mathbf{u}
  + \frac{9}{2}(\mathbf{c}_q \cdot \mathbf{u})^{2}
  - \frac{3}{2}\lVert \mathbf{u} \rVert^{2}
\right), \qquad q = 1,\dots,8.
$$

Here $\lVert \mathbf{u} \rVert^{2} = u_x^{2} + u_y^{2}$ and
$\mathbf{c}_q \cdot \mathbf{u} = c_{qx} u_x + c_{qy} u_y$.
The barycentric velocity $\mathbf{u}$ is shared by all species, so the
equilibrium forces them to relax towards a common bulk motion while
preserving the species densities.

Implementation steps:

1. broadcast $\phi_\sigma$ over the lattice;
2. compute $\lVert \mathbf{u} \rVert^{2}$ once;
3. precompute $\mathbf{c}_q \cdot \mathbf{u}$;
4. fill the equilibrium array with the formulas above.

## 4. `collision_multispecies_BGK`

This routine wires everything together:

1. call `compute_multispecies_moments` to get the macroscopic fields;
2. call `equilibrium_multispecies` to build $f^{\mathrm{eq}}$;
3. relax the populations with BGK.

For each species and direction

$$
F_{\sigma q}^{\text{new}} =
F_{\sigma q}^{\text{old}}
- \frac{\Delta t}{\tau_\sigma}
\left(
  F_{\sigma q}^{\text{old}} - f^{\mathrm{eq}}_{\sigma q}
\right)
= (1 - \omega_\sigma) F_{\sigma q}^{\text{old}}
 + \omega_\sigma f^{\mathrm{eq}}_{\sigma q},
$$

with $\omega_\sigma = \Delta t / \tau_\sigma$. The code accepts either a
single relaxation time or one per species and broadcasts as needed.
Returning $\rho_\sigma$, $\rho$, $u_x$, and $u_y$ avoids recomputation in the
caller.

## 5. Workflow Summary

Inside the time step the calls look like

```python
rho_s, ux_s, uy_s, rho_tot, ux_mix, uy_mix = compute_multispecies_moments(F, cx, cy)
feq = equilibrium_multispecies(rho_s, phi, cx, cy, w, ux_mix, uy_mix)
F_new = F * (1 - omega) + feq * omega
```

Compared to the single-species case, two conceptual changes appear:

1. macroscopic fields are tracked per species before being combined;
2. the equilibrium uses species-specific $ \phi_\sigma $ while sharing the same
   barycentric velocity.

## 6. Maxwellian Derivation (Sketch)

Starting from the continuous Maxwellian

$$
f_\sigma^{(M)}(\mathbf{v}) =
\frac{\rho_\sigma}{\left(2\pi R T_\sigma\right)^{D/2}}
\exp\!\left[
  -\frac{(\mathbf{v} - \mathbf{u})^{2}}{2 R T_\sigma}
\right],
$$

one expands the exponential in Hermite polynomials and matches the first
few velocity moments with the discrete D2Q9 lattice. The factor
$\phi_\sigma = R T_\sigma$ plays the role of a species-dependent sound speed
squared. The MATLAB script `LocalEquilibrium.m` carries out this matching
symbolically and yields the coefficients used above.

## 7. Cheat Sheet

| Quantity | Expression | Implemented in |
|----------|------------|----------------|
| species density | $ \rho_\sigma = \sum_q F_{\sigma q} $ | `compute_multispecies_moments` |
| species velocity | $ \mathbf{u}_\sigma = (\sum_q \mathbf{c}_q F_{\sigma q}) / \rho_\sigma $ | `compute_multispecies_moments` |
| barycentric velocity | $ \mathbf{u} = \sum_\sigma \rho_\sigma \mathbf{u}_\sigma / \sum_\sigma \rho_\sigma $ | `compute_multispecies_moments` |
| equilibrium (rest) | $ f^{\mathrm{eq}}_{\sigma 0} = (4/9)\,\rho_\sigma\left((9-5\phi_\sigma)/4 - 1.5\,\lVert \mathbf{u} \rVert^{2}\right) $ | `equilibrium_multispecies` |
| equilibrium (moving) | $ f^{\mathrm{eq}}_{\sigma q} = w_q \rho_\sigma (\phi_\sigma + 3\,\mathbf{c}_q \cdot \mathbf{u} + 4.5 (\mathbf{c}_q \cdot \mathbf{u})^{2} - 1.5 \lVert \mathbf{u} \rVert^{2}) $ | `equilibrium_multispecies` |
| BGK update | $ F^{\text{new}} = F^{\text{old}} + \omega (f^{\mathrm{eq}} - F^{\text{old}}) $ | `collision_multispecies_BGK` |

By keeping the workflow parallel to the single-species solver, these
helpers make it easy to reason about multicomponent LBM while still
supporting Maxwell-Stefan style transport models.

