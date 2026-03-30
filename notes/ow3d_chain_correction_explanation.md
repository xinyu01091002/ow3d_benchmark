# OW3D Chain Correction In Surface Horizontal Velocity

This note explains what the "chain correction" means in the repository diagnostics, how it appears in the OW3D source code, and why it matters for the subharmonic `u` comparison in [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m).

## 1. Short version

In OW3D, the stored horizontal velocity is not just the simple Cartesian derivative `\phi_x` at fixed physical `z`.

Instead, OW3D computes velocity on a sigma-coordinate grid. Because a constant-`sigma` surface is tilted in physical space when the free surface or bathymetry varies, differentiating at fixed `sigma` introduces an extra chain-rule term involving the vertical velocity `w`.

In the present 2D surface case, this reduces to

```math
u = \phi_x|_\sigma - \eta_x\, w
```

when:

- the evaluation is at the free surface, so `\sigma = 1`
- the bed is flat locally, so `h_x = 0`

That extra term,

```math
-\eta_x\,w,
```

is what the scripts call the **chain correction**.

## 2. Where it appears in the diagnostic script

The relevant decomposition is in [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m).

First, the script computes the direct constant-`sigma` derivative of the stored potential:

```matlab
phix_sigma_phase(ip, :) = spectral_derivative_x_local(phi_phase_periodic(ip, :).', x_vec).';
```

from [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L916).

Then it defines the sigma-coordinate chain coefficient:

```matlab
chain_coeff_phase = ((1 - sigma_value) .* hx_value - sigma_value .* etax_phases);
```

from [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L921).

The chain contribution itself is:

```matlab
chain_phase = chain_coeff_phase .* w_phase_slice;
```

from [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L922).

Finally, the script removes this from the raw OW3D velocity to get the "bare" part:

```matlab
bare_phase = u_phase_slice - chain_phase;
```

from [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L923).

So the decomposition used by the script is:

```math
u_{\text{OW3D}} = u_{\text{bare}} + u_{\text{chain}}
```

with

```math
u_{\text{chain}} = \left((1-\sigma)h_x - \sigma \eta_x\right) w.
```

## 3. Where it appears in the OW3D source

The same term is present directly in the OW3D Fortran source that writes the kinematics output.

In [`StoreKinematicData.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/IO/StoreKinematicData.f90#L185), OW3D first computes the derivative of `phi` in `x`:

```fortran
CALL DiffXEven(phi,U,1,...)
```

Then it explicitly adds the chain-rule contribution:

```fortran
U(k,i,j) = U(k,i,j) + ((1-z(k))/d(i,j)*hx(i,j)-z(k)/d(i,j)*etax(i,j))*W(k,i,j)
```

from:

- [`StoreKinematicData.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/IO/StoreKinematicData.f90#L191)
- [`StoreKinematicData.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/IO/StoreKinematicData.f90#L196)

This is the direct source-code origin of the correction used in the MATLAB decomposition.

## 4. Why this formula has that form

OW3D uses a vertical coordinate transformation. The sigma-coordinate derivatives are built in [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90).

In particular:

```fortran
dsigma(k,i,j,2) = ((one-dsigma(k,i,j,1))*hx(i,j)-dsigma(k,i,j,1)*Ex(i,j))/d
```

from:

- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L68)
- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L70)

and in the 2D/3D general form:

```fortran
dsigma(k,i,j,2) = ((one-dsigma(k,i,j,1))*hx(i,j)-dsigma(k,i,j,1)*Ex(i,j))/d
dsigma(k,i,j,4) = ((one-dsigma(k,i,j,1))*hy(i,j)-dsigma(k,i,j,1)*Ey(i,j))/d
dsigma(k,i,j,5) = one/d
```

from:

- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L106)
- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L108)
- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L115)
- [`DetermineTransformationConstantsArray.f90`](c:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/src/functions/DetermineTransformationConstantsArray.f90#L117)

Here:

- `h_x` is the bed slope
- `E_x` in the OW3D source plays the role of free-surface slope `\eta_x`
- `d = h + \eta`
- `z(k)` is the sigma level value used in storage

So the correction term in the stored velocity is the metric term produced by the curvilinear `sigma` transformation.

## 5. Surface simplification

In the diagnostic script, the main comparison is usually done at the free surface:

- `sigma_mode = 'surface'`
- so `sigma_value = 1`

see [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L21) and [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L107).

Also, the current decomposition call uses

```matlab
hx_value = 0.0
```

see [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L122).

Therefore

```math
u_{\text{chain}} = \left((1-\sigma)h_x - \sigma \eta_x\right) w
```

reduces at the surface to

```math
u_{\text{chain}} = -\eta_x w.
```

So the "bare" horizontal velocity used in the comparison becomes

```math
u_{\text{bare}} = u_{\text{OW3D}} - (-\eta_x w)
= u_{\text{OW3D}} + \eta_x w.
```

Equivalently,

```math
u_{\text{OW3D}} = \phi_x|_\sigma - \eta_x w
```

for the present flat-bottom surface case.

## 6. Connection to the surface-potential chain rule

This is closely related to the usual free-surface identity discussed in [`surface_velocity_from_surface_potential_note.md`](c:/Research/OW3D_benchmark/notes/surface_velocity_from_surface_potential_note.md):

```math
\widetilde{\Phi}_x
=
\widetilde{u}
+
\widetilde{w}\eta_x,
```

so

```math
\widetilde{u} = \widetilde{\Phi}_x - \widetilde{w}\eta_x.
```

The OW3D sigma-coordinate correction is the grid-coordinate version of the same geometric idea:

- moving in `x` at fixed physical `z` is not the same as moving in `x` at fixed `\sigma`
- when the sigma surfaces tilt, part of the vertical motion contributes to the stored horizontal velocity

## 7. Why it matters for the subharmonic comparison

The diagnostic script shows that this correction is especially important for the subharmonic `u_{20}` comparison.

The raw OW3D subharmonic velocity is taken from:

- [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L113)
- [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L117)

but the script often compares theory against the "bare" velocity with the chain correction removed:

- [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L154)
- [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m#L227)

This is because the strict MF12 comparison target is closer to the direct potential-gradient part than to the full sigma-coordinate stored OW3D velocity.

The practical finding in this repository has been:

- the superharmonic part is often much closer to the bare `\phi_x|_\sigma`
- the subharmonic part is much more strongly contaminated by the chain correction

which is why removing the chain term is central to the `u_{20}` diagnostic workflow.

## 8. What this note does **not** imply for `w20`

The same sigma-coordinate logic does **not** automatically mean that the vertical velocity `w` should be treated the same way as `u`.

The present repository diagnostics indicate:

- the `u20` mismatch is strongly tied to the horizontal sigma-coordinate chain correction
- the `w20` mismatch is not explained by the same mechanism

This is consistent with the OW3D storage path:

- for `u`, OW3D writes the x-derivative plus the sigma metric correction term
- for `w`, OW3D writes the vertical velocity as the transformed vertical derivative divided by `d = h + \eta`

So:

- `u20` requires a bare/raw decomposition
- `w20` does not currently show evidence for an analogous bare/raw split driven by the same chain term

The current `w20` issue appears to be different:

- in deep water, the comparison is very sensitive to which sigma layer is used
- matching MF12 against the physical height of a deeper sigma layer improves the amplitude agreement substantially
- near the top sigma layer, the remaining `w20` discrepancy is therefore better interpreted as a near-surface vertical-structure issue than as a missing `u`-style chain correction

So this note should be used as the explanation for `u20`, but not as a complete explanation for `w20`.
