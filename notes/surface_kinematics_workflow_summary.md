# Surface Kinematics Comparison Workflow

This note summarizes the workflow we developed for comparing OW3D surface kinematics against a VWA-like approximation based on the linear surface elevation.

The immediate test case has been the surface horizontal velocity `u_s`, but the same overall workflow can be reused for `v_s` and other surface kinematic quantities.

## 1. Overall goal

The purpose is to compare:

- an OW3D reference quantity extracted from the kinematics output and separated into harmonic orders
- a reduced approximation driven only by the linear surface elevation `\eta^{(1)}`

For the current work, the target quantity is the physical surface horizontal velocity

```math
u_s = \left.\Phi_x\right|_{z=\eta}.
```

## 2. OW3D reference path

The OW3D reference is generated in [`postprocess_ow3d_boundkinematics.m`](c:/Research/OW3D_benchmark/postprocess_ow3d_boundkinematics.m).

The workflow is:

1. Read the four phase-shifted OW3D kinematics files.
2. Choose a valid time index.
3. Extract the requested field at all `sigma` levels.
4. Apply four-phase separation to reconstruct first-, second-, and third-order components.
5. Optionally apply x-direction harmonic cleanup around `n k_p`.
6. Select the surface layer using `sigma_mode = 'surface'`.

For `u`, this gives the OW3D reference harmonics:

- first order: `u^{(1)}`
- second order: `u^{(2)}`
- third order: `u^{(3)}`

These are treated as the comparison target.

## 3. Reduced approximation path

The reduced approximation uses only the linear surface elevation extracted from the OW3D four-phase result:

```math
\eta^{(1)}(x).
```

The idea is:

1. Treat `\eta^{(1)}` as the driving linear signal.
2. Build order-dependent transfer coefficients for the target surface quantity.
3. Filter `\eta^{(1)}` in Fourier space to obtain auxiliary fields `\kappa_n`.
4. Combine the filtered fields with powers of the analytic signal of `\eta^{(1)}`.
5. Compare the resulting approximate harmonics against the OW3D-separated harmonics.

For `u_s`, this was implemented in a VWA-like form inside the same postprocessing script.

## 4. Harmonic convention

The working convention has been:

- first order: keep only the first harmonic
- second order: keep only the second harmonic
- third order: keep only the third harmonic

That means:

- ignore the second-order mean term
- ignore the third-order first-harmonic correction
- focus only on the dominant `\cos n\theta` contribution at order `n`

This is a deliberate truncated harmonic representation, not the full asymptotic expansion.

## 5. What worked well

For the present `u_s` test:

- the OW3D four-phase separation itself became reliable after fixing the reader and Hilbert-storage issues
- the surface extraction path is straightforward once `sigma = 1` is selected
- the first-order comparison behaved well
- the second-order comparison became reasonable once the VWA-like product handling was simplified

So the basic workflow is sound.

## 6. Where the main difficulty appeared

The main difficulty showed up in the third-order approximation for `u_s`.

The important finding was:

- the strange third-order behavior was not mainly caused by the OW3D reader
- it was not simply a plotting issue
- it was not just a small phase shift

Instead, the critical issue was in the reduced-kernel path:

- the third-order kernel `U_3(k)` became extremely large at very small wavenumber
- this made `\kappa_3` dominated by near-zero wavenumber content
- when multiplied by `(\eta_+)^2`, the resulting `u^{(3)}` tended to peak near `2k_p` instead of `3k_p`

This explained why the third-order spectrum looked wrong.

## 7. Important implementation lessons

Several implementation details matter a lot:

1. The input to the reduced approximation must be the linear surface elevation only.

2. The physical comparison target must be the actual surface kinematic quantity, not a derivative of the surface potential unless that is the intended quantity.

3. The analytic-signal product structure must be handled carefully.

4. Small-wavenumber behavior of the kernel must be controlled.

5. Spectral diagnostics are essential.

In practice, plotting the spectra of:

- the OW3D harmonic
- the approximate harmonic
- the intermediate `\kappa_n`

was very helpful for diagnosing why the third-order result was misbehaving.

## 8. Why this workflow is reusable

The same workflow should apply to other surface kinematics such as:

- transverse velocity `v_s`
- vertical velocity `w_s`
- horizontal acceleration
- pressure-like surface quantities

The structure is always similar:

1. define the physical surface quantity to compare
2. extract its OW3D reference harmonic components
3. derive or specify reduced-order transfer coefficients driven by `\eta^{(1)}`
4. build the approximate harmonic fields
5. compare in both physical space and spectral space

So the workflow is not specific to `u_s`; only the transfer coefficients and the final product convention need to be adapted.

## 9. Recommended workflow for future quantities

For a new surface quantity such as `v_s`, the recommended workflow is:

1. Confirm the exact physical definition of the target quantity.
2. Verify that the OW3D postprocessor is extracting the correct surface field.
3. Reuse the existing four-phase separation pipeline.
4. Start with first and second order before trusting third order.
5. Add spectra from the beginning, not only line plots.
6. Check the small-`k` behavior of the reduced kernel early.
7. Only after the reduced approximation is stable should it be generalized or moved into a standalone helper.

## 10. Current status

At the moment:

- the OW3D-based surface harmonic extraction workflow is in good shape
- the `u_s` reduced approximation is useful as a testing framework
- the first and second harmonics are much more trustworthy than the current third-order reduced approximation
- the third-order mismatch has been narrowed down mainly to the reduced-kernel/analytic-product side

So the main outcome so far is not just a partial `u_s` model, but a reusable workflow for diagnosing and comparing reduced surface-kinematics models against OW3D references.

## 11. Strict second-order `u20` note

The repository now also contains a narrower diagnostic workflow in [`quick_extract_ow3d_subharmonic_velocity.m`](c:/Research/OW3D_benchmark/quick_extract_ow3d_subharmonic_velocity.m) for the strict second-order difference-frequency/subharmonic horizontal velocity.

That script does four things in sequence:

1. Read the four phase-shifted OW3D kinematics files with the full field order.
2. Read the matching `EP_XXXXX.bin` snapshot with the full reader and align the `EP` x-grid to the kinematics x-grid by removing the two extra points.
3. Reconstruct `\eta^{(1)}` from the extracted linear spectrum and verify that it matches the OW3D first-order `\eta`.
4. Build the strict Appendix-A difference-frequency prediction `u^{(2-)}` and compare it against the OW3D-extracted subharmonic `u20`.

The present diagnostic outcome is:

- the first-order `\eta` reconstruction is excellent, so the linear-spectrum extraction step looks trustworthy
- the strict Appendix-A `u^{(2-)}` implementation agrees with the existing MF12 direct difference-frequency path
- however, the OW3D surface `u20` does not behave like the bare Appendix-A difference-frequency velocity term

The main reason appears to be the OW3D surface-velocity definition itself. In the OW3D source, the stored horizontal velocity contains a sigma-coordinate chain-rule correction in addition to the constant-sigma `\phi_x` part. The decomposition tests currently indicate:

- the surface subharmonic `u20` is dominated by the chain-rule correction contribution
- the superharmonic `u2` is much closer to the bare constant-sigma `\phi_x` contribution

So for subharmonic velocity, comparing the OW3D surface `u20` directly against the bare Appendix-A `u^{(2-)}` is not a like-for-like comparison. Any future theory comparison should either:

- reconstruct the OW3D-style chain-rule contribution explicitly, or
- compare at the level of a quantity whose definition matches the Appendix-A derivation more directly
