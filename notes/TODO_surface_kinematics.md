# TODO: Surface Kinematics Follow-Up

## Current status

The active surface-kinematics workflow in [`postprocess_ow3d_boundkinematics.m`](c:/Research/OW3D_benchmark/postprocess_ow3d_boundkinematics.m) is currently focused on:

- `u`
- `w`
- `phit`

The pressure-surrogate path has been intentionally removed for now to keep the code path clearer.

## Why pressure was removed for now

The previous pressure comparison depended strongly on the OW3D/reader definition of `phit`.

At the moment, the main unresolved issue is not the pressure algebra itself, but the target quantity used for `phit`.

So keeping pressure in the postprocessing script was making the workflow harder to interpret.

## Main suspicion

The OW3D reader currently allows two interpretations for the time derivative of the potential:

- `uncorrected`
  - `phit = Dt * phi`
- `sigma_corrected`
  - `phit = Dt * phi - w * sigma * etat`

The linear diagnostic in [`diagnose_ow3d_phit.m`](c:/Research/OW3D_benchmark/diagnose_ow3d_phit.m) suggests that:

- both versions have nearly the same first-harmonic shape as `-g eta`
- but the `sigma_corrected` version appears to amplify the first-harmonic magnitude too much at the surface
- the uncorrected `Dt * phi` appears closer to the expected linear scaling

So the current working suspicion is:

- the `sigma_corrected` formula for `phit` may be incorrect in sign, scaling, or interpretation for the present OW3D output

## Next recommended step

Before reintroducing any pressure-surrogate comparison:

1. confirm the intended meaning of OW3D kinematics-file `phi`
2. confirm whether `phit` should be interpreted as:
   - Eulerian `\Phi_t`
   - sigma-coordinate time derivative
   - or a transformed quantity requiring a different correction
3. only after that, reintroduce pressure based on the verified `phit`
