# Surface Velocity Versus Surface Velocity Potential

This note summarizes the proper way to compute the surface horizontal velocity from a perturbation expansion, and explains why it is generally incorrect to simply differentiate the surface velocity potential and call that the horizontal velocity.

## 1. Two different objects

Let the bulk velocity potential be

```math
\Phi(x,y,z,t).
```

The free surface is

```math
z=\eta(x,y,t).
```

The surface velocity potential is the bulk potential evaluated on the moving free surface:

```math
\widetilde{\Phi}(x,y,t) = \Phi(x,y,z=\eta(x,y,t),t).
```

These are not the same object:

- `\Phi(x,y,z,t)` is a field in the fluid domain.
- `\widetilde{\Phi}(x,y,t)` is a scalar field living only on the free surface.

## 2. The correct definition of surface horizontal velocity

The horizontal velocity is defined from the bulk potential:

```math
u(x,y,z,t)=\frac{\partial \Phi}{\partial x}.
```

Therefore the surface horizontal velocity is

```math
\widetilde{u}(x,y,t)=\left.\frac{\partial \Phi}{\partial x}\right|_{z=\eta(x,y,t)}.
```

This is the quantity to compute if the goal is the physical horizontal velocity at the free surface.

## 3. Why differentiating the surface potential is not the same thing

If instead we differentiate the surface potential itself,

```math
\frac{\partial \widetilde{\Phi}}{\partial x},
```

then the chain rule gives

```math
\frac{\partial \widetilde{\Phi}}{\partial x}
=
\left.\frac{\partial \Phi}{\partial x}\right|_{z=\eta}
+
\left.\frac{\partial \Phi}{\partial z}\right|_{z=\eta}\frac{\partial \eta}{\partial x}.
```

Using

```math
\widetilde{u}=\left.\Phi_x\right|_{z=\eta}, \qquad
\widetilde{w}=\left.\Phi_z\right|_{z=\eta},
```

this becomes

```math
\widetilde{\Phi}_x = \widetilde{u} + \widetilde{w}\,\eta_x.
```

Hence

```math
\widetilde{u} = \widetilde{\Phi}_x - \widetilde{w}\,\eta_x.
```

So the derivative of the surface potential is not equal to the horizontal velocity, except in the special case where the correction term vanishes.

## 4. What the difference means physically

The difference is

```math
\widetilde{\Phi}_x - \widetilde{u} = \widetilde{w}\,\eta_x.
```

This term appears because `\widetilde{\Phi}` is a quantity restricted to the moving surface. When differentiating it with respect to `x`, one differentiates along the surface, not purely horizontally through the bulk field.

Therefore:

- `\widetilde{u}` is the physical horizontal velocity component at the free surface.
- `\widetilde{\Phi}_x` is the derivative of the surface-evaluated potential along the moving surface.

They coincide only if the surface slope or vertical velocity is negligible.

## 5. Order-of-magnitude interpretation

For a weakly nonlinear wave,

- `\eta_x = O(\varepsilon)`
- `\widetilde{w} = O(\varepsilon)`

so

```math
\widetilde{w}\,\eta_x = O(\varepsilon^2).
```

Therefore:

- at first order, `\widetilde{\Phi}_x` and `\widetilde{u}` are the same
- from second order onward, they differ

This is why using `\widetilde{\Phi}_x` as the surface horizontal velocity is acceptable only in a purely linear setting, but is wrong for second- or third-order surface kinematics.

## 6. Proper perturbation expansion for surface horizontal velocity

Write the expansions

```math
\Phi = \varepsilon \Phi^{(1)} + \varepsilon^2 \Phi^{(2)} + \varepsilon^3 \Phi^{(3)} + \cdots
```

```math
\eta = \varepsilon \eta^{(1)} + \varepsilon^2 \eta^{(2)} + \varepsilon^3 \eta^{(3)} + \cdots
```

Then

```math
\widetilde{u}
=
\left.\Phi_x\right|_{z=\eta}.
```

To evaluate this consistently, Taylor-expand about the mean water level `z=0`:

```math
\widetilde{u}
=
\left.\Phi_x\right|_{z=0}
+
\eta\left.\Phi_{xz}\right|_{z=0}
+
\frac{1}{2}\eta^2\left.\Phi_{xzz}\right|_{z=0}
+
\cdots
```

Collecting terms by order gives:

```math
\widetilde{u}^{(1)} = \Phi_x^{(1)}\big|_{z=0}
```

```math
\widetilde{u}^{(2)} = \Phi_x^{(2)}\big|_{z=0} + \eta^{(1)}\Phi_{xz}^{(1)}\big|_{z=0}
```

```math
\widetilde{u}^{(3)} =
\Phi_x^{(3)}\big|_{z=0}
+
\eta^{(1)}\Phi_{xz}^{(2)}\big|_{z=0}
+
\eta^{(2)}\Phi_{xz}^{(1)}\big|_{z=0}
+
\frac{1}{2}\left(\eta^{(1)}\right)^2\Phi_{xzz}^{(1)}\big|_{z=0}.
```

This is the proper third-order perturbation-consistent surface horizontal velocity.

## 7. Why this is the preferred method

This approach is correct because:

- the velocity is defined from the bulk potential, not from the surface potential
- the evaluation at the nonlinear surface is handled consistently through a Taylor expansion
- the resulting terms can be organized cleanly by order in `\varepsilon`

This is the standard way to derive surface kinematics in Stokes-type perturbation theory.

## 8. Practical implication for implementation

If the goal is to compute third-order surface horizontal velocity, the recommended workflow is:

1. Compute `\Phi^{(1)}`, `\Phi^{(2)}`, `\Phi^{(3)}`.
2. Compute `\eta^{(1)}`, `\eta^{(2)}` as needed.
3. Form `\widetilde{u}` from the Taylor expansion of `\Phi_x` evaluated at `z=0`.
4. Keep all terms up to the target order in `\varepsilon`.

The following shortcut should generally be avoided:

```math
\widetilde{u} \stackrel{\text{wrong}}{=} \frac{\partial \widetilde{\Phi}}{\partial x}.
```

It is wrong beyond first order because it omits the correction

```math
-\widetilde{w}\,\eta_x.
```

## 9. Summary

The correct surface horizontal velocity is

```math
\widetilde{u} = \left.\Phi_x\right|_{z=\eta},
```

not

```math
\widetilde{\Phi}_x.
```

Their relation is

```math
\widetilde{\Phi}_x = \widetilde{u} + \widetilde{w}\,\eta_x.
```

Thus:

- they are identical at linear order
- they differ from second order onward
- for third-order surface kinematics, the correct route is to differentiate the bulk potential first, then evaluate at the nonlinear surface through a Taylor expansion

## 10. The same idea for a general surface-evaluated quantity

The same distinction applies to any bulk quantity

```math
Q(x,y,z,t)
```

that is later evaluated on the free surface:

```math
\widetilde{Q}(x,y,t)=Q(x,y,z=\eta(x,y,t),t).
```

If one differentiates `\widetilde{Q}` with respect to space or time, the chain rule gives extra terms from the moving surface:

```math
\widetilde{Q}_x = \left.Q_x\right|_{z=\eta} + \left.Q_z\right|_{z=\eta}\eta_x,
```

```math
\widetilde{Q}_y = \left.Q_y\right|_{z=\eta} + \left.Q_z\right|_{z=\eta}\eta_y,
```

```math
\widetilde{Q}_t = \left.Q_t\right|_{z=\eta} + \left.Q_z\right|_{z=\eta}\eta_t.
```

Therefore, differentiating a surface-restricted quantity is generally not the same as differentiating the bulk quantity first and then evaluating it at the surface.

This is exactly the same issue as for the surface potential.

## 11. Pressure

Pressure is a bulk quantity. In potential flow it is typically obtained from the Bernoulli equation in the fluid interior:

```math
p = -\rho \left(\Phi_t + \frac{1}{2}\left(\Phi_x^2+\Phi_y^2+\Phi_z^2\right) + gz\right) + C.
```

If one wants pressure at the free surface, the correct object is

```math
\widetilde{p} = p(x,y,z=\eta,t).
```

For perturbation theory, this should be expanded about `z=0`, just as for the velocity:

```math
\widetilde{p}
=
\left.p\right|_{z=0}
+
\eta \left.p_z\right|_{z=0}
+
\frac{1}{2}\eta^2 \left.p_{zz}\right|_{z=0}
+
\cdots
```

The important point is that the pressure should be derived from the bulk potential first and then evaluated at the nonlinear surface consistently.

It is generally not correct to build a surface-pressure formula only from `\widetilde{\Phi}` and then differentiate that without accounting for the moving-surface correction terms.

## 12. Acceleration

The same principle applies to acceleration.

For example, if the horizontal velocity is

```math
u=\Phi_x,
```

then the horizontal local acceleration in the bulk is

```math
u_t = \Phi_{xt}.
```

The surface value is

```math
\widetilde{u_t} = \left.\Phi_{xt}\right|_{z=\eta}.
```

To obtain a perturbation-consistent expression, Taylor-expand about `z=0`:

```math
\widetilde{u_t}
=
\left.\Phi_{xt}\right|_{z=0}
+
\eta \left.\Phi_{xzt}\right|_{z=0}
+
\frac{1}{2}\eta^2 \left.\Phi_{xzzt}\right|_{z=0}
+
\cdots
```

Likewise for vertical acceleration or any other derived quantity.

Again, if one differentiates a surface quantity directly, extra chain-rule terms appear because the free surface is moving.

## 13. Why this matters beyond first order

For weakly nonlinear waves, the correction terms generated by the moving surface are often one order higher than the leading quantity.

This means:

- at first order, many surface and bulk-evaluated derivatives coincide
- at second and third order, the distinction becomes essential

Hence for second- or third-order kinematics, pressure, acceleration, or other surface quantities, the correct workflow is always:

1. define the physical quantity in the bulk
2. evaluate it at `z=\eta`
3. Taylor-expand about `z=0`
4. collect terms order by order in `\varepsilon`

## 14. Unified practical rule

If the target quantity is a physical field quantity such as:

- velocity
- pressure
- acceleration
- time derivative of potential

then the safe rule is:

```math
\text{differentiate first in the bulk, then evaluate at the nonlinear surface.}
```

If instead one starts from a surface-restricted quantity and differentiates it directly, then one is computing a derivative along the moving surface, and this generally differs from the physical field component evaluated at the surface.
