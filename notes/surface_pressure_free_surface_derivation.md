# Surface Pressure From Bulk Quantities Evaluated At the Free Surface

This note explains the more complete perturbation-theory route for the surface pressure surrogate, consistent with the way we have been treating surface kinematics such as `u_s` and `w_s`.

The key idea is:

- do **not** start from a direct surface-pressure kernel
- first build the bulk quantities with their `z` dependence
- then evaluate them at the nonlinear free surface
- finally collect the retained harmonic contribution at each order

This is the same logic we used conceptually for `u_s` and `w_s`.

## 1. Quantity to be approximated

In the current OW3D postprocessing script, the pressure-like target is

```math
p(x,z,t)
=
-\left(\Phi_t + \frac12\left(u^2+w^2\right)\right),
```

with

```math
u=\Phi_x,
\qquad
w=\Phi_z.
```

The surface quantity of interest is therefore

```math
\widetilde{p}(x,t)
=
p(x,z=\eta(x,t),t).
```

So the correct object is

```math
\widetilde{p}(x,t)
=
-\left[
\Phi_t(x,\eta,t)
+
\frac12\left(
u(x,\eta,t)^2+w(x,\eta,t)^2
\right)
\right].
```

This already shows why a direct surface-`p` kernel is unnatural: `p` is not a primitive surface field, but a combination of other bulk quantities evaluated on the moving surface.

## 2. Perturbation expansions

Write the weakly nonlinear expansions

```math
\eta
=
\varepsilon\eta^{(1)}
+
\varepsilon^2\eta^{(2)}
+
\varepsilon^3\eta^{(3)}
+
\cdots,
```

```math
\Phi
=
\varepsilon\Phi^{(1)}
+
\varepsilon^2\Phi^{(2)}
+
\varepsilon^3\Phi^{(3)}
+
\cdots.
```

Then

```math
u
=
\Phi_x
=
\varepsilon u^{(1)}
+
\varepsilon^2 u^{(2)}
+
\varepsilon^3 u^{(3)}
+
\cdots,
```

```math
w
=
\Phi_z
=
\varepsilon w^{(1)}
+
\varepsilon^2 w^{(2)}
+
\varepsilon^3 w^{(3)}
+
\cdots,
```

```math
\Phi_t
=
\varepsilon \Phi_t^{(1)}
+
\varepsilon^2 \Phi_t^{(2)}
+
\varepsilon^3 \Phi_t^{(3)}
+
\cdots.
```

## 3. Evaluate bulk quantities at the nonlinear surface

For any bulk field `Q(x,z,t)`, the corresponding surface quantity is

```math
\widetilde{Q}(x,t)=Q(x,\eta(x,t),t).
```

To obtain an order-consistent perturbation expansion, expand around the mean surface `z=0`:

```math
Q(x,\eta,t)
=
Q(x,0,t)
+
\eta Q_z(x,0,t)
+
\frac12\eta^2 Q_{zz}(x,0,t)
+
\cdots.
```

Applying this to `u`, `w`, and `\Phi_t` gives

```math
\widetilde{u}
=
u|_{0}
+
\eta u_z|_{0}
+
\frac12\eta^2 u_{zz}|_{0}
+
\cdots,
```

```math
\widetilde{w}
=
w|_{0}
+
\eta w_z|_{0}
+
\frac12\eta^2 w_{zz}|_{0}
+
\cdots,
```

```math
\widetilde{\Phi_t}
=
\Phi_t|_{0}
+
\eta \Phi_{tz}|_{0}
+
\frac12\eta^2 \Phi_{tzz}|_{0}
+
\cdots.
```

Here `|_0` means evaluation at `z=0`.

## 4. Order-by-order surface expansion

Collecting terms in powers of `\varepsilon`, the surface-evaluated quantities become:

### First order

```math
\widetilde{u}^{(1)} = u^{(1)}|_0,
```

```math
\widetilde{w}^{(1)} = w^{(1)}|_0,
```

```math
\widetilde{\Phi_t}^{(1)} = \Phi_t^{(1)}|_0.
```

### Second order

```math
\widetilde{u}^{(2)}
=
u^{(2)}|_0
+
\eta^{(1)}u_z^{(1)}|_0,
```

```math
\widetilde{w}^{(2)}
=
w^{(2)}|_0
+
\eta^{(1)}w_z^{(1)}|_0,
```

```math
\widetilde{\Phi_t}^{(2)}
=
\Phi_t^{(2)}|_0
+
\eta^{(1)}\Phi_{tz}^{(1)}|_0.
```

### Third order

```math
\widetilde{u}^{(3)}
=
u^{(3)}|_0
+
\eta^{(1)}u_z^{(2)}|_0
+
\eta^{(2)}u_z^{(1)}|_0
+
\frac12(\eta^{(1)})^2u_{zz}^{(1)}|_0,
```

```math
\widetilde{w}^{(3)}
=
w^{(3)}|_0
+
\eta^{(1)}w_z^{(2)}|_0
+
\eta^{(2)}w_z^{(1)}|_0
+
\frac12(\eta^{(1)})^2w_{zz}^{(1)}|_0,
```

```math
\widetilde{\Phi_t}^{(3)}
=
\Phi_t^{(3)}|_0
+
\eta^{(1)}\Phi_{tz}^{(2)}|_0
+
\eta^{(2)}\Phi_{tz}^{(1)}|_0
+
\frac12(\eta^{(1)})^2\Phi_{tzz}^{(1)}|_0.
```

## 5. Surface pressure through third order

Now substitute the surface expansions into

```math
\widetilde{p}
=
-\left(
\widetilde{\Phi_t}
+
\frac12\left(\widetilde{u}^2+\widetilde{w}^2\right)
\right).
```

Then the perturbation orders are:

### First order

```math
\widetilde{p}^{(1)}
=
-\widetilde{\Phi_t}^{(1)}.
```

Since the quadratic velocity term starts at `O(\varepsilon^2)`, it does not enter at first order.

### Second order

```math
\widetilde{p}^{(2)}
=
-\left(
\widetilde{\Phi_t}^{(2)}
+
\frac12\left[
\left(\widetilde{u}^{(1)}\right)^2
+
\left(\widetilde{w}^{(1)}\right)^2
\right]
\right).
```

Equivalently, substituting the surface expansions,

```math
\widetilde{p}^{(2)}
=
-\left(
\Phi_t^{(2)}|_0
+
\eta^{(1)}\Phi_{tz}^{(1)}|_0
+
\frac12\left[
\left(u^{(1)}|_0\right)^2
+
\left(w^{(1)}|_0\right)^2
\right]
\right).
```

### Third order

```math
\widetilde{p}^{(3)}
=
-\left(
\widetilde{\Phi_t}^{(3)}
+
\widetilde{u}^{(1)}\widetilde{u}^{(2)}
+
\widetilde{w}^{(1)}\widetilde{w}^{(2)}
\right).
```

Substituting the explicit surface expansions,

```math
\widetilde{p}^{(3)}
=
-\Bigg(
\Phi_t^{(3)}|_0
+
\eta^{(1)}\Phi_{tz}^{(2)}|_0
+
\eta^{(2)}\Phi_{tz}^{(1)}|_0
+
\frac12(\eta^{(1)})^2\Phi_{tzz}^{(1)}|_0
+
u^{(1)}|_0
\left[
u^{(2)}|_0+\eta^{(1)}u_z^{(1)}|_0
\right]
+
w^{(1)}|_0
\left[
w^{(2)}|_0+\eta^{(1)}w_z^{(1)}|_0
\right]
\Bigg).
```

This is the consistent third-order surface-pressure expression under the present pressure definition.

## 6. Why this is better than a direct surface-pressure kernel

The direct-kernel idea for `p` is attractive computationally, but it hides the real structure of the quantity.

The expression for `p` contains:

- a time-derivative contribution
- a quadratic kinetic-energy contribution
- free-surface evaluation terms from the moving boundary

Those contributions do not naturally collapse into a single primitive surface-pressure kernel in the same simple way as some retained-harmonic `u_s` or `w_s` approximations.

So if the goal is consistency with the existing `u/w` workflow, the better route is:

1. derive `u(z)`, `w(z)`, and `\Phi_t(z)` order by order
2. evaluate them at the free surface
3. assemble `p` from those surface quantities

## 7. Harmonic truncation rule

If we apply the current retained-harmonic convention,

- first order: keep only the first harmonic
- second order: keep only the second harmonic
- third order: keep only the third harmonic

then after constructing `\widetilde{p}^{(1)}`, `\widetilde{p}^{(2)}`, and `\widetilde{p}^{(3)}` as above, we keep only the corresponding `n`th harmonic contribution at order `n`.

This means:

- lower-harmonic contamination created by Taylor-expansion cross terms is real
- but it can be intentionally discarded if the working convention is "order `n` keeps only harmonic `n`"

## 8. Practical implementation roadmap

For the code, this suggests the following workflow:

1. Build retained-harmonic kernels for the bulk quantities
   - `u(x,z,t)`
   - `w(x,z,t)`
   - `\Phi_t(x,z,t)`
2. Differentiate or evaluate those kernels with respect to `z` as needed
   - `u_z`, `u_{zz}`
   - `w_z`, `w_{zz}`
   - `\Phi_{tz}`, `\Phi_{tzz}`
3. Form
   - `\widetilde{u}^{(1:3)}`
   - `\widetilde{w}^{(1:3)}`
   - `\widetilde{\Phi_t}^{(1:3)}`
   at the free surface
4. Assemble
   - `\widetilde{p}^{(1)}`
   - `\widetilde{p}^{(2)}`
   - `\widetilde{p}^{(3)}`
5. Only then compare against the OW3D surface quantity

That is the consistent extension of the earlier surface-kinematics workflow to pressure.
