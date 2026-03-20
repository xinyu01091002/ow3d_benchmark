# Monochromatic Surface Pressure Surrogate: Explicit First-, Second-, and Third-Order Formulas

This note writes the monochromatic surface pressure surrogate in explicit form, in the same spirit as the earlier note for the surface horizontal velocity.

The goal here is not yet to build the full broadband implementation. The goal is to first write the clean monochromatic formulas that we would later generalize.

## 1. Quantity of interest

In the current OW3D postprocessing workflow, the surface pressure surrogate is

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

At the free surface, the dynamic boundary condition gives

```math
\Phi_t
+
\frac12\left(u^2+w^2\right)
+
g\eta
=
0
\qquad
\text{at } z=\eta.
```

Therefore

```math
\widetilde{p}(x,t)=g\eta(x,t).
```

So for the present quantity, once we evaluate at the nonlinear free surface, the surface pressure surrogate is exactly proportional to the free-surface elevation.

This is the key simplification.

## 2. Harmonic convention

As before, let

```math
\theta=\omega t-kx,
\qquad
\sigma=\tanh(kh).
```

We also keep the same retained-harmonic rule:

- first order: keep only `\cos\theta`
- second order: keep only `\cos 2\theta`
- third order: keep only `\cos 3\theta`

So the final formulas below are harmonic-truncated formulas.

## 3. First-order result

Take the first-order surface elevation as

```math
\eta^{(1)}=\varepsilon a\cos\theta.
```

Since `\widetilde{p}=g\eta`,

```math
\widetilde{p}^{(1)}
=
g\eta^{(1)}
=
\varepsilon g a\cos\theta.
```

This is the explicit first-order surface-pressure result.

## 4. Second-order result

For the monochromatic retained second harmonic, the surface elevation can be written as

```math
\eta^{(2)}
=
\varepsilon^2 a^2 k
\frac{3-\sigma^2}{8\sigma^3}
\cos 2\theta.
```

Using `\sigma=\tanh(kh)`, an equivalent hyperbolic-function form is

```math
\eta^{(2)}
=
\varepsilon^2
\frac{a^2 k}{4}
\left(2+\cosh 2kh\right)
\frac{\coth kh}{\sinh^2 kh}
\cos 2\theta.
```

Therefore the retained second-order surface pressure surrogate is

```math
\widetilde{p}^{(2)}
=
g\eta^{(2)}
=
\varepsilon^2 g a^2 k
\frac{3-\sigma^2}{8\sigma^3}
\cos 2\theta.
```

Equivalently,

```math
\widetilde{p}^{(2)}
=
\varepsilon^2
\frac{g a^2 k}{4}
\left(2+\cosh 2kh\right)
\frac{\coth kh}{\sinh^2 kh}
\cos 2\theta.
```

## 5. Third-order result

For the retained monochromatic third harmonic, a convenient explicit form is

```math
\eta^{(3)}
=
\varepsilon^3 a^3 k^2
\frac{27-9\sigma^2+9\sigma^4-3\sigma^6}{128\sigma^6}
\cos 3\theta.
```

Hence the retained third-order surface pressure surrogate is

```math
\widetilde{p}^{(3)}
=
g\eta^{(3)}
=
\varepsilon^3 g a^3 k^2
\frac{27-9\sigma^2+9\sigma^4-3\sigma^6}{128\sigma^6}
\cos 3\theta.
```

This is the explicit third-harmonic formula under the current retained-harmonic convention.

## 6. Final explicit summary

Putting the retained harmonics together,

```math
\widetilde{p}
\approx
\widetilde{p}^{(1)}
+
\widetilde{p}^{(2)}
+
\widetilde{p}^{(3)},
```

with

```math
\widetilde{p}^{(1)}
=
\varepsilon g a\cos\theta,
```

```math
\widetilde{p}^{(2)}
=
\varepsilon^2 g a^2 k
\frac{3-\sigma^2}{8\sigma^3}
\cos 2\theta,
```

```math
\widetilde{p}^{(3)}
=
\varepsilon^3 g a^3 k^2
\frac{27-9\sigma^2+9\sigma^4-3\sigma^6}{128\sigma^6}
\cos 3\theta.
```

Here again

```math
\sigma=\tanh(kh).
```

## 7. Why this matters for the code

These formulas show that, at the free surface, the present pressure surrogate is not an independent primitive quantity. It is the surface elevation multiplied by `g`.

So the natural workflow is:

1. derive or approximate `\eta^{(1)}`, `\eta^{(2)}`, `\eta^{(3)}`
2. map them to pressure with

```math
\widetilde{p}^{(n)} = g\,\eta^{(n)}
```

only after making sure the comparison target is truly the free-surface-evaluated quantity

```math
-\left(\Phi_t+\frac12(u^2+w^2)\right)\bigg|_{z=\eta}.
```

If the numerical target differs from this exact free-surface quantity, then discrepancies with the simple `g\eta` result are expected.
