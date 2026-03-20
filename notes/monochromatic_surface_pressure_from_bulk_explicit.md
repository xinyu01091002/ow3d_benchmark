# Monochromatic Surface Pressure From Bulk Quantities: Explicit Retained-Harmonic Formulas

This note gives explicit monochromatic formulas for the surface pressure surrogate by following the same route we used conceptually for the surface kinematics:

1. start from bulk quantities with their `z` dependence
2. evaluate them at the free surface through a perturbation expansion
3. keep only the retained harmonic at each order

Here we do **not** start by collapsing the result to `g\eta`. Instead, we keep the bulk-derived structure explicit.

## 1. Definitions and setup

The pressure surrogate of interest is

```math
\widetilde{p}(x,t)
=
-\left[
\Phi_t(x,\eta,t)
+
\frac12\left(
u(x,\eta,t)^2+w(x,\eta,t)^2
\right)
\right],
```

with

```math
u=\Phi_x,
\qquad
w=\Phi_z.
```

For a monochromatic wave, let

```math
\theta=\omega t-kx,
\qquad
Z=z+h.
```

Take

```math
\eta^{(1)}=\varepsilon a\cos\theta.
```

The retained monochromatic potentials are

```math
\Phi^{(1)}=\varepsilon F_1\cosh(kZ)\,a\sin\theta,
```

```math
\Phi^{(2)}=\varepsilon^2 F_2\cosh(2kZ)\,A_2\sin 2\theta,
```

```math
\Phi^{(3)}=\varepsilon^3 F_3\cosh(3kZ)\,A_3\sin 3\theta,
```

with

```math
F_1=-\frac{\omega}{k\sinh(kh)},
\qquad
F_2=-\frac{3}{4}\frac{h\omega}{\sinh^4(kh)},
\qquad
F_3=\frac{1}{32}\frac{h^2k\omega}{\sinh^7(kh)}\left(-11+2\cosh 2kh\right),
```

and

```math
A_2=\frac{a^2}{2h},
\qquad
A_3=\frac{a^3}{2h^2}.
```

The retained second-order free-surface elevation is

```math
\eta^{(2)}
=
\varepsilon^2 G_2A_2\cos 2\theta,
```

with

```math
G_2
=
\frac{1}{2}hk\left(2+\cosh 2kh\right)\frac{\coth kh}{\sinh^2 kh}.
```

## 2. Order-by-order surface expansion

The consistent free-surface expansion is

```math
\widetilde{p}^{(1)}=-\Phi_t^{(1)}\big|_{z=0},
```

```math
\widetilde{p}^{(2)}
=
-\left(
\Phi_t^{(2)}\big|_{z=0}
+
\eta^{(1)}\Phi_{tz}^{(1)}\big|_{z=0}
+
\frac12\left[
\left(u^{(1)}\big|_{z=0}\right)^2
+
\left(w^{(1)}\big|_{z=0}\right)^2
\right]
\right),
```

```math
\widetilde{p}^{(3)}
=
-\Bigg(
\Phi_t^{(3)}\big|_{z=0}
+
\eta^{(1)}\Phi_{tz}^{(2)}\big|_{z=0}
+
\eta^{(2)}\Phi_{tz}^{(1)}\big|_{z=0}
+
\frac12\left(\eta^{(1)}\right)^2\Phi_{tzz}^{(1)}\big|_{z=0}
+
\widetilde{u}^{(1)}\widetilde{u}^{(2)}
+
\widetilde{w}^{(1)}\widetilde{w}^{(2)}
\Bigg).
```

The retained-harmonic convention is:

- first order: keep only `\cos\theta`
- second order: keep only `\cos 2\theta`
- third order: keep only `\cos 3\theta`

## 3. First-order explicit result

From the first-order potential,

```math
\Phi_t^{(1)}\big|_{z=0}
=
\varepsilon a\omega F_1\cosh(kh)\cos\theta.
```

Using `F_1=-\omega/(k\sinh kh)`,

```math
\Phi_t^{(1)}\big|_{z=0}
=
-\varepsilon a\frac{\omega^2}{k}\coth(kh)\cos\theta.
```

Therefore

```math
\widetilde{p}^{(1)}
=
\varepsilon a\frac{\omega^2}{k}\coth(kh)\cos\theta.
```

Using the dispersion relation `\omega^2=gk\tanh(kh)`, this becomes

```math
\widetilde{p}^{(1)}=\varepsilon ag\cos\theta.
```

## 4. Second-order explicit result

The retained second-order contributions are:

### 4.1 Time-derivative term

```math
\Phi_t^{(2)}\big|_{z=0}
=
2\varepsilon^2\omega F_2\cosh(2kh)A_2\cos 2\theta
=
-\varepsilon^2 a^2\omega^2
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\cos 2\theta.
```

### 4.2 Surface-evaluation correction

```math
\eta^{(1)}\Phi_{tz}^{(1)}\big|_{z=0}
=
-\varepsilon^2 a^2\omega^2\cos^2\theta.
```

The retained second harmonic from this term is

```math
-\varepsilon^2 a^2\omega^2\frac12\cos 2\theta.
```

### 4.3 Quadratic velocity term

At `z=0`,

```math
u^{(1)}\big|_{z=0}
=
\varepsilon a\omega\coth(kh)\cos\theta,
```

```math
w^{(1)}\big|_{z=0}
=
-\varepsilon a\omega\sin\theta.
```

Hence

```math
\frac12\left[
\left(u^{(1)}\big|_{z=0}\right)^2
+
\left(w^{(1)}\big|_{z=0}\right)^2
\right]
```

has retained second harmonic

```math
\varepsilon^2 a^2\omega^2
\frac{1}{4\sinh^2(kh)}
\cos 2\theta.
```

### 4.4 Final second-order coefficient

Collecting the second-harmonic terms gives

```math
\widetilde{p}^{(2)}
=
\varepsilon^2 a^2\omega^2
\left[
\frac12
-\frac{1}{4\sinh^2(kh)}
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right]
\cos 2\theta.
```

This is the explicit bulk-derived retained second-order result.

## 5. Third-order explicit result

Now keep only the retained third harmonic `\cos 3\theta`.

### 5.1 Contribution from `-\Phi_t^{(3)}|_{z=0}`

```math
\widetilde{p}^{(3)}_{A}
=
-\varepsilon^3 a^3 k\omega^2
\frac{3}{64}
\frac{\left(-11+2\cosh 2kh\right)\cosh 3kh}{\sinh^7(kh)}
\cos 3\theta.
```

### 5.2 Contribution from `-\eta^{(1)}\Phi_{tz}^{(2)}|_{z=0}`

```math
\widetilde{p}^{(3)}_{B}
=
\varepsilon^3 a^3 k\omega^2
\frac{3}{2}
\frac{\cosh(kh)}{\sinh^3(kh)}
\cos 3\theta.
```

### 5.3 Contribution from `-\eta^{(2)}\Phi_{tz}^{(1)}|_{z=0}`

Using the retained `\eta^{(2)}` term,

```math
\widetilde{p}^{(3)}_{C}
=
\varepsilon^3 a^3 k\omega^2
\frac{\left(2+\cosh 2kh\right)\coth(kh)}{8\sinh^2(kh)}
\cos 3\theta.
```

### 5.4 Contribution from `-\frac12(\eta^{(1)})^2\Phi_{tzz}^{(1)}|_{z=0}`

```math
\widetilde{p}^{(3)}_{D}
=
\varepsilon^3 a^3 k\omega^2
\frac{1}{8}\coth(kh)
\cos 3\theta.
```

### 5.5 Contribution from `-\widetilde{u}^{(1)}\widetilde{u}^{(2)}`

With

```math
\widetilde{u}^{(1)}
=
\varepsilon a\omega\coth(kh)\cos\theta,
```

```math
\widetilde{u}^{(2)}
=
\varepsilon^2 a^2k\omega
\left[
\frac12
+
\frac{3\cosh 2kh}{4\sinh^4 kh}
\right]
\cos 2\theta,
```

the retained third harmonic gives

```math
\widetilde{p}^{(3)}_{E}
=
-\varepsilon^3 a^3k\omega^2
\frac{\coth(kh)}{2}
\left[
\frac12
+
\frac{3\cosh 2kh}{4\sinh^4(kh)}
\right]
\cos 3\theta.
```

### 5.6 Contribution from `-\widetilde{w}^{(1)}\widetilde{w}^{(2)}`

With

```math
\widetilde{w}^{(1)}
=
-\varepsilon a\omega\sin\theta,
```

```math
\widetilde{w}^{(2)}
=
-\varepsilon^2 a^2k\omega
\left[
\frac12\coth(kh)
+
\frac{3}{2}\frac{\cosh(kh)}{\sinh^3(kh)}
\right]
\sin 2\theta,
```

the retained third harmonic gives

```math
\widetilde{p}^{(3)}_{F}
=
\varepsilon^3 a^3k\omega^2
\frac12
\left[
\frac12\coth(kh)
+
\frac32\frac{\cosh(kh)}{\sinh^3(kh)}
\right]
\cos 3\theta.
```

### 5.7 Final third-order coefficient

Putting the retained third-harmonic pieces together,

```math
\widetilde{p}^{(3)}
=
\varepsilon^3 a^3 C_{p,3}\cos 3\theta,
```

with

```math
C_{p,3}
=
k\omega^2
\Bigg[
-\frac{3}{64}
\frac{\left(-11+2\cosh 2kh\right)\cosh 3kh}{\sinh^7(kh)}
+
\frac{3}{2}\frac{\cosh(kh)}{\sinh^3(kh)}
+
\frac{\left(2+\cosh 2kh\right)\coth(kh)}{8\sinh^2(kh)}
+
\frac18\coth(kh)
-\frac{\coth(kh)}{2}
\left(
\frac12+\frac{3\cosh 2kh}{4\sinh^4(kh)}
\right)
+
\frac12
\left(
\frac12\coth(kh)+\frac32\frac{\cosh(kh)}{\sinh^3(kh)}
\right)
\Bigg].
```

This is explicit in `kh`, with no `F_n` or `G_n` placeholders left.

## 6. Final summary

The retained-harmonic monochromatic surface pressure surrogate can therefore be written as

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
\varepsilon a\frac{\omega^2}{k}\coth(kh)\cos\theta,
```

```math
\widetilde{p}^{(2)}
=
\varepsilon^2 a^2\omega^2
\left[
\frac12
-\frac{1}{4\sinh^2(kh)}
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right]
\cos 2\theta,
```

```math
\widetilde{p}^{(3)}
=
\varepsilon^3 a^3 C_{p,3}\cos 3\theta,
```

where `C_{p,3}` is the explicit coefficient written above.

## 6.1 Compact implementation form

For implementation, the same result can be written in the compact retained-harmonic form

```math
\widetilde{p}^{(1)}=\varepsilon a\,C_{p,1}\cos\theta,
```

```math
\widetilde{p}^{(2)}=\varepsilon^2 a^2\,C_{p,2}\cos 2\theta,
```

```math
\widetilde{p}^{(3)}=\varepsilon^3 a^3\,C_{p,3}\cos 3\theta,
```

with

```math
C_{p,1}=\frac{\omega^2}{k}\coth(kh),
```

```math
C_{p,2}
=
\omega^2
\left[
\frac12
-\frac{1}{4\sinh^2(kh)}
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right],
```

and

```math
C_{p,3}
=
k\omega^2
\Bigg[
-\frac{3}{64}
\frac{\left(-11+2\cosh 2kh\right)\cosh 3kh}{\sinh^7(kh)}
+
\frac{3}{2}\frac{\cosh(kh)}{\sinh^3(kh)}
+
\frac{\left(2+\cosh 2kh\right)\coth(kh)}{8\sinh^2(kh)}
+
\frac18\coth(kh)
-\frac{\coth(kh)}{2}
\left(
\frac12+\frac{3\cosh 2kh}{4\sinh^4(kh)}
\right)
+
\frac12
\left(
\frac12\coth(kh)+\frac32\frac{\cosh(kh)}{\sinh^3(kh)}
\right)
\Bigg].
```

So if we later code this in the same style as the other surface quantities, the real task is just to evaluate `C_{p,1}`, `C_{p,2}`, and `C_{p,3}` as functions of `kh`.

## 7. Interpretation

This note gives the explicit bulk-derived retained-harmonic formulas before making any additional assumption about whether the numerical quantity being compared is exactly identical to the ideal Bernoulli free-surface quantity.

So this is the right starting point if we want to:

- stay consistent with the earlier `u_s` and `w_s` workflow
- keep the `z` dependence conceptually present
- and only later decide what simplifications are justified when comparing with the postprocessed OW3D quantity
