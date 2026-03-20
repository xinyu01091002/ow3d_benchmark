# Monochromatic Surface Potential Time Derivative From Bulk Quantities

This note derives the monochromatic surface potential time derivative from the bulk potential, following the same "bulk first, then evaluate at the free surface" logic used for the other surface quantities.

The target quantity is

```math
\widetilde{\Phi_t}(x,t)=\Phi_t(x,z=\eta(x,t),t).
```

We keep the same retained-harmonic convention:

- first order: keep only `\cos\theta`
- second order: keep only `\cos 2\theta`
- third order: keep only `\cos 3\theta`

## 1. Setup

Let

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

The retained second-order elevation is

```math
\eta^{(2)}=\varepsilon^2 G_2A_2\cos 2\theta,
```

where

```math
G_2
=
\frac{1}{2}hk\left(2+\cosh 2kh\right)\frac{\coth kh}{\sinh^2 kh}.
```

## 2. Free-surface expansion

Expanding `\Phi_t(x,\eta,t)` about `z=0`,

```math
\widetilde{\Phi_t}
=
\Phi_t|_{z=0}
+
\eta\Phi_{tz}|_{z=0}
+
\frac12\eta^2\Phi_{tzz}|_{z=0}
+
\cdots.
```

Collecting perturbation orders gives

```math
\widetilde{\Phi_t}^{(1)}=\Phi_t^{(1)}|_{z=0},
```

```math
\widetilde{\Phi_t}^{(2)}
=
\Phi_t^{(2)}|_{z=0}
+
\eta^{(1)}\Phi_{tz}^{(1)}|_{z=0},
```

```math
\widetilde{\Phi_t}^{(3)}
=
\Phi_t^{(3)}|_{z=0}
+
\eta^{(1)}\Phi_{tz}^{(2)}|_{z=0}
+
\eta^{(2)}\Phi_{tz}^{(1)}|_{z=0}
+
\frac12\left(\eta^{(1)}\right)^2\Phi_{tzz}^{(1)}|_{z=0}.
```

## 3. First-order explicit result

From `\Phi^{(1)}`,

```math
\Phi_t^{(1)}|_{z=0}
=
\varepsilon a\omega F_1\cosh(kh)\cos\theta.
```

Using `F_1=-\omega/(k\sinh kh)`,

```math
\widetilde{\Phi_t}^{(1)}
=
-\varepsilon a\frac{\omega^2}{k}\coth(kh)\cos\theta.
```

This is the retained first-harmonic result.

## 4. Second-order explicit result

### 4.1 Direct second-order term

```math
\Phi_t^{(2)}|_{z=0}
=
2\varepsilon^2\omega F_2\cosh(2kh)A_2\cos 2\theta
=
-\varepsilon^2 a^2\omega^2
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\cos 2\theta.
```

### 4.2 Surface-evaluation correction

Also,

```math
\Phi_{tz}^{(1)}|_{z=0}
=
\varepsilon a\omega F_1k\sinh(kh)\cos\theta
=
-\varepsilon a\omega^2\cos\theta.
```

Thus

```math
\eta^{(1)}\Phi_{tz}^{(1)}|_{z=0}
=
-\varepsilon^2 a^2\omega^2\cos^2\theta.
```

The retained second harmonic from this term is

```math
-\varepsilon^2 a^2\omega^2\frac12\cos 2\theta.
```

### 4.3 Final second-order result

Therefore

```math
\widetilde{\Phi_t}^{(2)}
=
-\varepsilon^2 a^2\omega^2
\left[
\frac12
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right]
\cos 2\theta.
```

## 5. Third-order explicit result

Now keep only the retained third harmonic.

### 5.1 Contribution from `\Phi_t^{(3)}|_{z=0}`

```math
\widetilde{\Phi_t}^{(3)}{}_A
=
\varepsilon^3 a^3
\left[
\frac{3}{64}
k\omega^2
\frac{\left(-11+2\cosh 2kh\right)\cosh 3kh}{\sinh^7(kh)}
\right]
\cos 3\theta.
```

### 5.2 Contribution from `\eta^{(1)}\Phi_{tz}^{(2)}|_{z=0}`

From the second-order potential,

```math
\Phi_{tz}^{(2)}|_{z=0}
=
4\varepsilon^2\omega kF_2\sinh(2kh)A_2\cos 2\theta
=
-3\varepsilon^2 a^2k\omega^2
\frac{\cosh(kh)}{\sinh^3(kh)}
\cos 2\theta.
```

Hence the retained third harmonic is

```math
\widetilde{\Phi_t}^{(3)}{}_B
=
-\varepsilon^3 a^3k\omega^2
\frac{3}{2}\frac{\cosh(kh)}{\sinh^3(kh)}
\cos 3\theta.
```

### 5.3 Contribution from `\eta^{(2)}\Phi_{tz}^{(1)}|_{z=0}`

Using `\Phi_{tz}^{(1)}|_{z=0}=-\varepsilon a\omega^2\cos\theta`, the retained third harmonic is

```math
\widetilde{\Phi_t}^{(3)}{}_C
=
-\varepsilon^3 a^3
\frac{k\omega^2G_2}{4h}
\cos 3\theta.
```

That is,

```math
\widetilde{\Phi_t}^{(3)}{}_C
=
-\varepsilon^3 a^3k\omega^2
\frac{\left(2+\cosh 2kh\right)\coth(kh)}{8\sinh^2(kh)}
\cos 3\theta.
```

### 5.4 Contribution from `\frac12(\eta^{(1)})^2\Phi_{tzz}^{(1)}|_{z=0}`

Since

```math
\Phi_{tzz}^{(1)}|_{z=0}
=
\varepsilon a\omega F_1k^2\cosh(kh)\cos\theta
=
-\varepsilon ak\omega^2\coth(kh)\cos\theta,
```

the retained third harmonic is

```math
\widetilde{\Phi_t}^{(3)}{}_D
=
-\varepsilon^3 a^3k\omega^2
\frac18\coth(kh)
\cos 3\theta.
```

### 5.5 Final third-order result

Collecting the retained third-harmonic pieces,

```math
\widetilde{\Phi_t}^{(3)}
=
\varepsilon^3 a^3 C_{\phi_t,3}\cos 3\theta,
```

with

```math
C_{\phi_t,3}
=
k\omega^2
\Bigg[
\frac{3}{64}
\frac{\left(-11+2\cosh 2kh\right)\cosh 3kh}{\sinh^7(kh)}
-\frac{3}{2}\frac{\cosh(kh)}{\sinh^3(kh)}
-\frac{\left(2+\cosh 2kh\right)\coth(kh)}{8\sinh^2(kh)}
-\frac18\coth(kh)
\Bigg].
```

## 6. Final retained-harmonic summary

The monochromatic retained-harmonic surface potential time derivative is

```math
\widetilde{\Phi_t}
\approx
\widetilde{\Phi_t}^{(1)}
+
\widetilde{\Phi_t}^{(2)}
+
\widetilde{\Phi_t}^{(3)},
```

with

```math
\widetilde{\Phi_t}^{(1)}
=
-\varepsilon a\frac{\omega^2}{k}\coth(kh)\cos\theta,
```

```math
\widetilde{\Phi_t}^{(2)}
=
-\varepsilon^2 a^2\omega^2
\left[
\frac12
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right]
\cos 2\theta,
```

```math
\widetilde{\Phi_t}^{(3)}
=
\varepsilon^3 a^3 C_{\phi_t,3}\cos 3\theta.
```

## 7. Compact implementation form

For later coding, write

```math
\widetilde{\Phi_t}^{(1)}=\varepsilon a\,C_{\phi_t,1}\cos\theta,
```

```math
\widetilde{\Phi_t}^{(2)}=\varepsilon^2 a^2\,C_{\phi_t,2}\cos 2\theta,
```

```math
\widetilde{\Phi_t}^{(3)}=\varepsilon^3 a^3\,C_{\phi_t,3}\cos 3\theta,
```

with

```math
C_{\phi_t,1}=-\frac{\omega^2}{k}\coth(kh),
```

```math
C_{\phi_t,2}
=
-\omega^2
\left[
\frac12
+
\frac{3\cosh(2kh)}{4\sinh^4(kh)}
\right],
```

and `C_{\phi_t,3}` given above.

This note is the direct bulk-derived target for the next implementation step if we want a more faithful `\phi_t` approximation than the simple `n\omega\mu_n` route.
