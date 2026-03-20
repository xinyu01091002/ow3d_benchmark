# Re-Derivation of the Third-Order Surface Horizontal Velocity Coefficient

This note re-derives the retained third-harmonic part of the monochromatic surface horizontal velocity from scratch, using the MF12 monochromatic formulas and the perturbation-consistent definition

```math
\widetilde{u}(x,t)=\left.\Phi_x\right|_{z=\eta}.
```

The purpose is narrow: check the algebra of the third-order `\cos 3\theta` coefficient.

We keep the same convention as before:

- first order: keep only `\cos\theta`
- second order: keep only `\cos 2\theta`
- third order: keep only `\cos 3\theta`
- ignore the third-order first-harmonic potential correction `F_{13}`

## 1. Monochromatic ingredients

Let

```math
\theta=\omega t-kx,
\qquad
Z=z+h.
```

Take

```math
\eta^{(1)}=\varepsilon a\cos\theta,
```

```math
\Phi^{(1)}=\varepsilon a F_1 \cosh(kZ)\sin\theta,
```

```math
\eta^{(2)}=\varepsilon^2 G_2 A_2 \cos 2\theta,
```

```math
\Phi^{(2)}=\varepsilon^2 F_2 A_2 \cosh(2kZ)\sin 2\theta,
```

```math
\Phi^{(3)}=\varepsilon^3 F_3 A_3 \cosh(3kZ)\sin 3\theta.
```

The monochromatic coefficients are

```math
F_1=-\frac{\omega}{k\sinh(kh)},
```

```math
G_2=\frac{1}{2}hk\left(2+\cosh 2kh\right)\frac{\coth(kh)}{\sinh^2(kh)},
```

```math
F_2=-\frac{3}{4}\frac{h\omega}{\sinh^4(kh)},
```

```math
F_3=\frac{1}{32}\frac{h^2 k\omega}{\sinh^7(kh)}\left(-11+2\cosh 2kh\right).
```

For the monochromatic self-interaction amplitudes:

```math
A_2=\frac{a^2}{2h},
\qquad
A_3=\frac{a^3}{2h^2}.
```

## 2. Third-order surface-velocity formula

The perturbation-consistent expansion is

```math
\widetilde{u}^{(3)}
=
\left.\Phi_x^{(3)}\right|_{z=0}
+
\eta^{(1)}\left.\Phi_{xz}^{(2)}\right|_{z=0}
+
\eta^{(2)}\left.\Phi_{xz}^{(1)}\right|_{z=0}
+
\frac{1}{2}\left(\eta^{(1)}\right)^2\left.\Phi_{xzz}^{(1)}\right|_{z=0}.
```

We now keep only the `\cos 3\theta` parts.

## 3. Term-by-term derivation

### 3.1 Direct third-order potential term

From

```math
\Phi^{(3)}=\varepsilon^3 F_3 A_3 \cosh(3kZ)\sin 3\theta,
```

we get

```math
\Phi_x^{(3)}
=
-3k\varepsilon^3 F_3 A_3 \cosh(3kZ)\cos 3\theta.
```

At `z=0`,

```math
\left.\Phi_x^{(3)}\right|_{z=0}
=
-3k\varepsilon^3 F_3 A_3 \cosh(3kh)\cos 3\theta.
```

Using `A_3=a^3/(2h^2)`, the contribution to the coefficient of `\varepsilon^3 a^3 \cos 3\theta` is

```math
C_{3,a}
=
-\frac{3k}{2h^2}F_3\cosh(3kh).
```

Substituting `F_3`:

```math
C_{3,a}
=
-\frac{3}{64}\omega k^2
\frac{-11+2\cosh(2kh)}{\sinh^7(kh)}
\cosh(3kh).
```

### 3.2 The term `\eta^{(1)} \Phi_{xz}^{(2)}`

First,

```math
\Phi_x^{(2)}
=
-2k\varepsilon^2 F_2 A_2 \cosh(2kZ)\cos 2\theta.
```

Hence

```math
\Phi_{xz}^{(2)}
=
-4k^2\varepsilon^2 F_2 A_2 \sinh(2kZ)\cos 2\theta.
```

At `z=0`,

```math
\left.\Phi_{xz}^{(2)}\right|_{z=0}
=
-4k^2\varepsilon^2 F_2 A_2 \sinh(2kh)\cos 2\theta.
```

Multiply by `\eta^{(1)}=\varepsilon a\cos\theta`:

```math
\eta^{(1)}\left.\Phi_{xz}^{(2)}\right|_{z=0}
=
-4k^2\varepsilon^3 a F_2 A_2 \sinh(2kh)\cos\theta\cos 2\theta.
```

Use

```math
\cos\theta\cos 2\theta=\frac{1}{2}\left(\cos\theta+\cos 3\theta\right).
```

So the `\cos 3\theta` part is

```math
-2k^2\varepsilon^3 a F_2 A_2 \sinh(2kh)\cos 3\theta.
```

Using `A_2=a^2/(2h)`, the coefficient contribution is

```math
C_{3,b}
=
-\frac{k^2 F_2}{h}\sinh(2kh).
```

Substituting `F_2`:

```math
C_{3,b}
=
\frac{3}{4}k^2\omega\frac{\sinh(2kh)}{\sinh^4(kh)}
=
\frac{3}{2}k^2\omega\frac{\cosh(kh)}{\sinh^3(kh)}.
```

### 3.3 The term `\eta^{(2)} \Phi_{xz}^{(1)}`

From the first-order potential,

```math
\Phi_x^{(1)}
=
-\varepsilon a k F_1 \cosh(kZ)\cos\theta,
```

so

```math
\Phi_{xz}^{(1)}
=
-\varepsilon a k^2 F_1 \sinh(kZ)\cos\theta.
```

At `z=0`,

```math
\left.\Phi_{xz}^{(1)}\right|_{z=0}
=
-\varepsilon a k^2 F_1 \sinh(kh)\cos\theta
=
\varepsilon a k\omega \cos\theta.
```

Multiply by

```math
\eta^{(2)}=\varepsilon^2 G_2 A_2 \cos 2\theta,
```

to get

```math
\eta^{(2)}\left.\Phi_{xz}^{(1)}\right|_{z=0}
=
\varepsilon^3 a k\omega G_2 A_2 \cos\theta\cos 2\theta.
```

Again using

```math
\cos\theta\cos 2\theta=\frac{1}{2}\left(\cos\theta+\cos 3\theta\right),
```

the `\cos 3\theta` part gives

```math
\frac{1}{2}\varepsilon^3 a k\omega G_2 A_2 \cos 3\theta.
```

Using `A_2=a^2/(2h)`, the coefficient contribution is

```math
C_{3,c}
=
\frac{k\omega G_2}{4h}.
```

Substituting `G_2`:

```math
C_{3,c}
=
\frac{1}{8}k^2\omega
\left(2+\cosh 2kh\right)\frac{\coth(kh)}{\sinh^2(kh)}
=
\frac{1}{8}k^2\omega
\frac{\left(2+\cosh 2kh\right)\cosh(kh)}{\sinh^3(kh)}.
```

### 3.4 The term `\frac{1}{2}(\eta^{(1)})^2 \Phi_{xzz}^{(1)}`

Differentiate once more:

```math
\Phi_{xzz}^{(1)}
=
-\varepsilon a k^3 F_1 \cosh(kZ)\cos\theta.
```

At `z=0`,

```math
\left.\Phi_{xzz}^{(1)}\right|_{z=0}
=
\varepsilon a \omega k^2 \coth(kh)\cos\theta.
```

Then

```math
\frac{1}{2}\left(\eta^{(1)}\right)^2\left.\Phi_{xzz}^{(1)}\right|_{z=0}
=
\frac{1}{2}\varepsilon^3 a^3 \omega k^2 \coth(kh)\cos^3\theta.
```

Use

```math
\cos^3\theta=\frac{3\cos\theta+\cos 3\theta}{4}.
```

Hence the `\cos 3\theta` part gives

```math
C_{3,d}
=
\frac{1}{8}\omega k^2 \coth(kh).
```

## 4. Final re-derived coefficient

Collecting the retained `\cos 3\theta` terms,

```math
\widetilde{u}^{(3)}_{\text{harmonic}}
=
\varepsilon^3 a^3 \widehat{C}_3 \cos 3\theta,
```

with

```math
\widehat{C}_3=C_{3,a}+C_{3,b}+C_{3,c}+C_{3,d},
```

i.e.

```math
\widehat{C}_3
=
-\frac{3}{64}\omega k^2
\frac{-11+2\cosh(2kh)}{\sinh^7(kh)}
\cosh(3kh)
+
\frac{3}{2}k^2\omega\frac{\cosh(kh)}{\sinh^3(kh)}
+
\frac{1}{8}k^2\omega
\frac{(2+\cosh 2kh)\cosh(kh)}{\sinh^3(kh)}
+
\frac{1}{8}\omega k^2 \coth(kh).
```

## 5. Conclusion of the algebra check

This independent re-derivation reproduces the same corrected `\widehat{C}_3` currently used in the code and in the updated monochromatic derivation note.

So if the third-order VWA-like comparison is still badly off, the remaining problem is unlikely to be the simple hand algebra in the monochromatic `\widehat{C}_3` expression itself. More likely candidates are:

- the way the monochromatic coefficient is being promoted to a broadband VWA-like kernel
- the analytic-signal convention used in the product form
- the fact that the OW3D third-order `u_s` for a packet is not represented well by this very truncated monochromatic-harmonic closure
