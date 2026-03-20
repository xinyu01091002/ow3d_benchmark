# Third-order theory for multi-directional irregular waves

Per A. Madsen† and David R. Fuhrman 

Department of Mechanical Engineering, Technical University of Denmark, DK-2800 Kgs Lyngby, Denmark 

(Received 1 July 2011; revised 31 October 2011; accepted 9 February 2012; first published online 2 April 2012) 

A new third-order solution for multi-directional irregular water waves in finite water depth is presented. The solution includes explicit expressions for the surface elevation, the amplitude dispersion and the vertical variation of the velocity potential. Expressions for the velocity potential at the free surface are also provided, and the formulation incorporates the effect of an ambient current with the option of specifying zero net volume flux. Harmonic resonance may occur at third order for certain combinations of frequencies and wavenumber vectors, and in this situation the perturbation theory breaks down due to singularities in the transfer functions. We analyse harmonic resonance for the case of a monochromatic short-crested wave interacting with a plane wave having a different frequency, and make long-term simulations with a high-order Boussinesq formulation in order to study the evolution of wave trains exposed to harmonic resonance. 

Key words: surface gravity waves, waves/free-surface flows 

# 1. Introduction

Nonlinear irregular multi-directional water waves are commonly described to second order using the formulation by Sharma & Dean (1981). This is based on a double summation over all possible pairs of wave components, utilizing a second-order solution for bi-directional bichromatic waves as the kernel in the summation. In this work, we present a third-order analytical solution for trichromatic tri-directional waves, which is then used as the kernel in a triple summation over all triplets, plus a double summation over the relevant doublets. The outcome of this work is thus a third-order theory for multi-directional irregular waves in finite water depth, including the effect of ambient or wave-induced currents. The formulation is an extension of the work by Zhang & Chen (1999) from deep water to arbitrary depth, and from collinear interactions to multi-directional interactions. It is also an extension of the work by Madsen & Fuhrman (2006) from third-order bi-directional bichromatic waves to multi-directional irregular waves. 

Fuhrman & Madsen (2006) demonstrated the importance of accounting for third-order effects in boundary conditions for monochromatic short-crested waves in connection with numerical or laboratory experiments. They made numerical simulations using a high-order Boussinesq-type formulation and performed an analysis 

of the consequence of neglecting the second and third harmonics at the boundary. On this basis they concluded that, for wave fields with two-dimensional surface patterns, the neglect of third-harmonic terms at the boundary can cause a release of spurious free first harmonics, which lead to unsteady features such as slow modulations of amplitudes in the main wave direction in addition to curving of crestlines with local dips and peaks in the transverse direction. These unsteady phenomena were found to be at least as significant as the well-known spurious features resulting from the neglect of second-harmonic terms. Their conclusions were experimentally confirmed by Henderson, Patterson & Segur (2006), even though the third-order boundary conditions for their paddle motion were only approximate and without any compensation for evanescent modes. Fuhrman & Madsen (2006) and Henderson et al. (2006) found that the third-order short-crested theory was useful for providing steady solutions even for grazing angles as small as $1 0 ^ { \circ }$ and wave steepness as high as $k a = 0 . 3 0$ . The new theory presented in this paper provides a method for nonlinear analysis of a general sea state, and may also be useful for the generation of more complex wave trains involving several wave components of different frequencies and directions. 

Indeed, phenomena related to short-crested waves have been extensively studied in the literature with theoretical, numerical and experimental methods. The reason is clearly that this case is of fundamental importance to the understanding of multidirectional irregular waves, and at the same time it has a relatively simple form and mathematical description even at high orders of nonlinearity. Experimentally, shortcrested waves have been studied by, for example, Hammack, Scheffner & Segur (1989), Kimmoun, Ioualalen & Kharif (1999), Hammack, Henderson & Segur (2005) and Henderson et al. (2006). Analytical solutions have been provided by Chappelear (1961) and Hsu, Tsuchiya & Silvester (1979) to third order and by Ioualalen (1993) to fourth order. Roberts (1983), Roberts & Schwartz (1983) and Marchant & Roberts (1987) computed short-crested waves using a numerical perturbation method up to 27th and 35th order for deep water and finite depth, respectively. Roberts & Peregrine (1983) treated the important limit of grazing angles, where the short-crested deepwater waves become long-crested. Fully numerical computations have been performed by, for example, Bryant (1985), Craig & Nicholls (2002) and Fuhrman & Madsen (2006). In the present work we study the interaction of a monochromatic short-crested wave with a plane wave of another frequency, which is of similar fundamental interest, being the simplest extension beyond the monochromatic short-crested case. For this purpose we use a high-order Boussinesq-type model with third-order boundary conditions determined by the new theory. 

One of the important properties of short-crested waves is the occurrence of harmonic resonance. This phenomenon was first discussed and defined by Phillips (1960), and he concluded that it occurs whenever the frequency of a bound higher harmonic becomes identical to the frequency of a free infinitesimal wave of the same wavenumber. While non-resonating components are bounded in time and receive only a limited pulse of energy from the primary components, the resonating components can receive a continuing energy transfer until they reach the level of the primary components. This is the key mechanism in the long-term evolution of wave spectra in the ocean, and the seminal work by Phillips (1960) had a major impact on the formulation of wave evolution models, such as Hasselmann (1962), Zakharov (1968), Krasitskii (1994) and Janssen (2003). 

Roberts (1981, 1983) identified the resonating higher harmonics for short-crested waves, and he found the resonances to be densely distributed but of relatively high 

order (the lowest order was shown to be sixth order). Ioualalen & Kharif (1993, 1994) investigated the stability of deep-water short-crested waves in order to estimate the time scales of these resonances, and they identified the associated instabilities as McLean’s (1982) class I. However, they found them to be very sporadic and weak, with a long time scale evolution. Other stability aspects of short-crested waves were investigated by, for example, Badulin et al. (1995), Ioualalen, Roberts & Kharif (1996), Ioualalen, Kharif & Roberts (1999) and Fuhrman, Madsen & Bingham (2006). Smith & Roberts (1999) and Ioualalen et al. (2006) concluded that harmonic resonance manifests in the occurrence of multiple solutions (non-uniqueness) and associated bifurcation behaviour. 

For more complicated wave patterns involving at least two different frequencies, harmonic resonance can occur at third order, as shown by Phillips (1960). As an example, Phillips (1960) and Longuet-Higgins (1962) considered the case of two interacting bichromatic bi-directional waves, which are resonating for certain wavenumber combinations. They predicted the linear growth rate of the resonating component, and a simple laboratory experiment with two perpendicular primary waves was proposed in order to demonstrate and investigate this phenomenon. In this work, we present a similar, but slightly more complicated, numerical experiment with a combination of a short-crested wave interacting with a plane wave of another frequency. Our objective is to study the evolution of the wave pattern in the vicinity of third-order resonance. 

Finally, we acknowledge that the Zakharov formulation, e.g. as given by Zakharov (1968, 1999) and Krasitskii (1994), in principle allows for a determination of thirdorder expressions for the surface elevation and surface potential of steady trichromatic tri-directional waves in finite depth. However, until recently, it has been a problem to determine a unique solution to the Zakharov kernel function in finite depth, an issue only recently resolved by Janssen & Onorato (2007) for the special case of $T ( k _ { 1 } , k _ { 1 } , k _ { 1 } , k _ { 1 } )$ and by Stiassnie & Gramstad (2009) for the case of $T ( k _ { 1 } , k _ { 2 } , k _ { 1 } , k _ { 2 } )$ . Rather than pursuing a theory based on the Zakharov formulation, the present work is founded on classical Stokes-type perturbation expansions. 

The paper is organized as follows. In $\ S 2$ we summarize the governing equations and the perturbation method necessary to obtain higher-order analytical solutions for irregular waves. The starting point is a linear wave group consisting of a superposition of at least three different wave components (a trichromatic wave with different frequencies and different wavenumber vectors). The classical first- and second-order solutions to this problem are given in $\ S \ S 3 . 1$ and 3.2, respectively. The new third-order solution for irregular multi-directional waves is given in $\ S 3 . 3$ . This solution is defined in terms of the surface elevation and the vertical variation of the velocity potential, including general relations for amplitudes, wavenumbers, frequencies, transfer functions and amplitude dispersion. In $\ S 3 . 4$ the velocity potential is evaluated directly on the free surface, and this may be used to provide boundary conditions for numerical models. In $\ S 3 . 5$ we illustrate the analytical solution for a particular case of a short-crested wave interacting with a plane wave of a different frequency. Then, $\ S 4$ is devoted to the discussion of harmonic resonance. First, the resonance curves are identified and depicted for a monochromatic wave interacting with a plane wave of another frequency. Second, long-term simulations are made with a high-order Boussinesq-type formulation to investigate the evolution of wave fields exposed to harmonic resonance. Finally, $\ S 5$ contains a summary and conclusions. 

# 2. The governing equations and the perturbation method

We consider the classical problem of an irrotational flow of an incompressible inviscid fluid with a free surface and a horizontal bottom. A Cartesian coordinate system is adopted with the $x$ -axis and the $y$ -axis located on the mean water plane (MWP) and with the z-axis pointing vertically upwards. Hence, the domain is bounded by the horizontal sea bed at $z = - h$ and by the free surface $z = \eta ( x , y , t )$ . The irrotationality of the flow is expressed through the introduction of the velocity potential $\phi$ defined by 

$$
u (x, y, z, t) \equiv \frac {\partial \Phi}{\partial x}, \quad v (x, y, z, t) \equiv \frac {\partial \Phi}{\partial y}, \quad w (x, y, z, t) \equiv \frac {\partial \Phi}{\partial z}, \tag {2.1}
$$

where u, v and $w$ are the components of the particle velocity in the $x , \ y$ and z directions, respectively. In this framework, the governing equations for the fully nonlinear wave problem include two linear equations (the Laplace equation and the kinematic bottom condition), 

$$
\frac {\partial^ {2} \Phi}{\partial x ^ {2}} + \frac {\partial^ {2} \Phi}{\partial y ^ {2}} + \frac {\partial^ {2} \Phi}{\partial z ^ {2}} = 0 \quad \text {f o r} - h \leqslant z \leqslant \eta (x, y, t), \tag {2.2}
$$

$$
\frac {\partial \Phi}{\partial z} = 0 \quad \text {a t} z = - h, \tag {2.3}
$$

and two nonlinear equations (the kinematic and dynamic surface conditions), 

$$
\frac {\partial \eta}{\partial t} - \widetilde {w} + \widetilde {u} \frac {\partial \eta}{\partial x} + \widetilde {v} \frac {\partial \eta}{\partial y} = 0 \quad \text {a t} z = \eta (x, y, t), \tag {2.4}
$$

$$
\widetilde {\Psi} + g \eta + \frac {1}{2} \left(\widetilde {u} ^ {2} + \widetilde {v} ^ {2} + \widetilde {w} ^ {2}\right) = 0 \quad \text {a t} z = \eta (x, y, t). \tag {2.5}
$$

Note that (2.4) and (2.5) have been expressed in terms of surface variables defined by 

$$
\widetilde {u} \equiv u (x, y, \eta , t) \equiv \left(\frac {\partial \Phi}{\partial x}\right) _ {z = \eta}, \quad \widetilde {v} \equiv v (x, y, \eta , t) \equiv \left(\frac {\partial \Phi}{\partial y}\right) _ {z = \eta}, \tag {2.6}
$$

$$
\widetilde {w} \equiv w (x, y, \eta , t) \equiv \left(\frac {\partial \Phi}{\partial z}\right) _ {z = \eta}, \quad \widetilde {\Psi} \equiv \left(\frac {\partial \Phi}{\partial t}\right) _ {z = \eta}. \tag {2.7}
$$

We emphasize that, in connection with the perturbation method introduced in the following, it is necessary to express the velocity variables at the free surface in terms of Taylor series expansions from the mean water datum $z = 0$ . In this process, we include the first three terms in the Taylor expansions and obtain 

$$
\widetilde {u} \simeq \left(\frac {\partial \Phi}{\partial x} + \eta \frac {\partial^ {2} \Phi}{\partial x \partial z} + \frac {1}{2} \eta^ {2} \frac {\partial^ {3} \Phi}{\partial x \partial z ^ {2}}\right) _ {z = 0}, \quad \widetilde {v} \simeq \left(\frac {\partial \Phi}{\partial y} + \eta \frac {\partial^ {2} \Phi}{\partial y \partial z} + \frac {1}{2} \eta^ {2} \frac {\partial^ {3} \Phi}{\partial y \partial z ^ {2}}\right) _ {z = 0}, \tag {2.8}
$$

$$
\widetilde {w} \simeq \left(\frac {\partial \Phi}{\partial z} + \eta \frac {\partial^ {2} \Phi}{\partial z ^ {2}} + \frac {1}{2} \eta^ {2} \frac {\partial^ {3} \Phi}{\partial z ^ {3}}\right) _ {z = 0}, \quad \widetilde {\Psi} \simeq \left(\frac {\partial \Phi}{\partial t} + \eta \frac {\partial^ {2} \Phi}{\partial t \partial z} + \frac {1}{2} \eta^ {2} \frac {\partial^ {3} \Phi}{\partial t \partial z ^ {2}}\right) _ {z = 0}. \tag {2.9}
$$

The analytical solutions derived in the following chapters are of the form 

$$
\eta (x, y, t) = \eta^ {(1)} + \eta^ {(2)} + \eta^ {(3)} + \dots , \tag {2.10}
$$

$$
\Phi (x, y, z, t) = \Phi^ {(1)} + \Phi^ {(2)} + \Phi^ {(3)} + \dots , \tag {2.11}
$$

$$
\omega = \omega^ {(1)} + \omega^ {(3)} + \dots . \tag {2.12}
$$

These are based on the classical perturbation method, which assumes that some parameter (ε) naturally appearing in the governing equations, in our case the 

nonlinearity, is small. Typically, a perturbation analysis is performed in dimensionless variables, in which case $\varepsilon$ represents a given physical quantity such as the wavenumber times the wave amplitude $( k a )$ or the wave amplitude divided by the water depth $( a / h )$ . We prefer, however, to perform the analysis in dimensional variables and emphasize that in this case $\varepsilon$ has no physical meaning, but merely appears as a marker convenient for collecting terms of various orders of magnitude. Once the hierarchy of solutions have been obtained, $\varepsilon$ can be set to unity. 

# 3. The third-order formulation for irregular waves

# 3.1. The first-order solution

The starting point for the perturbation method is to choose a specific form of the first-order solution $\eta ^ { ( 1 ) }$ . Once this has been done, we can determine the equivalent form of the first-order velocity potential using the linear relationship $g \eta \simeq - \varPhi _ { t }$ at $z = 0$ . In the present work, we consider a first-order progressive wave group made up of $N$ frequencies $\omega _ { n }$ , where $n = 1 , 2 , \ldots , N$ . It should be emphasized that, to achieve a third-order formulation valid for irregular waves, we need to require that $N \geqslant 3$ , i.e. as a minimum the first-order wave train should consist of a trichromatic wave group. With the wavenumber vectors defined by $\pmb { k } _ { n } \equiv ( k _ { n x } , k _ { n y } )$ , and the phase functions given by 

$$
\theta_ {n} \equiv \omega_ {n} t - k _ {n x} x - k _ {n y} y, \tag {3.1}
$$

the first-order solution now reads 

$$
\eta^ {(1)} = \varepsilon \sum_ {n = 1} ^ {N} \left(a _ {n} \cos \theta_ {n} + b _ {n} \sin \theta_ {n}\right), \tag {3.2}
$$

$$
\boldsymbol {\Phi} ^ {(1)} = \boldsymbol {U} \cdot \boldsymbol {x} + \varepsilon \sum_ {n = 1} ^ {N} F _ {n} \cosh (\kappa_ {n} Z) (a _ {n} \sin \theta_ {n} - b _ {n} \cos \theta_ {n}), \tag {3.3}
$$

where $\boldsymbol { x } \equiv ( x , y )$ , $Z \equiv z + h$ and $U$ is an ambient current vector (constant in time and space). Note that the formulation has been given in terms of cosine and sine components as an alternative to cosine components with random phase shifts. Hence, (3.2) and (3.3) describe a linear multidirectional irregular wave train, and, provided that $N \geqslant 3$ , this can be lifted to third order as demonstrated in $\ S 3 . 3$ . 

By inserting (3.3) into the Laplace equation, we obtain 

$$
\kappa_ {n} = \left| \boldsymbol {k} _ {n} \right| = \sqrt {k _ {n x} ^ {2} + k _ {n y} ^ {2}}. \tag {3.4}
$$

Next, we insert (3.2) and (3.3) into (2.8) and (2.9), which again are inserted in (2.4) and (2.5). We collect terms proportional to $O ( \varepsilon )$ and obtain two independent homogeneous equations. These are satisfied by the linear dispersion relation (including a Doppler shift from the ambient current) 

$$
\omega_ {n} ^ {(1)} = \boldsymbol {k} _ {n} \cdot \boldsymbol {U} + \omega_ {1 n}, \quad \omega_ {1 n} \equiv \sqrt {g \kappa_ {n} \tanh  (h \kappa_ {n})}, \tag {3.5}
$$

and by the velocity coefficients 

$$
F _ {n} = \frac {- \omega_ {1 n}}{\kappa_ {n} \sinh (h \kappa_ {n})}. \tag {3.6}
$$

This defines the first-order solution in the perturbation method. 

# 3.2. The second-order solution

In the derivation of the second-order surface elevation, it is convenient to first preassess its form by calculating the square of the first-order solution divided by $h$ . As discussed by Madsen & Fuhrman (2006), this procedure makes sense because the governing equations contain quadratic nonlinearities, and it defines the form (but not the magnitude) of $\eta ^ { ( 2 ) }$ . Similarly, the form (but not the magnitude) of the second-order velocity potential can be pre-assessed by using $g \eta ^ { ( 2 ) } \simeq - \bar { \phi _ { t } ^ { ( 2 ) } }$ at $z = 0$ . This procedure leads to the following form of the second-order bound waves: 

$$
\begin{array}{l} \eta^ {(2)} = \varepsilon^ {2} \sum_ {n = 1} ^ {N} G _ {2 n} \left(A _ {2 n} \cos 2 \theta_ {n} + B _ {2 n} \sin 2 \theta_ {n}\right) \\ + \varepsilon^ {2} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} G _ {n \pm m} \left(A _ {n \pm m} \cos \theta_ {n \pm m} + B _ {n \pm m} \sin \theta_ {n \pm m}\right), \tag {3.7} \\ \end{array}
$$

$$
\begin{array}{l} \Phi^ {(2)} = \varepsilon^ {2} \sum_ {n = 1} ^ {N} F _ {2 n} \cosh (\kappa_ {2 n} Z) (A _ {2 n} \sin 2 \theta_ {n} - B _ {2 n} \cos 2 \theta_ {n}) \\ + \varepsilon^ {2} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} F _ {n \pm m} \cosh (\kappa_ {n \pm m} Z) (A _ {n \pm m} \sin \theta_ {n \pm m} - B _ {n \pm m} \cos \theta_ {n \pm m}). \tag {3.8} \\ \end{array}
$$

With $N = 3$ , the first summation in (3.7) and (3.8) contains three terms, which represent the self–self interactions $2 \theta _ { n }$ (known from Stokes theory for monochromatic waves). Then follows a double summation containing three sum-interaction terms (using the upper signs) and three difference-interaction terms (using the lower signs). These six terms represent pairwise bichromatic interactions with phase functions defined by 

$$
\theta_ {n \pm m} \equiv \theta_ {n} \pm \theta_ {m}. \tag {3.9}
$$

The amplitude coefficients for the bichromatic interactions are given by 

$$
A _ {n \pm m} \equiv \frac {1}{h} \left(a _ {n} a _ {m} \mp b _ {n} b _ {m}\right), \quad B _ {n \pm m} \equiv \frac {1}{h} \left(a _ {m} b _ {n} \pm a _ {n} b _ {m}\right), \tag {3.10}
$$

while the coefficients for self–self interaction read 

$$
A _ {2 n} \equiv \frac {1}{2} A _ {n + n} = \frac {1}{2 h} \left(a _ {n} ^ {2} - b _ {n} ^ {2}\right), \quad B _ {2 n} \equiv \frac {1}{2} B _ {n + n} = \frac {1}{2 h} a _ {n} b _ {n}. \tag {3.11}
$$

Next, we determine the wavenumbers $\kappa _ { 2 n }$ and $\kappa _ { n \pm m }$ by requiring that (3.8) should satisfy the Laplace equation, and this leads to 

$$
\kappa_ {n \pm m} = | \boldsymbol {k} _ {n} \pm \boldsymbol {k} _ {m} | = \sqrt {\left(k _ {n x} \pm k _ {m x}\right) ^ {2} + \left(k _ {n y} \pm k _ {m y}\right) ^ {2}}, \tag {3.12}
$$

$$
\kappa_ {2 n} = \left| 2 \boldsymbol {k} _ {n} \right| = 2 \kappa_ {n}. \tag {3.13}
$$

The corresponding frequency relations read 

$$
\omega_ {n \pm m} \equiv \omega_ {1 n} \pm \omega_ {1 m}, \quad \omega_ {2 n} \equiv 2 \omega_ {1 n}, \tag {3.14}
$$

where $\omega _ { 1 n }$ and $\omega _ { 1 m }$ satisfy (3.5). 

# 3.2.1. The second-order transfer functions

At this stage the solution specified in (3.7) and (3.8) automatically satisfies the Laplace equation and the kinematic bottom condition. The remaining problem is 

to satisfy the nonlinear surface conditions. First, we insert $\eta = \eta ^ { ( 1 ) } + \eta ^ { ( 2 ) }$ and $\phi = \phi ^ { ( 1 ) } + \phi ^ { ( 2 ) }$ into (2.8) and (2.9), which again are inserted in (2.4) and (2.5). Next, terms of order $O ( \varepsilon ^ { 2 } )$ are collected and we obtain a set of algebraic equations for the determination of the second-order transfer functions $G _ { 2 n }$ , $G _ { n \pm m }$ , $F _ { 2 n }$ and $F _ { n \pm m }$ . In order to achieve a compact formulation of the result, we introduce the definitions 

$$
\alpha_ {n \pm m} \equiv \omega_ {n \pm m} \cosh \left(h \kappa_ {n \pm m}\right), \tag {3.15}
$$

$$
\gamma_ {n \pm m} \equiv \kappa_ {n \pm m} \sinh (h \kappa_ {n \pm m}), \tag {3.16}
$$

$$
\beta_ {n \pm m} \equiv \omega_ {n \pm m} ^ {2} \cosh \left(h \kappa_ {n \pm m}\right) - g \kappa_ {n \pm m} \sinh \left(h \kappa_ {n \pm m}\right). \tag {3.17}
$$

The solution for $G _ { n + m }$ , which is the super-harmonic transfer function for the surface elevation, becomes 

$$
\begin{array}{l} G _ {n + m} = \Lambda_ {2} [ \bullet ] \equiv \frac {h}{2 \omega_ {1 n} \omega_ {1 m} \beta_ {n + m}} \left(g \alpha_ {n + m} \left(\omega_ {1 n} \left(\kappa_ {m} ^ {2} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right) + \omega_ {1 m} \left(\kappa_ {n} ^ {2} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right)\right) \right. \\ + \gamma_ {n + m} \left(g ^ {2} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \omega_ {1 n} ^ {2} \omega_ {1 m} ^ {2} - \omega_ {1 n} \omega_ {1 m} \omega_ {n + m} ^ {2}\right)). \tag {3.18} \\ \end{array}
$$

We note that it is a function of the following arguments: 

$$
\Lambda_ {2} [ \bullet ] \equiv \Lambda_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {1 m}, \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\omega_ {n + m}, \alpha_ {n + m}, \gamma_ {n + m}, \beta_ {n + m} \right\} \right]. \tag {3.19}
$$

The solution for $G _ { n - m }$ , which is the sub-harmonic transfer function for the surface −elevation, can be determined by switching the arguments in (3.19) using 

$$
G _ {n - m} = \Lambda_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{- \omega_ {1 m}, - \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\omega_ {n - m}, \alpha_ {n - m}, \gamma_ {n - m}, \beta_ {n - m} \right\} \right]. \tag {3.20}
$$

Similarly, the super-harmonic velocity potential transfer function $F _ { n + m }$ becomes 

$$
\begin{array}{l} F _ {n + m} = \Gamma_ {2} [ \bullet ] \equiv \frac {h}{2 \omega_ {1 n} \omega_ {1 m} \beta_ {n + m}} (\omega_ {1 n} \omega_ {1 m} \omega_ {n + m} (\omega_ {n + m} ^ {2} - \omega_ {1 n} \omega_ {1 m}) \\ - g ^ {2} \omega_ {1 n} \left(\kappa_ {m} ^ {2} + 2 \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right) - g ^ {2} \omega_ {1 m} \left(\kappa_ {n} ^ {2} + 2 \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right), \tag {3.21} \\ \end{array}
$$

where $\varGamma$ is a function of the following arguments: 

$$
\Gamma_ {2} [ \bullet ] \equiv \Gamma_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {1 m}, \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\omega_ {n + m}, \beta_ {n + m} \right\} \right]. \tag {3.22}
$$

The sub-harmonic function $F _ { n - m }$ can be determined by 

$$
F _ {n - m} = \Gamma_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{- \omega_ {1 m}, - \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\omega_ {n - m}, \beta_ {n - m} \right\} \right]. \tag {3.23}
$$

We emphasize that the expressions for $G _ { n \pm m }$ and $F _ { n \pm m }$ agree with the original derivation by Sharma & Dean (1981). 

Also the transfer functions for the self–self interactions can be determined from (3.18) and (3.21) using, for example, 

$$
G _ {2 n} \equiv G _ {n + n} = \Lambda_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {n + n}, \alpha_ {n + n}, \gamma_ {n + n}, \beta_ {n + n} \right\} \right], \tag {3.24}
$$

$$
F _ {2 n} \equiv F _ {n + n} = \Gamma_ {2} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {n + n}, \beta_ {n + n} \right\} \right], \tag {3.25}
$$

and these expressions simplify to the classical solution by Stokes (1847) for monochromatic waves: 

$$
G _ {2 n} = \frac {1}{2} h \kappa_ {n} \left(2 + \cosh 2 h \kappa_ {n}\right) \frac {\coth h \kappa_ {n}}{\sinh^ {2} h \kappa_ {n}}, \quad F _ {2 n} = - \frac {3}{4} \frac {h \omega_ {1 n}}{\sinh^ {4} h \kappa_ {n}}. \tag {3.26}
$$

In this connection we must emphasize that an unfortunate error appears in Madsen & Fuhrman (2006), where we claimed that $\begin{array} { r } { G _ { 2 n } \equiv \frac { 1 } { 2 } G _ { n + n } } \end{array}$ and $F _ { 2 n } \equiv \overset { \cdot } { \underset { 2 } { \cdot } } F _ { n + n }$ . The correct 

connection is given by (3.24) and (3.25). Finally, we note that $G _ { n \pm m }$ and $F _ { n \pm m }$ are even/odd functions, respectively, i.e. 

$$
G _ {m \pm n} = G _ {n \pm m}, \quad F _ {m \pm n} = \pm F _ {n \pm m}. \tag {3.27}
$$

# 3.3. The third-order solution

The derivation of the third-order solution follows the principles outlined in connection with the second-order solution. First, we pre-assess the form of the surface elevation by calculating the square of $\eta ^ { ( 1 ) } + \eta ^ { ( 2 ) }$ and dividing by $h$ , while collecting terms of order $\varepsilon ^ { 3 }$ . Next, we pre-assess the form of the third-order velocity potential by using $g \eta ^ { ( 3 ) } \simeq - \phi _ { t } ^ { ( 3 ) }$ at $z = 0$ . This procedure leads to the following form of the third-order bound waves: 

$$
\begin{array}{l} \eta^ {(3)} = \varepsilon^ {3} \sum_ {n = 1} ^ {N} G _ {1 3 n} \left(a _ {n} \cos \theta_ {n} + b _ {n} \sin \theta_ {n}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} G _ {3 n} \left(A _ {3 n} \cos 3 \theta_ {n} + B _ {3 n} \sin 3 \theta_ {n}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} G _ {n \pm 2 m} \left(A _ {n \pm 2 m} \cos \theta_ {n \pm 2 m} + B _ {n \pm 2 m} \sin \theta_ {n \pm 2 m}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} G _ {2 n \pm m} \left(A _ {2 n \pm m} \cos \theta_ {2 n \pm m} + B _ {2 n \pm m} \sin \theta_ {2 n \pm m}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \sum_ {p = m + 1} ^ {N} G _ {n \pm m \pm p} \left(A _ {n \pm m \pm p} \cos \theta_ {n \pm m \pm p} + B _ {n \pm m \pm p} \sin \theta_ {n \pm m \pm p}\right), \tag {3.28} \\ \end{array}
$$

$$
\begin{array}{l} \Phi^ {(3)} = \varepsilon^ {3} \sum_ {n = 1} ^ {N} F _ {1 3 n} \cosh (\kappa_ {n} Z) (a _ {n} \sin \theta_ {n} - b _ {n} \cos \theta_ {n}) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} F _ {3 n} \cosh (\kappa_ {3 n} Z) (A _ {3 n} \sin 3 \theta_ {n} - B _ {3 n} \cos 3 \theta_ {n}) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} F _ {n \pm 2 m} \cosh (\kappa_ {n \pm 2 m} Z) (A _ {n \pm 2 m} \sin \theta_ {n \pm 2 m} - B _ {n \pm 2 m} \cos \theta_ {n \pm 2 m}) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} F _ {2 n \pm m} \cosh (\kappa_ {2 n \pm m} Z) (A _ {2 n \pm m} \sin \theta_ {2 n \pm m} - B _ {2 n \pm m} \cos \theta_ {2 n \pm m}) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \sum_ {p = m + 1} ^ {N} F _ {n \pm m \pm p} \cosh (\kappa_ {n \pm m \pm p} Z) \\ \times \left(A _ {n \pm m \pm p} \sin \theta_ {n \pm m \pm p} - B _ {n \pm m \pm p} \cos \theta_ {n \pm m \pm p}\right). \tag {3.29} \\ \end{array}
$$

With $N = 3$ , the first summation in (3.28) and (3.29) contains three terms, representing the third-order correction to the first-order velocity potential. These terms, which have phase functions identical to the first-order solution, play a role in the necessary removal of secular terms. The next summation also contains three terms, which represent the self–self–self interactions $3 \theta _ { n }$ (known from Stokes theory for 

monochromatic waves). Then follows two double summations each containing three sum-interaction terms (using the upper signs) and three difference-interaction terms (using the lower signs). These 12 terms represent pairwise bichromatic interactions with phase functions defined by 

$$
\theta_ {n \pm 2 m} \equiv \theta_ {n} \pm 2 \theta_ {m}, \quad \theta_ {2 n \pm m} \equiv 2 \theta_ {n} \pm \theta_ {m}. \tag {3.30}
$$

Finally, there follows the triple summation containing four terms, which represent the trichromatic interactions with phase functions defined by 

$$
\theta_ {n \pm m \pm p} \equiv \theta_ {n} \pm \theta_ {m} \pm \theta_ {p}. \tag {3.31}
$$

Hence with $N = 3$ , the third-order contribution consists of 22 interaction terms. 

The pre-assessment also provides the amplitude coefficients for the trichromatic interactions in (3.28) and (3.29), which can be expressed by 

$$
A _ {n \pm m \pm p} = \Theta_ {A} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, \pm b _ {m} \right\}, \left\{a _ {p}, \pm b _ {p} \right\} \right], \tag {3.32}
$$

$$
B _ {n \pm m \pm p} = \Theta_ {B} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, \pm b _ {m} \right\}, \left\{a _ {p}, \pm b _ {p} \right\} \right], \tag {3.33}
$$

where 

$$
\Theta_ {A} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, b _ {m} \right\}, \left\{a _ {p}, b _ {p} \right\} \right] \equiv \frac {a _ {n} a _ {m} a _ {p} - b _ {n} b _ {m} a _ {p} - b _ {n} a _ {m} b _ {p} - a _ {n} b _ {m} b _ {p}}{h ^ {2}}, \tag {3.34}
$$

$$
\Theta_ {B} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, b _ {m} \right\}, \left\{a _ {p}, b _ {p} \right\} \right] \equiv \frac {b _ {n} a _ {m} a _ {p} + a _ {n} b _ {m} a _ {p} + a _ {n} a _ {m} b _ {p} - b _ {n} b _ {m} b _ {p}}{h ^ {2}}. \tag {3.35}
$$

It turns out that the pairwise bichromatic interactions and the self–self–self interactions can also be expressed in terms of the trichromatic coefficients and we obtain, for example, 

$$
A _ {n \pm 2 m} = \frac {1}{2} \Theta_ {A} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, \pm b _ {m} \right\}, \left\{a _ {m}, \pm b _ {m} \right\} \right], \tag {3.36}
$$

$$
B _ {n \pm 2 m} = \frac {1}{2} \Theta_ {B} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {m}, \pm b _ {m} \right\}, \left\{a _ {m}, \pm b _ {m} \right\} \right] \tag {3.37}
$$

and 

$$
A _ {3 n} = \frac {1}{2} \Theta_ {A} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {n}, b _ {n} \right\} \right], \tag {3.38}
$$

$$
B _ {3 n} = \frac {1}{2} \Theta_ {B} \left[ \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {n}, b _ {n} \right\}, \left\{a _ {n}, b _ {n} \right\} \right]. \tag {3.39}
$$

# 3.3.1. Wavenumber and frequency relations

The next step is to satisfy the Laplace equation by (3.29). This leads to the following wavenumber combinations: 

$$
\kappa_ {n \pm m \pm p} = | \boldsymbol {k} _ {n} \pm \boldsymbol {k} _ {m} \pm \boldsymbol {k} _ {p} | = \sqrt {\left(k _ {n x} \pm k _ {m x} \pm k _ {p x}\right) ^ {2} + \left(k _ {n y} \pm k _ {m y} \pm k _ {p y}\right) ^ {2}}, \tag {3.40}
$$

$$
\kappa_ {n \pm 2 m} = | \boldsymbol {k} _ {n} \pm 2 \boldsymbol {k} _ {m} | = \sqrt {\left(k _ {n x} \pm 2 k _ {m x}\right) ^ {2} + \left(k _ {n y} \pm 2 k _ {m y}\right) ^ {2}}, \tag {3.41}
$$

$$
\kappa_ {2 n \pm m} = | 2 \boldsymbol {k} _ {n} \pm \boldsymbol {k} _ {m} | = \sqrt {\left(2 k _ {n x} \pm k _ {m x}\right) ^ {2} + \left(2 k _ {n y} \pm k _ {m y}\right) ^ {2}}, \tag {3.42}
$$

$$
\kappa_ {3 n} = \left| 3 \boldsymbol {k} _ {n} \right| = 3 \kappa_ {n}. \tag {3.43}
$$

The corresponding frequency relations read 

$$
\left. \begin{array}{l} \omega_ {n \pm m \pm p} \equiv \omega_ {1 n} \pm \omega_ {1 m} \pm \omega_ {1 p}, \quad \omega_ {n \pm 2 m} \equiv \omega_ {1 n} \pm 2 \omega_ {1 m}, \\ \omega_ {2 n \pm m} \equiv 2 \omega_ {1 n} \pm \omega_ {1 m}, \quad \omega_ {3 n} \equiv 3 \omega_ {1 n}, \end{array} \right\} \tag {3.44}
$$

where $\omega _ { 1 n }$ , $\omega _ { 1 m }$ and $\omega _ { 1 p }$ satisfy (3.5). 

# 3.3.2. The third-order transfer functions

The remaining problem is to satisfy the nonlinear kinematic and dynamic surface conditions. In order to do so, we insert $\eta = \eta ^ { ( 1 ) } + \eta ^ { ( 2 ) } + \eta ^ { ( 3 ) }$ and $\phi = \phi ^ { ( 1 ) } + \phi ^ { ( 2 ) } + \phi ^ { ( 3 ) }$ into (2.8) and (2.9), which are again inserted in (2.4) and (2.5). In this process the phase functions should incorporate the frequencies 

$$
\omega_ {n} = \boldsymbol {k} _ {n} \cdot \boldsymbol {U} + \omega_ {1 n} \left(1 + \varepsilon^ {2} \omega_ {3 n}\right), \tag {3.45}
$$

where $\omega _ { 3 n }$ defines the third-order amplitude dispersion, which is necessary in order to remove secular terms (i.e. terms that resonate with the first-order solution). Terms of order $O ( \varepsilon ^ { 3 } )$ are collected, and we obtain a set of algebraic equations for the determination of the third-order transfer functions $G _ { n \pm m \pm p }$ , $G _ { n \pm 2 m }$ , $G _ { 2 n \pm m }$ , $G _ { 3 n }$ and $F _ { n \pm m \pm p }$ , $F _ { n \pm 2 m }$ , $F _ { 2 n \pm m }$ , $F _ { 3 n }$ . In order to achieve a compact formulation of the result, we introduce the following definitions: 

$$
\beta_ {n \pm m \pm p} \equiv \omega_ {n \pm m \pm p} ^ {2} \cosh \left(h \kappa_ {n \pm m \pm p}\right) - g \kappa_ {n \pm m \pm p} \sinh \left(h \kappa_ {n \pm m \pm p}\right), \tag {3.46}
$$

$$
\beta_ {n \pm 2 m} \equiv \omega_ {n \pm 2 m} ^ {2} \cosh (h \kappa_ {n \pm 2 m}) - g \kappa_ {n \pm 2 m} \sinh (h \kappa_ {n \pm 2 m}), \tag {3.47}
$$

$$
\beta_ {2 n \pm m} \equiv \omega_ {2 n \pm m} ^ {2} \cosh (h \kappa_ {2 n \pm m}) - g \kappa_ {2 n \pm m} \sinh (h \kappa_ {2 n \pm m}), \tag {3.48}
$$

$$
\beta_ {3 n} \equiv \omega_ {3 n} ^ {2} \cosh (h \kappa_ {3 n}) - g \kappa_ {3 n} \sinh (h \kappa_ {3 n}) \tag {3.49}
$$

and 

$$
\alpha_ {n \pm m \pm p} \equiv \omega_ {n \pm m \pm p} \cosh \left(h \kappa_ {n \pm m \pm p}\right), \quad \gamma_ {n \pm m \pm p} \equiv \kappa_ {n \pm m \pm p} \sinh \left(h \kappa_ {n \pm m \pm p}\right), \tag {3.50}
$$

$$
\alpha_ {n \pm 2 m} \equiv \omega_ {n \pm 2 m} \cosh (h \kappa_ {n \pm 2 m}), \quad \gamma_ {n \pm 2 m} \equiv \kappa_ {n \pm 2 m} \sinh (h \kappa_ {n \pm 2 m}), \tag {3.51}
$$

$$
\alpha_ {2 n \pm m} \equiv \omega_ {2 n \pm m} \cosh (h \kappa_ {2 n \pm m}), \quad \gamma_ {2 n \pm m} \equiv \kappa_ {2 n \pm m} \sinh (h \kappa_ {2 n \pm m}). \tag {3.52}
$$

The solution for $G _ { n + m + p }$ becomes 

$$
\begin{array}{l} G _ {n + m + p} = \Lambda_ {3} [ \bullet ] \equiv \frac {h ^ {2}}{4 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \left(\omega_ {1 n} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {n} ^ {2}\right) \right. \right. \\ + \omega_ {1 m} (\boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {n} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {m} ^ {2}) + \omega_ {1 p} (\boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {n} + \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {m} + \kappa_ {p} ^ {2})) \\ + \gamma_ {n + m + p} \left(\frac {g}{\omega_ {1 n}} \left(\omega_ {1 m} \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {n} + \omega_ {1 p} \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {n} - \omega_ {n + m + p} \kappa_ {n} ^ {2}\right) \right. \\ + \frac {g}{\omega_ {1 m}} \left(\omega_ {1 n} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \omega_ {1 p} \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {m} - \omega_ {n + m + p} \kappa_ {m} ^ {2}\right) \\ \left. \left. + \frac {g}{\omega_ {1 p}} \left(\omega_ {1 n} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \omega_ {1 m} \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} - \omega_ {n + m + p} \kappa_ {p} ^ {2}\right)\right)\right) \\ - \frac {h F _ {n + m}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \cosh (h \kappa_ {n + m}) \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {n + m} ^ {2}\right) \right. \\ \left. + \gamma_ {n + m + p} \left(\frac {g}{\omega_ {1 p}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p}\right) \cosh (h \kappa_ {n + m}) - \gamma_ {n + m} \omega_ {n + m + p}\right)\right) \\ - \frac {h F _ {n + p}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \cosh (h \kappa_ {n + p}) \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {n + p} ^ {2}\right) \right. \\ \left. + \gamma_ {n + m + p} \left(\frac {g}{\omega_ {1 m}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p}\right) \cosh \left(h \kappa_ {n + p}\right) - \gamma_ {n + p} \omega_ {n + m + p}\right)\right) \\ - \frac {h F _ {m + p}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \cosh (h \kappa_ {m + p}) \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {m + p} ^ {2}\right) \right. \\ \end{array}
$$

$$
\begin{array}{l} \left. + \gamma_ {n + m + p} \left(\frac {g}{\omega_ {1 n}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p}\right) \cosh (h \kappa_ {m + p}) - \gamma_ {m + p} \omega_ {n + m + p}\right)\right) \\ + \frac {h G _ {n + m}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \frac {g}{\omega_ {1 p}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {p} ^ {2}\right) - \gamma_ {n + m + p} \omega_ {1 p} ^ {2}\right) \\ + \frac {h G _ {n + p}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \frac {g}{\omega_ {1 m}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {m} ^ {2}\right) - \gamma_ {n + m + p} \omega_ {1 m} ^ {2}\right) \\ + \frac {h G _ {m + p}}{2 \beta_ {n + m + p}} \left(\alpha_ {n + m + p} \frac {g}{\omega_ {1 n}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {n} ^ {2}\right) - \gamma_ {n + m + p} \omega_ {1 n} ^ {2}\right). \tag {3.53} \\ \end{array}
$$

This is one of four trichromatic transfer functions for the surface elevation, and we note that it is a function of the following arguments: 

$$
\begin{array}{l} \Lambda_ {3} [ \bullet ] \equiv \Lambda_ {3} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\omega_ {1 m}, \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\omega_ {1 p}, \boldsymbol {k} _ {p}, \kappa_ {p} \right\}, \left\{\kappa_ {n + m}, \gamma_ {n + m}, G _ {n + m}, F _ {n + m} \right\}, \right. \\ \{\kappa_ {n + p}, \gamma_ {n + p}, G _ {n + p}, F _ {n + p} \}, \{\kappa_ {m + p}, \gamma_ {m + p}, G _ {m + p}, F _ {m + p} \}, \\ \left. \left\{\omega_ {n + m + p}, \alpha_ {n + m + p}, \gamma_ {n + m + p}, \beta_ {n + m + p} \right\} \right]. \tag {3.54} \\ \end{array}
$$

The solutions for $G _ { n - m + p }$ , $G _ { n + m - p }$ and $G _ { n - m - p }$ can all be determined by switching the arguments in (3.54) and we generally obtain 

$$
\begin{array}{l} G _ {n \pm m \pm p} \equiv \Lambda_ {3} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{\pm \omega_ {1 m}, \pm \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{\pm \omega_ {1 p}, \pm \boldsymbol {k} _ {p}, \kappa_ {p} \right\}, \right. \\ \{\kappa_ {n \pm m}, \gamma_ {n \pm m}, G _ {n \pm m}, F _ {n \pm m} \}, \{\kappa_ {n \pm p}, \gamma_ {n \pm p}, G _ {n \pm p}, F _ {n \pm p} \}, \\ \left. \left\{\kappa_ {\pm m \pm p}, \gamma_ {\pm m \pm p}, G _ {\pm m \pm p}, F _ {\pm m \pm p} \right\}, \left\{\omega_ {n \pm m \pm p}, \alpha_ {n \pm m \pm p}, \gamma_ {n \pm m \pm p}, \beta_ {n \pm m \pm p} \right\} \right]. \tag {3.55} \\ \end{array}
$$

Similarly, we determine $F _ { n + m + p }$ to be 

$$
\begin{array}{l} F _ {n + m + p} = \Gamma_ {3} [ \bullet ] \equiv \frac {- g h ^ {2}}{4 \beta_ {n + m + p}} \left(\omega_ {1 n} (\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {n} ^ {2}) \right. \\ + \omega_ {1 m} (\boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {n} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {m} ^ {2}) + \omega_ {1 p} (\boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {n} + \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {m} + \kappa_ {p} ^ {2}) \\ + \frac {\omega_ {n + m + p}}{\omega_ {1 n}} \left(\omega_ {1 m} \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {n} + \omega_ {1 p} \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {n} - \omega_ {n + m + p} \kappa_ {n} ^ {2}\right) \\ + \frac {\omega_ {n + m + p}}{\omega_ {1 m}} \left(\omega_ {1 n} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \omega_ {1 p} \boldsymbol {k} _ {p} \cdot \boldsymbol {k} _ {m} - \omega_ {n + m + p} \kappa_ {m} ^ {2}\right) \\ \left. + \frac {\omega_ {n + m + p}}{\omega_ {1 p}} \left(\omega_ {1 n} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \omega_ {1 m} \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} - \omega_ {n + m + p} \kappa_ {p} ^ {2}\right)\right) \\ + \frac {h F _ {n + m}}{2 \beta_ {n + m + p}} \left(g \cosh \left(h \kappa_ {n + m}\right) \left(\left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {n + m} ^ {2}\right) \right. \right. \\ \left. + \frac {\omega_ {n + m + p}}{\omega_ {1 p}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p}\right)\right) - \gamma_ {n + m} \omega_ {n + m + p} ^ {2}) \\ + \frac {h F _ {n + p}}{2 \beta_ {n + m + p}} \left(g \cosh (h \kappa_ {n + p}) \left(\left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {n + p} ^ {2}\right) \right. \right. \\ \left. + \frac {\omega_ {n + m + p}}{\omega_ {1 m}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p}\right)\right) - \gamma_ {n + p} \omega_ {n + m + p} ^ {2}) \\ + \frac {h F _ {m + p}}{2 \beta_ {n + m + p}} \left(g \cosh \left(h \kappa_ {m + p}\right) \left(\left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {m + p} ^ {2}\right) \right. \right. \\ \end{array}
$$

Third-order theory for multi-directional irregular waves 

$$
\begin{array}{l} \left. + \frac {\omega_ {n + m + p}}{\omega_ {1 n}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p}\right)\right) - \gamma_ {m + p} \omega_ {n + m + p} ^ {2}) \\ + \frac {h G _ {n + m}}{2 \beta_ {n + m + p}} \left(\omega_ {1 p} ^ {2} \omega_ {n + m + p} - \frac {g ^ {2}}{\omega_ {1 p}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {p} ^ {2}\right)\right) \\ + \frac {h G _ {n + p}}{2 \beta_ {n + m + p}} \left(\omega_ {1 m} ^ {2} \omega_ {n + m + p} - \frac {g ^ {2}}{\omega_ {1 m}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {m} \cdot \boldsymbol {k} _ {p} + \kappa_ {m} ^ {2}\right)\right) \\ + \frac {h G _ {m + p}}{2 \beta_ {n + m + p}} \left(\omega_ {1 n} ^ {2} \omega_ {n + m + p} - \frac {g ^ {2}}{\omega_ {1 n}} \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {p} + \kappa_ {n} ^ {2}\right)\right). \tag {3.56} \\ \end{array}
$$

We note that ${ \varGamma } _ { 3 }$ is a function of the same arguments as $\varLambda _ { 3 }$ and consequently we can determine the four trichromatic transfer functions by 

$$
\begin{array}{l} F _ {n \pm m \pm p} \equiv \Gamma_ {3} [ \{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \}, \{\pm \omega_ {1 m}, \pm \boldsymbol {k} _ {m}, \kappa_ {m} \}, \{\pm \omega_ {1 p}, \pm \boldsymbol {k} _ {p}, \kappa_ {p} \}, \\ \{\kappa_ {n \pm m}, \gamma_ {n \pm m}, G _ {n \pm m}, F _ {n \pm m} \}, \{\kappa_ {n \pm p}, \gamma_ {n \pm p}, G _ {n \pm p}, F _ {n \pm p} \}, \\ \left\{\kappa_ {\pm m \pm p}, \gamma_ {\pm m \pm p}, G _ {\pm m \pm p}, F _ {\pm m \pm p} \right\}, \left\{\omega_ {n \pm m \pm p}, \beta_ {n \pm m \pm p} \right\} ]. \tag {3.57} \\ \end{array}
$$

Next, we focus on the third-order bichromatic transfer functions such as $G _ { n \pm 2 m }$ and $G _ { 2 n \pm m }$ . They have previously been determined by Madsen & Fuhrman (2006), but it turns out that they can be expressed in terms of the trichromatic transfer functions using 

$$
G _ {n + 2 m} = G _ {n + m + m}, \quad F _ {n + 2 m} = F _ {n + m + m}, \tag {3.58}
$$

$$
G _ {n - 2 m} = G _ {n - m - m}, \quad F _ {n - 2 m} = F _ {n - m - m}, \tag {3.59}
$$

$$
G _ {2 n + m} = G _ {n + n + m}, \quad F _ {2 n + m} = F _ {n + n + m}, \tag {3.60}
$$

$$
G _ {2 n - m} = G _ {n + n - m}, \quad F _ {2 n - m} = F _ {n + n - m}. \tag {3.61}
$$

As an example we get 

$$
\begin{array}{l} G _ {n - m - m} = \Lambda_ {3} \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n}, \kappa_ {n} \right\}, \left\{- \omega_ {1 m}, - \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \left\{- \omega_ {1 m}, - \boldsymbol {k} _ {m}, \kappa_ {m} \right\}, \right. \\ \{\kappa_ {n - m}, \gamma_ {n - m}, G _ {n - m}, F _ {n - m} \}, \{\kappa_ {n - m}, \gamma_ {n - m}, G _ {n - m}, F _ {n - m} \}, \\ \left\{\kappa_ {2 m}, \gamma_ {2 m}, G _ {2 m}, - F _ {2 m} \right\}, \left\{\omega_ {n - 2 m}, \alpha_ {n - 2 m}, \gamma_ {n - 2 m}, \beta_ {n - 2 m} \right\} ], \tag {3.62} \\ \end{array}
$$

where we have used that $G _ { - m - m } = G _ { 2 m }$ and $F _ { - m - m } = - F _ { 2 m }$ according to (3.27). 

Finally, the self–self–self transfer functions can be determined by 

$$
G _ {3 n} = \frac {1}{3} G _ {n + n + n}, \quad F _ {3 n} = \frac {1}{3} F _ {n + n + n}, \tag {3.63}
$$

which simplify to the classical solution by Stokes (1847) for monochromatic waves, i.e. 

$$
G _ {3 n} = \frac {3}{1 2 8} \frac {h ^ {2} \kappa_ {n} ^ {2}}{\sinh^ {6} h \kappa_ {n}} (1 4 + 1 5 \cosh 2 h \kappa_ {n} + 6 \cosh 4 h \kappa_ {n} + \cosh 6 h \kappa_ {n}), \tag {3.64}
$$

$$
F _ {3 n} = \frac {1}{3 2} \frac {h ^ {2} \kappa_ {n} \omega_ {1 n}}{\sinh^ {7} h \kappa_ {n}} (- 1 1 + 2 \cosh 2 h \kappa_ {n}). \tag {3.65}
$$

# 3.3.3. Third-order correction to the first-order potential

In order to remove secular terms at third order, it is necessary (in addition to the frequency expansion) to include a third-order correction to either the first-order elevation or the first-order velocity potential. We choose to correct the velocity 

potential, i.e. $G _ { 1 3 n } \equiv 0$ , and obtain the result 

$$
F _ {1 3 n} = c _ {n} ^ {2} \Upsilon_ {n n} + \sum_ {\substack {m = 1 \\ m \neq n}} ^ {N} c _ {m} ^ {2} \Upsilon_ {n m}, \tag{3.66}
$$

where $c _ { n } \equiv \sqrt { a _ { n } ^ { 2 } + b _ { n } ^ { 2 } }$ , and where 

$$
\Upsilon_ {n n} = \omega_ {1 n} \kappa_ {n} \left(\frac {- 1 3 + 2 4 \cosh 2 h \kappa_ {n} + \cosh 4 h \kappa_ {n}}{6 4 \sinh^ {5} h \kappa_ {n}}\right), \tag {3.67}
$$

$$
\begin{array}{l} \Upsilon_ {n m} = \Upsilon [ \bullet ] _ {n \neq m} \equiv \frac {g}{4 \omega_ {1 n} \omega_ {1 m} \cosh h \kappa_ {n}} \left(\omega_ {1 m} \left(\kappa_ {n} ^ {2} - \kappa_ {m} ^ {2}\right) - \omega_ {1 n} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right) \\ + \frac {G _ {n + m} + G _ {n - m}}{4 h \omega_ {1 n} ^ {2} \omega_ {1 m} \cosh h \kappa_ {n}} \left(g ^ {2} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \omega_ {1 m} ^ {3} \omega_ {1 n}\right) \\ - \frac {1}{4 h \cosh h \kappa_ {n}} \left(F _ {n + m} \kappa_ {n + m} \sinh h \kappa_ {n + m} + F _ {n - m} \kappa_ {n - m} \sinh h \kappa_ {n - m}\right) \\ + \frac {g F _ {n + m} \cosh h \kappa_ {n + m}}{4 h \omega_ {1 n} ^ {2} \omega_ {1 m} \cosh h \kappa_ {n}} \left(\left(\omega_ {1 n} + \omega_ {1 m}\right) \left(\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \kappa_ {m} ^ {2}\right) - \omega_ {1 m} \kappa_ {n + m} ^ {2}\right) \\ + \frac {g F _ {n - m} \cosh h \kappa_ {n - m}}{4 h \omega_ {1 n} ^ {2} \omega_ {1 m} \cosh h \kappa_ {n}} \left((\omega_ {1 n} - \omega_ {1 m}) (\boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} - \kappa_ {m} ^ {2}) - \omega_ {1 m} \kappa_ {n - m} ^ {2}\right). \tag {3.68} \\ \end{array}
$$

3.3.4. The time-averaged volume flux and the wave-induced return current The time-averaged volume flux vector (Eulerian drift) is defined by 

$$
\boldsymbol {M} \equiv \overline {{\int_ {- h} ^ {0} \nabla \Phi \mathrm {d} z}} + \overline {{\int_ {0} ^ {\eta} \nabla \Phi \mathrm {d} z}}, \tag {3.69}
$$

where $\mathbf { V }$ is the horizontal gradient operator, while the overbar represents the timeaveraging process. As discussed in Madsen & Fuhrman (2006), the second integral in (3.69) is evaluated by using Taylor series expansions from $z = 0$ ; and by substituting the third-order expressions for $\phi$ and $\eta$ , we obtain the result 

$$
\boldsymbol {M} = h \boldsymbol {U} + \varepsilon^ {2} \sum_ {n = 1} ^ {N} \left(\frac {c _ {n} ^ {2} \omega_ {1 n}}{2 \kappa_ {n}} \coth h \kappa_ {n}\right) \boldsymbol {k} _ {n} + O \left(\varepsilon^ {4}\right), \tag {3.70}
$$

which is a straightforward generalization of the Stokes second-order result for monochromatic unidirectional waves. 

Under certain conditions, e.g. in closed wave tanks, $M$ must be zero, and as a result a wave-induced return current will appear. According to (3.70) this return current can be determined by 

$$
\boldsymbol {U} = - \varepsilon^ {2} \sum_ {n = 1} ^ {N} \left(\frac {c _ {n} ^ {2} \omega_ {1 n}}{2 h \kappa_ {n}} \coth h \kappa_ {n}\right) \boldsymbol {k} _ {n}. \tag {3.71}
$$

# 3.3.5. The third-order dispersion relation

Finally, we need to determine $\omega _ { 3 n }$ , which represents the third-order correction to the dispersion relation 

$$
\omega_ {n} = \omega_ {1 n} \left(1 + \varepsilon^ {2} \omega_ {3 n}\right) + \boldsymbol {k} _ {n} \cdot \boldsymbol {U}. \tag {3.72}
$$

It turns out that the dispersion relation of each wave component will be influenced by the pairwise bichromatic interaction with all other components, and we obtain 

$$
\omega_ {3 n} = c _ {n} ^ {2} \kappa_ {n} ^ {2} \Omega_ {n n} + \sum_ {\substack {m = 1 \\ m \neq n}} ^ {N} c _ {m} ^ {2} \kappa_ {m} ^ {2} \Omega_ {n m}. \tag{3.73}
$$

The self–self interaction term (known from Stokes third-order theory) reads 

$$
\Omega_ {n n} = \frac {8 + \cosh \left(4 h \kappa_ {n}\right)}{1 6 \sinh^ {4} \left(h \kappa_ {n}\right)}, \tag {3.74}
$$

while the mutual interaction terms are given by 

$$
\begin{array}{l} \Omega_ {n m} = \Omega [ \bullet ] _ {n \neq m} = \frac {1}{\kappa_ {m} ^ {2}} \left(\frac {\left(2 \omega_ {1 m} ^ {2} + \omega_ {1 n} ^ {2}\right)}{4 \omega_ {1 n} \omega_ {1 m}} \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m} + \frac {1}{4} \kappa_ {m} ^ {2}\right) \\ + \left(G _ {n + m} + G _ {n - m}\right) \frac {1}{\kappa_ {m} ^ {2}} \left(\frac {g \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}}{4 h \omega_ {1 n} \omega_ {1 m}} - \frac {\omega_ {1 m} ^ {2}}{4 g h}\right) \\ + \frac {\omega_ {1 n}}{4 g h} \frac {1}{\kappa_ {m} ^ {2}} \left(F _ {n + m} \kappa_ {n + m} \sinh h \kappa_ {n + m} + F _ {n - m} \kappa_ {n - m} \sinh h \kappa_ {n - m}\right) \\ - \frac {F _ {n + m} \cosh h \kappa_ {n + m}}{4 h \omega_ {1 n} \omega_ {1 m}} \frac {1}{\kappa_ {m} ^ {2}} \left(\left(\omega_ {1 n} - \omega_ {1 m}\right) \left(\kappa_ {m} ^ {2} + \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right) + \omega_ {1 m} \kappa_ {n + m} ^ {2}\right) \\ + \frac {F _ {n - m} \cosh h \kappa_ {n - m}}{4 h \omega_ {1 n} \omega_ {1 m}} \frac {1}{\kappa_ {m} ^ {2}} \left(\left(\omega_ {1 n} + \omega_ {1 m}\right) \left(\kappa_ {m} ^ {2} - \boldsymbol {k} _ {n} \cdot \boldsymbol {k} _ {m}\right) - \omega_ {1 m} \kappa_ {n - m} ^ {2}\right). \tag {3.75} \\ \end{array}
$$

The nonlinear dispersion relation for pairwise bichromatic interactions in finite depth was analysed and discussed extensively by Madsen & Fuhrman (2006). Furthermore, we used this theory to explain numerical results for nonlinear Bragg scattering (see Madsen, Fuhrman & Wang 2006). In the present paper we have nothing further to add on this particular issue, although we emphasize that the amplitude dispersion due to mutual interactions between different wave components is a very important and useful result of the present theory. 

# 3.4. The velocity potential at the free surface

For numerical modelling of nonlinear waves, it is typically necessary to specify the velocity potential evaluated directly on the free surface. Hence we provide expressions for this quantity valid for multi-directional irregular waves to third order. The potential at the free surface is given by 

$$
\begin{array}{l} \widetilde {\boldsymbol {\Phi}} = \boldsymbol {U} \cdot \boldsymbol {x} + \varepsilon \sum_ {n = 1} ^ {N} \left(\mu_ {n} + \mu_ {n} ^ {*}\right) \left(a _ {n} \sin \theta_ {n} - b _ {n} \cos \theta_ {n}\right) \\ + \varepsilon^ {2} \sum_ {n = 1} ^ {N} \mu_ {2 n} \left(A _ {2 n} \sin 2 \theta_ {n} - B _ {2 n} \cos 2 \theta_ {n}\right) \\ + \varepsilon^ {2} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \mu_ {n \pm m} \left(A _ {n \pm m} \sin \theta_ {n \pm m} - B _ {n \pm m} \cos \theta_ {n \pm m}\right) \\ \end{array}
$$

$$
\begin{array}{l} + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \mu_ {3 n} \left(A _ {3 n} \sin 3 \theta_ {n} - B _ {3 n} \cos 3 \theta_ {n}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \mu_ {n \pm 2 m} \left(A _ {n \pm 2 m} \sin \theta_ {n \pm 2 m} - B _ {n \pm 2 m} \cos \theta_ {n \pm 2 m}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \mu_ {2 n \pm m} \left(A _ {2 n \pm m} \sin \theta_ {2 n \pm m} - B _ {2 n \pm m} \cos \theta_ {2 n \pm m}\right) \\ + \varepsilon^ {3} \sum_ {n = 1} ^ {N} \sum_ {m = n + 1} ^ {N} \sum_ {p = m + 1} ^ {N} \mu_ {n \pm m \pm p} \left(A _ {n \pm m \pm p} \sin \theta_ {n \pm m \pm p} - B _ {n \pm m \pm p} \cos \theta_ {n \pm m \pm p}\right), \tag {3.76} \\ \end{array}
$$

where 

$$
\mu_ {n} = F _ {n} \cosh h \kappa_ {n}, \tag {3.77}
$$

$$
\mu_ {n \pm m} = F _ {n \pm m} \cosh h \kappa_ {n \pm m} - \frac {1}{2} h \left(\omega_ {1 n} \pm \omega_ {1 m}\right), \tag {3.78}
$$

$$
\mu_ {2 n} = F _ {2 n} \cosh h \kappa_ {2 n} - h \omega_ {1 n}, \tag {3.79}
$$

$$
\mu_ {3 n} = F _ {3 n} \cosh h \kappa_ {3 n} - \frac {g h ^ {2}}{4} \frac {\kappa_ {n} ^ {2}}{\omega_ {1 n}} + \frac {h}{2} \left(F _ {2 n} \gamma_ {2 n} - \omega_ {1 n} G _ {2 n}\right), \tag {3.80}
$$

$$
\begin{array}{l} \mu_ {n + m + p} = \Pi [ \bullet ] = F _ {n + m + p} \cosh h \kappa_ {n + m + p} - \frac {g h ^ {2}}{4} \left(\frac {\kappa_ {n} ^ {2}}{\omega_ {1 n}} + \frac {\kappa_ {m} ^ {2}}{\omega_ {1 m}} + \frac {\kappa_ {p} ^ {2}}{\omega_ {1 p}}\right) \\ - \frac {h}{2} \left(\omega_ {1 n} G _ {m + p} + \omega_ {1 m} G _ {n + p} + \omega_ {1 p} G _ {n + m}\right) \\ + \frac {h}{2} \left(F _ {n + m} \gamma_ {n + m} + F _ {n + p} \gamma_ {n + p} + F _ {m + p} \gamma_ {m + p}\right). \tag {3.81} \\ \end{array}
$$

We note that $\Pi [ \bullet ]$ has the arguments 

$$
\begin{array}{l} \Pi [ \bullet ] \equiv \Pi [ \{\omega_ {1 n}, \boldsymbol {k} _ {n} \}, \{\omega_ {1 m}, \boldsymbol {k} _ {m} \}, \{\omega_ {1 p}, \boldsymbol {k} _ {p} \}, \{G _ {n + m}, F _ {n + m} \}, \\ \left\{G _ {n + p}, F _ {n + p} \right\}, \left\{G _ {m + p}, F _ {m + p} \right\}, \left\{F _ {n + m + p} \right\} ]. \tag {3.82} \\ \end{array}
$$

The solutions for $\mu _ { n - m + p }$ , $\mu _ { n + m - p }$ and $\mu _ { n - m - p }$ are determined by switching the corresponding arguments in $\Pi [ \bullet ]$ . This is also true for the third-order bichromatic interactions, which are determined by, for example, 

$$
\begin{array}{l} \mu_ {n \pm 2 m} = \Pi \left[ \left\{\omega_ {1 n}, \boldsymbol {k} _ {n} \right\}, \left\{\pm \omega_ {1 m}, \pm \boldsymbol {k} _ {m} \right\}, \left\{\pm \omega_ {1 m}, \pm \boldsymbol {k} _ {m} \right\}, \left\{G _ {n \pm m}, F _ {n \pm m} \right\}, \right. \\ \left\{G _ {n \pm p}, F _ {n \pm p} \right\}, \left\{G _ {2 m}, \pm F _ {2 m} \right\}, \left\{F _ {n \pm 2 m} \right\} ], \tag {3.83} \\ \end{array}
$$

and similar expressions for $\mu _ { 2 n \pm m }$ 

Finally, the third-order correction to the first-order potential leads to 

$$
\mu_ {n} ^ {*} = F _ {n} ^ {*} \cosh h \kappa_ {n} + c _ {n} ^ {2} \Xi_ {n n} + c _ {m} ^ {2} \Xi_ {n m} + c _ {p} ^ {2} \Xi_ {n p}, \tag {3.84}
$$

where 

$$
\Xi_ {n n} = \frac {1}{4 h} \left(\omega_ {1 n} G _ {2 n} + F _ {2 n} \kappa_ {2 n} \sinh h \kappa_ {2 n} - \frac {g h \kappa_ {n} ^ {2}}{2 \omega_ {1 n}}\right) \tag {3.85}
$$

and 

$$
\Xi_ {n m} = \Xi [ \bullet ] _ {n \neq m} = \frac {1}{2 h} \left(\omega_ {1 m} \left(G _ {n + m} - G _ {n - m}\right) + F _ {n + m} \gamma_ {n + m} + F _ {n - m} \gamma_ {n - m} - \frac {g h \kappa_ {n} ^ {2}}{2 \omega_ {1 n}}\right). \tag {3.86}
$$

This completes the third-order theory, which provides explicit expressions for the surface elevation, the amplitude dispersion and the vertical variation of the velocity potential for irregular multi-directional water waves in finite depth (expressed in terms of the mean water depth h). Although the theory involves a large number of terms and lengthy expressions, we have attempted to express it elegantly in the form of multi-purpose functional expressions. 

# 3.5. A simple example involving trichromatic interactions

We now illustrate the new theory on a third-order trichromatic interaction defined by the wavenumber vectors 

$$
\boldsymbol {k} _ {n} = \kappa_ {n} (\cos \varphi_ {n}, \sin \varphi_ {n}), \quad \boldsymbol {k} _ {m} = \kappa_ {m} (\cos \varphi_ {m}, \sin \varphi_ {m}), \quad \boldsymbol {k} _ {p} = \kappa_ {p} (\cos \varphi_ {p}, \sin \varphi_ {p}). \tag {3.87}
$$

We consider the following specifications: $h = 1 . 0 ~ \mathrm { m }$ , $g = 9 . 8 1 \mathrm { ~ m ~ s } ^ { - 2 }$ , $\omega _ { n } = 4 . 0 ~ \mathrm { s ^ { - 1 } }$ , $a _ { n } = 0 . 0 2 5 ~ \mathrm { m }$ , $\varphi _ { n } = 3 0 ^ { \circ }$ , $\omega _ { m } ^ { - } = 4 . 2 ~ \mathrm { s } ^ { - 1 }$ , $a _ { m } = 0 . 0 2 5 ~ \mathrm { m }$ , $\varphi _ { m } = - 2 0 ^ { \circ }$ , $\omega _ { p } = 4 . 4 ~ \mathrm { s ^ { - 1 } }$ , $a _ { p } = 0 . 0 5 0 ~ \mathrm { m }$ and $\varphi _ { p } = 0 ^ { \circ }$ . Additionally, we assume that the Eulerian drift must be zero (wave tank conditions) and consequently the wave-induced return current vector is determined by (3.71). Now the third-order dispersion relation leads to the wavenumbers $h \kappa _ { n } = 1 . 7 1 6 6 6$ , $h \kappa _ { m } = 1 . 8 5 7 3 7$ and $h \kappa _ { p } = 2 . 0 2 0 1 0$ . Consequently, the wave steepnesses are $a _ { n } \kappa _ { n } = 0 . 0 4 3$ , $a _ { m } \kappa _ { m } = 0 . 0 4 7$ and $a _ { p } \kappa _ { p } = 0 . 1 0 1$ . 

Figure $1 ( a )$ shows a perspective plot of the surface elevation, figure $1 ( b )$ shows the first- and third-order surface elevations along the centreline (i.e. $y = 0$ ), while figure $1 ( c )$ shows the first- and third-order velocity profiles at the point $( x , y ) = ( 0 , 0 )$ . The relevant coefficients corresponding to figure 1 are given to aid in checking any implementation of the theory (see tables 2 and 3 in the Appendix). 

# 4. Harmonic resonance and singularities at third order

Harmonic resonance may occur for certain combinations of frequencies and wavenumber vectors, and this heralds the breakdown of the perturbation theory owing to the inherent singularities in the transfer functions. In this section, we identify and discuss these combinations. 

As seen from $\ S 3$ , the third-order transfer functions for a trichromatic interaction generally have denominators 

$$
\beta_ {n \pm m \pm p} \equiv \omega_ {n \pm m \pm p} ^ {2} \cosh \left(h \kappa_ {n \pm m \pm p}\right) - g \kappa_ {n \pm m \pm p} \sinh \left(h \kappa_ {n \pm m \pm p}\right), \tag {4.1}
$$

where 

$$
\kappa_ {n \pm m \pm p} \equiv | \boldsymbol {k} _ {n} \pm \boldsymbol {k} _ {m} \pm \boldsymbol {k} _ {p} |, \quad \omega_ {n \pm m \pm p} \equiv \omega_ {1 n} \pm \omega_ {1 m} \pm \omega_ {1 p}, \tag {4.2}
$$

and where $\omega _ { 1 n }$ , $\omega _ { 1 m }$ and $\omega _ { 1 p }$ satisfy the linear dispersion relation (3.5). While $\omega _ { n \pm m \pm p }$ represents the frequency of the bound wave with wavenumber $\kappa _ { n \pm m \pm p }$ , we may introduce the frequency of the corresponding free wave as 

$$
\widetilde {\omega} _ {n \pm m \pm p} \equiv \sqrt {g \kappa_ {n \pm m \pm p} \tanh  h \kappa_ {n \pm m \pm p}}. \tag {4.3}
$$

Now the denominators (4.1) can be expressed as 

$$
\beta_ {n \pm m \pm p} \equiv \cosh \left(h \kappa_ {n \pm m \pm p}\right) \left(\omega_ {n \pm m \pm p} ^ {2} - \widetilde {\omega} _ {n \pm m \pm p} ^ {2}\right). \tag {4.4}
$$


(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/2f4ea453bd34c2aa831d0264353d2314316c972b96581e9d2bdd472d6972a489.jpg)



(b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/47a5002888e446d557db1aedc3064b526af569bcdf667b4a35030523b86a0791.jpg)



(c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/aefce56a2f9d501eab6ae26dc343b04aa3e5747ac900d1161a5f552d43615854.jpg)



FIGURE 1. (Colour online available at journals.cambridge.org/flm) Third-order solution for a trichromatic interaction in finite depth (specifications defined in $\ S 3 . 5 )$ ). (a) Perspective plot of the surface elevation. $( b )$ Surface elevation along the centreline $( y = 0$ ). (c) Velocity profile at the centre point $( x , y ) = ( 0 , 0 )$ . First-order theory (dashed line); third-order theory (full line).


For certain combinations of wave modes, these denominators may go to zero, and according to (4.2) and (4.4) this happens if 

$$
\widetilde {\omega} _ {n \pm m \pm p} \rightarrow \pm \omega_ {n \pm m \pm p} = \pm \left(\omega_ {1 n} \pm \omega_ {1 m} \pm \omega_ {1 p}\right), \quad \boldsymbol {k} _ {n \pm m \pm p} = \boldsymbol {k} _ {n} \pm \boldsymbol {k} _ {m} \pm \boldsymbol {k} _ {p}. \tag {4.5}
$$

This is satisfied when the frequency of the bound third-order component becomes identical to the frequency of a free surface wave with the same wavenumber. The condition (4.5) agrees with the general resonance condition first derived by Phillips (1960). 

# 4.1. Linear resonance curves for a short-crested wave interacting with a plane wave

In the following we shall focus on the canonical case of a monochromatic shortcrested wave interacting with a plane wave of another frequency, representing perhaps the simplest genuine three-wave configuration. We express this situation by defining the wavenumber vectors 

$$
\boldsymbol {k} _ {n} = \kappa (p, q), \quad \boldsymbol {k} _ {m} = \kappa (p, - q), \quad \boldsymbol {k} _ {p} = \kappa (1, 0), \tag {4.6}
$$

where 

$$
p = \rho \sin \varphi , \quad q = \rho \cos \varphi . \tag {4.7}
$$

In deep water, the corresponding frequencies satisfying the linear dispersion relation read 

$$
\omega_ {1 n} = \sqrt {g \rho \kappa}, \quad \omega_ {1 m} = \sqrt {g \rho \kappa}, \quad \omega_ {1 p} = \sqrt {g \kappa}. \tag {4.8}
$$

The advantage of this relatively simple case is that it can easily be studied in a physical or numerical wave tank with fully reflective lateral boundaries as long as the width of the tank $( w )$ is adjusted according to the wavenumber components in the transverse $y$ direction, e.g. $\kappa q w = \pi$ . It turns out that singularities occur for certain combinations of $p$ and $q$ , and these solutions will be pursued in the following. 

First of all, harmonic resonance can occur for the following third-order interactions: 

$$
3 a: \quad \boldsymbol {k} _ {3 a} = 2 \boldsymbol {k} _ {n} - \boldsymbol {k} _ {p}, \quad \kappa_ {3 a} = | \boldsymbol {k} _ {3 a} |, \quad \omega_ {3 a} = 2 \omega_ {1 n} - \omega_ {1 p}, \quad \widetilde {\omega} _ {3 a} = \sqrt {g \kappa_ {3 a}}, \tag {4.9}
$$

$$
3 b: \quad \boldsymbol {k} _ {3 b} = 2 \boldsymbol {k} _ {p} - \boldsymbol {k} _ {n}, \quad \kappa_ {3 b} = | \boldsymbol {k} _ {3 b} |, \quad \omega_ {3 b} = 2 \omega_ {1 p} - \omega_ {1 n}, \quad \widetilde {\omega} _ {3 b} = \sqrt {g \kappa_ {3 b}}, \tag {4.10}
$$

$$
3 c: \quad \boldsymbol {k} _ {3 c} = \boldsymbol {k} _ {n} + \boldsymbol {k} _ {m} - \boldsymbol {k} _ {p}, \quad \kappa_ {3 c} = | \boldsymbol {k} _ {3 c} |, \quad \omega_ {3 c} = \omega_ {1 n} + \omega_ {1 m} - \omega_ {1 p}, \quad \widetilde {\omega} _ {3 c} = \sqrt {g \kappa_ {3 c}}. \tag {4.11}
$$

The corresponding resonance conditions read: 

$$
3 a: \quad r _ {3 a} \equiv - 1 + 4 (p ^ {2} + q ^ {2}) ^ {1 / 4} - 4 (p ^ {2} + q ^ {2}) ^ {1 / 2} + ((1 - 2 p) ^ {2} + 4 q ^ {2}) ^ {1 / 2} = 0, \tag {4.12}
$$

$$
3 b: \quad r _ {3 b} \equiv - 4 + 4 \left(p ^ {2} + q ^ {2}\right) ^ {1 / 4} - \left(p ^ {2} + q ^ {2}\right) ^ {1 / 2} + \left((- 2 + p) ^ {2} + q ^ {2}\right) ^ {1 / 2} = 0, \tag {4.13}
$$

$$
3 c: \quad r _ {3 c} \equiv \left((- 1 + 2 p) ^ {2}\right) ^ {1 / 2} + \left(1 - 2 \left(p ^ {2} + q ^ {2}\right) ^ {1 / 4}\right) ^ {2} = 0. \tag {4.14}
$$

The solutions of (4.12)–(4.14) are shown in figure 2 in terms of $( p , q )$ . These are the $x$ and $y$ components of the $\pmb { k } _ { n }$ vector, which starts at the origin and ends on one of the three curves 3a, 3b or $3 c$ . By definition $k _ { m }$ is always the mirror of $\pmb { k } _ { n }$ , while $k _ { p }$ is represented by the unit vector $( 1 , 0 )$ in figure 2. We note that the curve for case $_ { 3 b }$ is identical to the ‘figure-eight’ as shown by Phillips (1960) and Longuet-Higgins (1962). 

An alternative way of presenting the deep-water result is to insert (4.7) in (4.12)–(4.14), which leads to the resonance solutions 

$$
3 a: \quad \sin \varphi = - 6 - 3 \rho + 8 \rho^ {1 / 2} + 2 \rho^ {- 1 / 2}, \tag {4.15}
$$

$$
3 b: \quad \sin \varphi = - 6 - 3 \rho^ {- 1} + 8 \rho^ {- 1 / 2} + 2 \rho^ {1 / 2}, \tag {4.16}
$$

$$
3 c: \quad \sin \varphi = - 2 + 2 \rho^ {- 1 / 2}. \tag {4.17}
$$

These solutions are shown in figure 3 in terms of $\rho$ and $\varphi$ . 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/967cedf75f53a0eaab0b093a40d7e1d3db4f4c371843b28bdb847bfe96710296.jpg)



FIGURE 2. (Colour online) Harmonic resonance curves in deep water for the case of a monochromatic short-crested wave interacting with a plane wave of a different frequency. Wavenumber vectors: $k _ { n } = \kappa ( p , q ) , k _ { m } = \kappa ( p , \bar { - } q )$ $\pmb { k } _ { n } = \kappa ( p , q )$ and $\pmb { k } _ { p } = \kappa ( 1 , 0 )$ . Curves in terms of $( p , q )$ for $\pmb { k } _ { 3 a } = 2 \pmb { k } _ { n } - \pmb { k } _ { p }$ , $\pmb { k } _ { 3 b } = 2 \pmb { k } _ { p } - \pmb { k } _ { n }$ and $k _ { 3 c } = k _ { n } + k _ { m } - k _ { p }$ .


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/074d186b0f66cf434238bef6809d6be6bde2e7e82ee36b609d8f4c259affce73.jpg)



FIGURE 3. (Colour online) Specifications as in figure 2, but solution shown in terms of $\rho$ and $\varphi$ , where $p = \rho \sin \varphi$ and $q = \rho \cos \varphi$ .


Finally, the case of intermediate water depth $( \kappa h = 1 . 0 )$ ) is treated in figure 4. Now (4.8) and the last column of (4.9)–(4.11) are replaced by the general linear dispersion relation, and the resonance condition solved numerically. We notice that the resonance curves are similar to the deep-water case shown in figure 2, except for a general shrinking of the domain. 

# 4.2. Wave trains exposed to harmonic resonance

The third-order perturbation theory presented in this paper will obviously fail whenever we are in the vicinity of harmonic resonance. This begs the question: Is 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/e32bb70f6fd3c9af75612f5762a9a9d98e520f3d8f2e4b05dfff00755c086b4f.jpg)



FIGURE 4. (Colour online) Harmonic resonance curves in finite depth for the case of a monochromatic short-crested wave interacting with a plane wave of a different frequency. Wavenumber vectors: $\pmb { k } _ { n } = \kappa ( p , q )$ , $\pmb { k } _ { m } = \kappa ( p , - q )$ and $\bar { \pmb { k } _ { p } } = \kappa ( 1 , 0 )$ . Curves in terms of $( p , q )$ for $\pmb { k } _ { 3 a } = 2 \pmb { k } _ { n } - \pmb { k } _ { p }$ , $\pmb { k } _ { 3 b } = 2 \pmb { k } _ { p } - \pmb { k } _ { n }$ and $k _ { 3 c } = k _ { n } + k _ { m } - \dot { k } _ { p }$ . Finite depth: $\kappa h = 1 . 0$ .


this a mathematical failure of the theory or is it physically evident? Moreover, what will actually happen to a wave train exposed to harmonic resonance? In order to study these issues, we make numerical simulations of waves propagating up to 1000 water depths (in the $x$ direction) over a flat bottom. As mentioned earlier, the width of the numerical wave tank $( w )$ is adjusted according to the wavenumber components in the transverse y direction, e.g. $\kappa q w = \pi$ , so that the lateral boundaries can be fully reflective. Again, the test cases are well suited for investigation in laboratory facilities, although it does require a relatively long wave flume. 

In the following we shall test the response to various second- or third-order boundary conditions. If these conditions are adequate, the input waves will travel down the flume as steady waves with constant harmonic amplitudes. Although the simulations will be conducted in the time domain and the natural output is the time and space variation of, for example, the surface elevation, we shall focus on the spatial evolution of the relevant harmonics for a range of wave conditions, obtained through regression techniques. The vehicle for this investigation is a numerical model, which solves the high-order Boussinesq-type formulation by Madsen, Bingham & Liu (2002) and Madsen, Bingham & Schaffer (¨ 2003). Linear and nonlinear properties are accurately represented up to wavenumber times water depth $k h \simeq 2 5$ , while the interior velocity field is accurate up to $k h \simeq 1 2$ even for highly nonlinear waves. A detailed description of the numerical scheme can be found in Madsen et al. (2002) for one horizontal dimension, and in Fuhrman & Bingham (2004) for two horizontal dimensions. 

For the interaction between a monochromatic short-crested wave and a plane wave at another frequency (as defined by (4.6) and (4.7)), we can expect energy on two first-order frequencies, four second-order frequencies and six third-order frequencies. We introduce the notation $a _ { i j }$ and $\omega _ { i j }$ for the relevant amplitudes and frequencies, where $j$ covers the interval from unity to the maximum number of frequencies at the 

order i. The relevant amplitude and frequency functions are given by 

$$
\left. \begin{array}{l} a _ {1 1} \equiv \left| a _ {n} + a _ {m} + G _ {2 n - m} A _ {2 n - m} + G _ {n - 2 m} A _ {n - 2 m} \right| \quad \text {a t} \omega_ {1 1} \equiv \omega_ {n}, \\ a _ {1 2} \equiv \left| a _ {p} + G _ {n - m + p} A _ {n - m + p} + G _ {n - m - p} A _ {n - m - p} \right| \quad \text {a t} \omega_ {1 2} \equiv \omega_ {p}, \\ a _ {2 1} \equiv \left| G _ {2 n} A _ {2 n} + G _ {2 m} A _ {2 m} + G _ {n + m} A _ {n + m} \right| \quad \text {a t} \omega_ {2 1} \equiv 2 \omega_ {n}, \\ a _ {2 2} \equiv \left| G _ {n + p} A _ {n + p} + G _ {m + p} A _ {m + p} \right| \quad \text {a t} \omega_ {2 2} \equiv \omega_ {n} + \omega_ {p}, \\ a _ {2 3} \equiv \left| G _ {2 p} A _ {2 p} \right| \quad \text {a t} \omega_ {2 3} \equiv 2 \omega_ {p}, \\ a _ {2 4} \equiv \left| G _ {n - p} A _ {n - p} + G _ {m - p} A _ {m - p} \right| \quad \text {a t} \omega_ {2 4} \equiv \omega_ {n} - \omega_ {p}, \\ a _ {3 1} \equiv \left| G _ {n + m + p} A _ {n + m + p} + G _ {2 n + p} A _ {2 n + p} + G _ {2 m + p} A _ {2 m + p} \right| \quad \text {a t} \omega_ {3 1} \equiv 2 \omega_ {n} + \omega_ {p}, \\ a _ {3 2} \equiv \left| G _ {n + m - p} A _ {n + m - p} + G _ {2 n - p} A _ {2 n - p} + G _ {2 m - p} A _ {2 m - p} \right| \quad \text {a t} \omega_ {3 2} \equiv 2 \omega_ {n} - \omega_ {p}, \\ a _ {3 3} \equiv \left| G _ {3 n} A _ {3 n} + G _ {3 m} A _ {3 m} + G _ {n + 2 m} A _ {n + 2 m} + G _ {2 n + m} A _ {2 n + m} \right| \quad \text {a t} \omega_ {3 3} \equiv 3 \omega_ {n}, \\ a _ {3 4} \equiv \left| G _ {n + 2 p} A _ {n + 2 p} + G _ {m + 2 p} A _ {m + 2 p} \right| \quad \text {a t} \omega_ {3 4} \equiv \omega_ {n} + 2 \omega_ {p}, \\ a _ {3 5} \equiv \left| G _ {n - 2 p} A _ {n - 2 p} + G _ {m - 2 p} A _ {m - 2 p} \right| \quad \text {a t} \omega_ {3 5} \equiv \omega_ {n} - 2 \omega_ {p}, \\ a _ {3 6} \equiv | G _ {3 p} A _ {3 p} | \quad \text {a t} \omega_ {3 6} \equiv 3 \omega_ {p}. \end{array} \right\} (4. 1 8)
$$

Throughout this investigation, we consider the case defined by $h = 1 . 0 ~ \mathrm { m }$ , $h \kappa _ { n } =$ $h \kappa _ { m } \simeq 1 . 4 8$ and $h \kappa _ { p } \simeq 2 . 0$ . The first-order amplitudes are defined at a relatively low nonlinearity $( a _ { n } / h = a _ { m } / h = 0 . 0 1 5$ , $a _ { p } / h = 0 . 0 3 0 )$ ), while the angle $\varphi$ defining the short-crested waves is varied from 50 to $7 2 ^ { \circ }$ (with $9 0 ^ { \circ }$ corresponding to the direction of the plane wave). We emphasize that the current vector $U$ is set to zero throughout the calculations (i.e. an Eulerian drift is allowed). 

Whenever it makes physical sense, we use third-order boundary conditions. However, close to resonance we use second-order boundary conditions. The reason is that the perturbation theory will produce infinite third-order transfer functions (i.e. amplitudes) close to their singularities, while the second-order solution will remain finite. Hence, while we can use the third-order theory to predict when resonance can be expected, the third-order transfer functions are not valid in this region. The second-order boundary conditions therefore provide the highest-order basis to investigate the resulting nonlinear harmonic evolution. With second-order conditions, the dispersion relation for the input is linear, which leads to the input frequencies $\omega _ { n } = \omega _ { m } = 3 . 6 1 7 7 7 \ s ^ { - 1 }$ and $\omega _ { p } = \bar { 4 } . 3 4 9 0 5 \ \mathrm { s ^ { - 1 } }$ . Nonlinearities occurring in the system may then slightly modify $\kappa _ { n }$ and $\kappa _ { p }$ but approximately we have that $\rho \equiv \kappa _ { n } / \kappa _ { p } \simeq 0 . 7 4$ . In this case (4.15)–(4.17) predict resonance to occur at approximately $\varphi \simeq 7 0 ^ { \circ }$ in the component $\omega _ { 3 5 }$ . 

Figure 5 depicts computed harmonic amplitudes based on $\varphi = 5 0 ^ { \circ }$ , which is far from resonance. In this case the theoretical second harmonics are 20–60 times smaller than the first harmonics, while the third harmonics are another 10–20 times smaller than the second harmonics, except for $a _ { 3 5 }$ , which is only three times smaller than the largest second harmonic. We note that all 12 harmonics generated at the western boundary remain constant throughout the domain, i.e. showing no sign of harmonic modulation. Hence, for this case, the third-order boundary condition is in good agreement with the nonlinearity inherent in the model system. 

Figure 6 depicts computed harmonic amplitudes based on $\varphi = 6 5 ^ { \circ }$ , which brings us closer to the theoretical resonance in $\omega _ { 3 5 }$ . This change increases the theoretical magnitude of the third harmonic $a _ { 3 5 }$ by a factor of five, while the remaining third harmonics and the second harmonics are practically unchanged. As this is by far the most significant third-order component, it will be the only one shown on the remaining plots. It is interesting to note that the amplitude of $a _ { 3 5 }$ has now become somewhat 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/39b42d372ae1f3b9eb6a510f410983813b30a279d75f72daea176adbbef2e00f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/56d1f200ca140103b9c463f19f4ecb5c1a7e3eaadf36b7f21ce1ec6b968b0104.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/7ce31fdb23cf1f1829796a0e7bfbbf3f159970779b4b02d72645e3f2bef3d20c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/d9f95e6ca34d006a074b59c8167dcc8540ec6753e5964efbe8780b77e1e45fa7.jpg)



FIGURE 5. Spatial evolution of harmonic amplitudes computed by a high-order Boussinesq model; $\varphi = 5 0 ^ { \circ }$ (third-order input): $( a )$ first harmonics; $( b )$ second harmonics; $( c )$ third harmonic $\omega _ { 3 5 } = \omega _ { n } - 2 \omega _ { p }$ ; and $( d )$ remaining third harmonics. Specifications: $h = 1 . 0 \mathrm { m }$ , $a _ { n } = a _ { m } = 0 . 0 1 5 \mathrm { m }$ , $a _ { p } \doteq 0 . 0 3 0 \mathrm { m }$ , $h \kappa _ { n } = h \bar { \kappa _ { m } } = 1 . 4 8$ , $h \kappa _ { p } = 2 . \bar { 0 0 }$ , $\pmb { k } _ { n } = \kappa _ { n } ( \sin \varphi , \cos \varphi )$ , $k _ { m } = \kappa _ { m } ( \sin \varphi , - \cos \varphi )$ ), $\pmb { k } _ { p } = \kappa _ { p } ( 1 , 0 )$ .


larger than the largest of the second harmonics, and strictly speaking this violates the basic assumptions for the perturbation solution. Nevertheless, we notice only a moderate modulation of the amplitude of $a _ { 3 5 }$ throughout the wave tank, indicating that the third-order input is still in reasonable agreement with the nonlinearity inherent in the model system, despite the fact that the transfer function for $G _ { n - 2 p }$ is close to a singularity at this point. The spatial variation of the first and second harmonics shows little sign of modulation. 

Figure 7 also covers $\varphi = 6 5 ^ { \circ }$ , but now only second-order boundary conditions have been applied. This triggers a modulation in all third harmonics (again only $a _ { 3 5 }$ is shown) and the effect carries over to the second harmonics, while the first harmonics 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/4d89ccd57a65600612551c8be9e5c7bdaabfe5d773ba17b00422be1f75337884.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/d3eb34ca4a2956664bd384412ad36adef4225fd6079f01eeae899c29087edd51.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/4b6b6c54c9906bd387026c32b837e09c560b304707e5d7079460442a5ef2c12f.jpg)



FIGURE 6. Spatial evolution of harmonic amplitudes computed by a high-order Boussinesq model; $\varphi = 6 5 ^ { \circ }$ (third-order input): $( a )$ first harmonics; $( b )$ second harmonics; and $( c )$ third harmonic $\omega _ { 3 5 } = \omega _ { n } - 2 \omega _ { p }$ .


remain almost constant. In comparison with figure 6, it can be concluded that the third-order boundary condition leads to a much more stable harmonic evolution than if we apply only second-order conditions. This conclusion supports the perturbation theory even though its basic assumptions have been violated (with $a _ { 3 5 }$ being larger than the largest of the second harmonics). 

Figure 8 depicts computed harmonic amplitudes based on $\varphi = 6 9 ^ { \circ }$ , i.e. very close to the theoretical harmonic resonance. In this case the theoretical prediction of $a _ { 3 5 }$ is unrealistically high, and the third-order theory is now useless for providing boundary conditions. Instead, we run the model with second-order boundary conditions, and this expectedly triggers a strong modulation in $a _ { 3 5 }$ . The pattern is not fundamentally different from that shown in figure 7 except that the local peak values and beat lengths of $a _ { 3 5 }$ have grown significantly. Consequently, this effect now carries over to the second and first harmonics, which are no longer stable and constant. 

Figure 9 depicts the computed harmonic amplitudes based on the three cases (a) $\varphi = 7 0 ^ { \circ }$ , $( b )$ $\varphi = 7 1 ^ { \circ }$ and (c) $\varphi = 7 2 ^ { \circ }$ , surrounding the theoretical harmonic resonance. Again we utilize second-order boundary conditions. The plots show the resulting spatial variation of $a _ { 1 1 }$ , $a _ { 1 2 }$ and $a _ { 3 5 }$ , the latter being zero at the boundary. In all three cases we notice a strong modulation in $a _ { 3 5 }$ , which clearly carries over to $a _ { 1 1 }$ and $a _ { 1 2 }$ . The most extreme case appears to be $\varphi = 7 1 ^ { \circ }$ (figure $9 b$ ), which we shall discuss in further detail in the following. Near the input boundary, the wave 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/3a4faede6a4e015dc645220278c118a8eb40addca442d81e82f4d14d74bcab38.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/4efc5ac6f727151a78142dd03870b72ccd211c0e92b953720092b4f0fbdcd437.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/0862a71214268f186a03e302157d7309c654c9fe3bfc43574c7d3ef0202dee8b.jpg)



FIGURE 7. Spatial evolution of harmonic amplitudes computed by a high-order Boussinesq model; $\varphi = 6 5 ^ { \circ }$ (second-order input): $( a )$ first harmonics; $( b )$ second harmonics; and $( c )$ third harmonic $\omega _ { 3 5 } = \omega _ { n } - 2 \omega _ { p }$ .


train is dominated by the plane wave $a _ { 1 2 }$ and the short-crested wave $a _ { 1 1 }$ , while $a _ { 3 5 }$ is almost zero. Figure $1 0 ( a )$ shows a perspective plot of the corresponding surface elevation for the segment $0 \leqslant x / h \leqslant 1 5$ . According to figure $9 ( b )$ , $a _ { 3 5 } / h$ gradually grows, owing to harmonic resonance, to reach a maximum of approximately 0.03 in the region $x / h \simeq 4 0 0$ . Simultaneously, the amplitude of the plane wave $a _ { 1 2 } / h$ decays to a minimum of 0.01 in the same region, while $a _ { 1 1 } / h$ grows slightly to reach 0.038. This implies that, in the region near $x / h \simeq 4 0 0$ , the surface elevation will be dominated by the original short-crested wave $a _ { 1 1 }$ and the new resonating short-crested wave $a _ { 3 5 }$ , while the original plane wave $a _ { 1 2 }$ significantly decays. Figure $1 0 ( b )$ shows a perspective plot of the corresponding surface elevation covering the segment of $3 8 5 \leqslant x / h \leqslant 4 0 0$ . The difference between figures $1 0 ( a )$ and $1 0 ( b )$ is remarkable. 

# 4.3. Estimating the beat length of $a _ { 3 5 }$

In all cases shown in figures 6–9, the pattern of $a _ { 3 5 }$ looks like harmonic modulation rather than the growth of a singularity, and we shall therefore try to estimate the beat lengths in the following and compare them with the numerical results. 

First of all, the bound third harmonic is defined by 

$$
\omega_ {3 5} = 2 \omega_ {p} - \omega_ {n}, \quad \boldsymbol {k} _ {3 5} = 2 \boldsymbol {k} _ {p} - \boldsymbol {k} _ {n}, \tag {4.19}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/bd5504c69ada528a1562159a42563d802e9e0e3d26c4cf1351195e270434fbe7.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/378f137ae2386f028d935a6e7aa5af4c53d952a1d17bd525f811ea7eb526edf0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/0cf288fcdc3aacb990dfbe57cd5c8a36a8f86bef3ef17d3f590582f1c5e28cae.jpg)



FIGURE 8. Spatial evolution of harmonic amplitudes computed by a high-order Boussinesq model; $\varphi = 6 9 ^ { \circ }$ (second-order input): $( a )$ first harmonics; $( b )$ second harmonics; and $( c )$ third harmonic $\omega _ { 3 5 } = \omega _ { n } - 2 \omega _ { p }$ .


where the primary waves are given by 

$$
\boldsymbol {k} _ {p} \equiv \kappa_ {p} (1, 0), \quad \boldsymbol {k} _ {n} \equiv \kappa_ {n} (\sin \varphi , \cos \varphi), \quad \boldsymbol {k} _ {m} \equiv \kappa_ {n} (\sin \varphi , - \cos \varphi). \tag {4.20}
$$

Hence it is defined entirely in terms of the free waves $( \omega _ { p } , k _ { p } )$ and $( \omega _ { n } , k _ { n } )$ , which both satisfy the dispersion relation. By forcing the component $\omega _ { 3 5 }$ to zero at the input boundary, a spurious free wave $( \omega _ { s } , k _ { s } )$ of the same frequency $\omega _ { s } = \omega _ { 3 5 }$ will be released. Since the spurious wave is a free wave, its wavenumber $k _ { s }$ will satisfy the dispersion relation, in contrast to $\pmb { k } _ { 3 5 }$ , which will not. 

Next, the y component of $k _ { s }$ will be identical to that of the bound wave because of the the physical constraints in the $y$ direction (with fully reflective lateral boundaries). Hence we get 

$$
k _ {s y} = 2 k _ {p y} - k _ {n y}. \tag {4.21}
$$

This implies that the $x$ component of $k _ { s }$ can be determined by 

$$
k _ {s x} = \sqrt {\kappa_ {s} ^ {2} - k _ {s y} ^ {2}} \quad \text {w h e r e} \kappa_ {s} \equiv | \boldsymbol {k} _ {s} |. \tag {4.22}
$$

Typically, $k _ { s x }$ (the spurious wave) will be different from the $x$ component of $\pmb { k } _ { 3 5 }$ (the bound wave), and therefore a spatial modulation in the $x$ direction can be expected 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/8a8c98720a366e22b88cb5b6b1290b0c98793473c8a7ea17ab577c964d20f7f6.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/b72ab564b1bacfb77f177c051935798163322f6fbe5c2fb7ce5e0cea9b7b1df0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/ca028adf117b3d92bdd6c5913adc0f318e0cf8c505400375784e2d043fc76587.jpg)



FIGURE 9. Spatial evolution of the harmonic amplitudes $a _ { 1 1 }$ , $a _ { 1 2 }$ and $a _ { 3 5 }$ . Computed by a high-order Boussinesq model with second-order input: (a $) \varphi = 7 0 ^ { \circ }$ ; (b $) \varphi = { \bar { 7 } } 1 ^ { \circ }$ ; and (c) $\varphi = 7 2 ^ { \circ }$ .


of $a _ { 3 5 }$ . The resulting beat length can thus be estimated by 

$$
x _ {\text {b e a t}} = \frac {2 \pi}{\left| k _ {s x} - k _ {3 5 x} \right|} \quad \text {w h e r e} k _ {3 5 x} \equiv 2 \kappa_ {p} - \kappa_ {n} \sin \varphi . \tag {4.23}
$$

To utilize (4.23) and compare it with numerical simulations, we may consider various options for determining the relevant wavenumber vectors. 

1. Option 1 is to assume that $\kappa _ { n }$ , $\kappa _ { p }$ and $\kappa _ { s }$ can all be estimated by the linear dispersion relation (3.5). 

2. Option 2 is to assume that $\kappa _ { n }$ and $\kappa _ { p }$ obey the third-order amplitude dispersion relation (3.45) with (3.73) involving the interaction between the primary components $\pmb { k } _ { n }$ , $k _ { m }$ and $k _ { p }$ , while the spurious wave is treated as linear. 

3. Option 3 is to assume that the spurious wave is infinitesimal, but that its wavenumber $k _ { s }$ participates in the nonlinear interaction with the three primary waves $\pmb { k } _ { n }$ , $k _ { m }$ and $k _ { p }$ . In comparison with option 2, this will have no effects on $\pmb { k } _ { n }$ , $k _ { m }$ and $k _ { p }$ but it could have a strong effect on $k _ { s }$ and therefore also on the beat length. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/2f37e80f2adad83038d7f7c99e5603f06118659e37a44c371092f90fa812ed3a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-20/06102069-08f2-4baf-8a9b-e193299921e2/4d5e34f1c2aed2b5eca8f7a6ff428c29be93076164f221644e4bee5a62e9ef57.jpg)



FIGURE 10. Segments of the computed free surface elevation for the case of $\varphi = 7 1 ^ { \circ }$ : (a) $0 \leqslant x / h \leqslant 1 5$ ; and (b) $3 8 5 \leqslant x / h \leqslant 4 0 0$ .


4. Option 4 is to allow $a _ { s }$ to be finite, which would influence also the determination of $\pmb { k } _ { n }$ , $k _ { m }$ and $k _ { p }$ . However, this option would require that the magnitude of the spurious wave is known or can be estimated from the beginning. 

We have utilized options 1, 2 and 3 and determined the corresponding beat lengths from (4.23). Table 1 shows the predictions compared to the numerical simulations for a range of wave angles (previously studied in figures 4–8). The individual predictions will typically become very large when $\varphi$ approaches the point of singularity. We notice that options 1, 2 and 3 reach this point for $\varphi = 7 0$ , 69 and $7 1 ^ { \circ }$ , respectively. As expected, option 3 generally leads to the most accurate prediction of the beat length, and except for $\varphi = 7 0 – 7 1 ^ { \circ }$ it is relatively accurate. This implies that to predict the beat length most accurately near harmonic resonance, the wavenumbers should generally be determined by the third-order dispersion relation, allowing for the nonlinear interaction between $\pmb { k } _ { n }$ , $k _ { m }$ , $k _ { p }$ and $k _ { s }$ . 

# 5. Summary and conclusions

A new third-order solution for multi-directional irregular water waves in finite water depth has been presented. The solution provides explicit expressions for the surface elevation, the amplitude dispersion and the vertical variation of the velocity potential. Specific expressions for the velocity potential at the free surface are also given. The formulation incorporates the effect of an ambient arbitrarily directed current with the option of specifying zero net volume flux. 

<table><tr><td>φ (deg.)</td><td>Simulated</td><td>Option 1</td><td>Option 2</td><td>Option 3</td></tr><tr><td>60</td><td></td><td>38</td><td>42</td><td>36</td></tr><tr><td>65</td><td>78</td><td>82</td><td>103</td><td>72</td></tr><tr><td>66</td><td></td><td>105</td><td>141</td><td>89</td></tr><tr><td>67</td><td></td><td>142</td><td>217</td><td>114</td></tr><tr><td>68</td><td>147</td><td>215</td><td>449</td><td>156</td></tr><tr><td>69</td><td>200</td><td>425</td><td>15722</td><td>243</td></tr><tr><td>70</td><td>352</td><td>6780</td><td>443</td><td>523</td></tr><tr><td>71</td><td>725</td><td>511</td><td>230</td><td>5317</td></tr><tr><td>72</td><td>400</td><td>252</td><td>157</td><td>457</td></tr><tr><td>75</td><td></td><td>107</td><td>86</td><td>132</td></tr></table>


TABLE 1. Predicted and simulated beat lengths $L _ { b e a t } / h$ for component $a _ { 3 5 }$ (Eulerian drift is allowed).


The new third-order transfer functions for trichromatic interactions basically require the input of the amplitude, the frequency and the direction of (at least) three independent linear waves on a given water depth. By setting one of these amplitudes to zero, the formulation automatically simplifies to the third-order formulation for bichromatic waves given by Madsen & Fuhrman (2006). If we additionally use the same frequency for the remaining two first-order waves, the formulation further simplifies to describe the third-order monochromatic short-crested waves considered by Chappelear (1961) and Hsu et al. (1979). Through the elegant use of multipurpose functions, the new formulation is of the same complexity as the bichromatic formulation by Madsen & Fuhrman (2006). 

Harmonic resonance may occur at third order for certain combinations of frequencies and wavenumber vectors. This situation has been analysed for a monochromatic short-crested wave interacting with a plane wave with a different frequency, and the resulting resonance curves have been shown for deep water as well as finite-depth interactions. In the case of harmonic resonance, the perturbation theory will break down due to inherent singularities in the transfer functions. However, even when the theory fails to predict the wave amplitudes near harmonic resonance, it can still provide a fairly accurate estimate of the interacting wavenumbers via the third-order amplitude dispersion. 

Utilizing generation boundary conditions based on the new theory, we have made long-term simulations with a high-order Boussinesq model to study the evolution of wave trains both away from and near harmonic resonance. Far from resonance, the amplitudes are steady, providing confirmation of the theory. Near resonance, secondorder boundary conditions have been applied, and consequently harmonic modulation will occur in the third and higher harmonics. The modulation is particularly strong for the third-order components, which are close to resonance, and in extreme cases it has been demonstrated that their amplitudes may grow to become as large as the first-order components. However, all the observed patterns show a recurrence from zero to a peak and then back to zero, and the corresponding beat lengths can be estimated by the third-order theory using various options. The best option has been to assume that the spurious wave (the resonating component) is infinitesimal, but that its wavenumber $k _ { s }$ participates in the nonlinear interaction with the three primary finite-amplitude waves with wavenumbers $\pmb { k } _ { n }$ , $k _ { m }$ and $k _ { p }$ . This estimate is fairly accurate as long as the spurious wave is not too large. 

<table><tr><td>Gn+m</td><td>1.21227</td><td>Fn+m</td><td>0.06725</td><td>μn+m</td><td>-3.2100</td></tr><tr><td>Gn-m</td><td>0.23110</td><td>Fn-m</td><td>-0.15114</td><td>μn-m</td><td>-0.2644</td></tr><tr><td>Gn+p</td><td>1.84576</td><td>Fn+p</td><td>-0.00177</td><td>μn+p</td><td>-4.2063</td></tr><tr><td>Gn-p</td><td>-0.02151</td><td>Fn-p</td><td>-1.06833</td><td>μn-p</td><td>-1.4616</td></tr><tr><td>Gm+p</td><td>2.08857</td><td>Fm+p</td><td>-0.01257</td><td>μm+p</td><td>-4.5563</td></tr><tr><td>Gm-p</td><td>-0.05519</td><td>Fm-p</td><td>-1.39583</td><td>μm-p</td><td>-1.6396</td></tr><tr><td>G2n</td><td>2.2099</td><td>F2n</td><td>-0.05665</td><td>μ2n</td><td>-4.8517</td></tr><tr><td>G2m</td><td>2.2496</td><td>F2m</td><td>-0.03275</td><td>μ2m</td><td>-4.8383</td></tr><tr><td>G2p</td><td>2.3213</td><td>F2p</td><td>-0.01744</td><td>μ2p</td><td>-4.8699</td></tr><tr><td>Gn+m+p</td><td>5.7970</td><td>Fn+m+p</td><td>0.01363</td><td>μn+m+p</td><td>-14.6932</td></tr><tr><td>Gn-m+p</td><td>-19.4274</td><td>Fn-m+p</td><td>8.3119</td><td>μn-m+p</td><td>43.6845</td></tr><tr><td>Gn+m-p</td><td>11.0085</td><td>Fn+m-p</td><td>-17.635</td><td>μn+m-p</td><td>-31.1801</td></tr><tr><td>Gn-m-p</td><td>-20.1654</td><td>Fn-m-p</td><td>-5.2989</td><td>μn-m-p</td><td>-42.1139</td></tr><tr><td>Gn+2m</td><td>4.4984</td><td>Fn+2m</td><td>0.02081</td><td>μn+2m</td><td>-12.3356</td></tr><tr><td>Gn-2m</td><td>-10.0604</td><td>Fn-2m</td><td>-2.3172</td><td>μn-2m</td><td>-22.1368</td></tr><tr><td>G2n+m</td><td>4.3213</td><td>F2n+m</td><td>0.02734</td><td>μ2n+m</td><td>-12.0861</td></tr><tr><td>G2n-m</td><td>-6.2325</td><td>F2n-m</td><td>2.1495</td><td>μ2n-m</td><td>15.4132</td></tr><tr><td>Gn+2p</td><td>8.2574</td><td>Fn+2p</td><td>0.008435</td><td>μn+2p</td><td>-19.0664</td></tr><tr><td>Gn-2p</td><td>-46.2707</td><td>Fn-2p</td><td>-12.2461</td><td>μn-2p</td><td>-93.2046</td></tr><tr><td>G2n+p</td><td>7.7497</td><td>F2n+p</td><td>0.01389</td><td>μ2n+p</td><td>-18.3511</td></tr><tr><td>G2n-p</td><td>-14.9513</td><td>F2n-p</td><td>10.212</td><td>μ2n-p</td><td>37.873</td></tr><tr><td>Gm+2p</td><td>9.7798</td><td>Fm+2p</td><td>0.005915</td><td>μm+2p</td><td>-21.6745</td></tr><tr><td>Gm-2p</td><td>-84.3545</td><td>Fm-2p</td><td>-32.1242</td><td>μm-2p</td><td>-177.486</td></tr><tr><td>G2m+p</td><td>9.4711</td><td>F2m+p</td><td>0.007542</td><td>μ2m+p</td><td>-21.2328</td></tr><tr><td>G2m-p</td><td>-45.472</td><td>F2m-p</td><td>30.1335</td><td>μ2m-p</td><td>108.768</td></tr><tr><td>G3n</td><td>3.2570</td><td>F3n</td><td>0.004151</td><td>μ3n</td><td>-7.3561</td></tr><tr><td>G3m</td><td>3.4665</td><td>F3m</td><td>0.002496</td><td>μ3m</td><td>-7.6360</td></tr><tr><td>G3p</td><td>3.7803</td><td>F3p</td><td>0.001326</td><td>μ3p</td><td>-8.0819</td></tr></table>


TABLE 2. Transfer functions for the third-order trichromatic interaction. Note that $G$ is dimensionless, while $F$ and $\mu$ have the dimension $\mathrm { ~ m ~ s ~ } ^ { - 1 }$ .


<table><tr><td>\( \Upsilon_{nn} \)</td><td>0.6312</td><td>\( \Xi_{nn} \)</td><td>0.5333</td><td>\( \Omega_{nn} \)</td><td>0.5796</td></tr><tr><td>\( \Upsilon_{mm} \)</td><td>0.5361</td><td>\( \Xi_{mm} \)</td><td>0.7036</td><td>\( \Omega_{mm} \)</td><td>0.5571</td></tr><tr><td>\( \Upsilon_{pp} \)</td><td>0.4529</td><td>\( \Xi_{pp} \)</td><td>0.8938</td><td>\( \Omega_{pp} \)</td><td>0.5395</td></tr><tr><td>\( \Upsilon_{nm} \)</td><td>0.2540</td><td>\( \Xi_{nm} \)</td><td>1.3647</td><td>\( \Omega_{nm} \)</td><td>0.5486</td></tr><tr><td>\( \Upsilon_{mn} \)</td><td>0.1406</td><td>\( \Xi_{mn} \)</td><td>1.5546</td><td>\( \Omega_{mn} \)</td><td>0.6422</td></tr><tr><td>\( \Upsilon_{np} \)</td><td>0.8033</td><td>\( \Xi_{np} \)</td><td>1.5626</td><td>\( \Omega_{np} \)</td><td>0.7386</td></tr><tr><td>\( \Upsilon_{pn} \)</td><td>0.4317</td><td>\( \Xi_{pn} \)</td><td>2.0059</td><td>\( \Omega_{pn} \)</td><td>1.0228</td></tr><tr><td>\( \Upsilon_{mp} \)</td><td>0.8573</td><td>\( \Xi_{mp} \)</td><td>1.7494</td><td>\( \Omega_{mp} \)</td><td>0.9000</td></tr><tr><td>\( \Upsilon_{pm} \)</td><td>0.6402</td><td>\( \Xi_{pm} \)</td><td>1.9925</td><td>\( \Omega_{pm} \)</td><td>1.0646</td></tr></table>


TABLE 3. Important coefficients for the third-order trichromatic interaction. Note that $\varOmega$ is dimensionless, while $\boldsymbol { \Upsilon }$ and $\boldsymbol { \Xi }$ have the dimension $\mathrm { m ^ { - 1 } ~ s ^ { - 1 } }$ .


# Acknowledgement

We acknowledge the Danish Center for Scientific Computing for providing supercomputing resources. 

# Appendix

Tables 2 and 3 provide the relevant coefficients corresponding to figure 1 and they are given to aid in checking any implementation of the theory. 

# R E F E R E N C E S



BADULIN, S. I., SHRIRA, V. I., KHARIF, C. & IOUALALEN, M. 1995 On two approaches to the problem of instability of short-crested waves. J. Fluid Mech. 303, 297–325. 





BRYANT, P. J. 1985 Doubly periodic progressive permanent waves in deep water. J. Fluid Mech. 161, 27–42. 





CHAPPELEAR, J. E. 1961 On the description of short-crested waves. Tech. Memo no. 125, Beach Erosion Board, U.S. Army Corps of Engineers. 





CRAIG, W. & NICHOLLS, D. P. 2002 Travelling gravity water waves in two and three dimensions. Eur. J. Mech. (B/Fluids) 21, 615–641. 





FUHRMAN, D. R. & BINGHAM, H. B. 2004 Numerical solutions of fully nonlinear and highly dispersive Boussinesq equations in two horizontal dimensions. Intl J. Numer. Meth. Fluids 44, 231–255. 





FUHRMAN, D. R. & MADSEN, P. A. 2006 Short-crested waves in deep water: a numerical investigation of recent laboratory experiments. J. Fluid Mech. 559, 391–411. 





FUHRMAN, D. R., MADSEN, P. A. & BINGHAM, H. B. 2006 Numerical simulation of lowest-order short-crested wave instabilities. J. Fluid Mech. 563, 415–441. 





HAMMACK, J. L., HENDERSON, D. M. & SEGUR, H. 2005 Progressive waves with persistent, two-dimensional surface patterns in deep water. J. Fluid Mech. 532, 1–51. 





HAMMACK, J. L., SCHEFFNER, N. & SEGUR, H. 1989 Two-dimensional periodic waves in shallow water. J. Fluid Mech. 209, 567–589. 





HASSELMANN, K. 1962 On the nonlinear energy transfer in a gravity-wave spectrum. Part 1. General theory. J. Fluid Mech. 12, 481. 





HENDERSON, D. M., PATTERSON, M. S. & SEGUR, H. 2006 On the laboratory generation of two-dimensional, progressive, surface waves of nearly permanent form in deep water. J. Fluid Mech. 559, 413–427. 





HSU, J. R. C., TSUCHIYA, Y. & SILVESTER, R. 1979 Third-order approximation to short-crested waves. J. Fluid Mech. 90 (1), 179–196. 





IOUALALEN, M. 1993 Fourth order approximation of short-crested waves. C. R. Acad. Sci. Paris 316 (II), 1193–1200. 





IOUALALEN, M. & KHARIF, C. 1993 Stability of three-dimensional progressive gravity waves on deep water to superharmonic disturbances. Eur. J. Mech. B 12 (3), 401–414. 





IOUALALEN, M. & KHARIF, C. 1994 On the subharmonic instabilities of steady three-dimensional deep water waves. J. Fluid Mech. 262, 265–291. 





IOUALALEN, M., KHARIF, C. & ROBERTS, A. J. 1999 Stability regimes of finite depth short-crested water waves. J. Phys. Oceanogr. 29, 2318–2331. 





IOUALALEN, M., OKAMURA, M., CORNIER, S., KHARIF, C. & ROBERTS, A. J. 2006 Computation of short-crested deepwater waves. ASCE J. Waterway Port Coastal Ocean Engng 132 (3), 157–165. 





IOUALALEN, M., ROBERTS, A. J. & KHARIF, C. 1996 On the observability of finite-depth short-crested water waves. J. Fluid Mech. 322, 1–19. 





JANSSEN, P. A. E. M. 2003 Nonlinear four-wave interactions and freak waves. J. Phys. Oceanogr. 33, 863–884. 





JANSSEN, P. A. E. M. & ONORATO, M. 2007 The intermediate water depth limit of the Zakharov equation and consequences for wave prediction. J. Phys. Oceanogr. 37, 2389–2400. 





KIMMOUN, O., IOUALALEN, M. & KHARIF, C. 1999 Instabilities of steep short-crested surface waves in deep water. Phys. Fluids 11 (6), 1679–1681. 





KRASITSKII, V. P. 1994 On reduced equations in the Hamiltonian theory of weakly nonlinear surface waves. J. Fluid Mech. 272, 1–20. 





LONGUET-HIGGINS, M. S. 1962 Resonant interactions between two trains of gravity waves. J. Fluid Mech. 12, 321–332. 





MADSEN, P. A., BINGHAM, H. B. & LIU, H. 2002 A new Boussinesq method for fully nonlinear waves from shallow to deep water. J. Fluid Mech. 462, 1–30. 





MADSEN, P. A., BINGHAM, H. B. & SCHAFFER ¨ , H. A. 2003 Boussinesq-type formulations for fully nonlinear and extremely dispersive water waves – derivation and analysis. Proc. R. Soc. Lond. A 459, 1075–1104. 





MADSEN, P. A. & FUHRMAN, D. R. 2006 Third-order theory for bichromatic bi-directional water waves. J. Fluid Mech. 557, 369–397. 





MADSEN, P. A., FUHRMAN, D. R. & WANG, B. H. 2006 A Boussinesq-type method for fully nonlinear waves interacting with a rapidly varying bathymetry. Coast. Engng. 53, 487–504. 





MARCHANT, T. R. & ROBERTS, A. J. 1987 Properties of short-crested waves in water of finite depth. J. Austral. Math. Soc. B 29, 103–125. 





MCLEAN, J. W. 1982 Instabilities of finite-amplitude gravity waves on water of finite depth. J. Fluid Mech. 114, 331–341. 





PHILLIPS, O. M. 1960 On the dynamics of unsteady gravity waves of finite amplitude. Part I. J. Fluid Mech. 9, 193–217. 





ROBERTS, A. J. 1981 The behaviour of harmonic resonant steady solutions to a model differential equation. Q. J. Mech. Appl. Maths 34, 287–310. 





ROBERTS, A. J. 1983 Highly nonlinear short-crested water waves. J. Fluid Mech. 135, 301–321. 





ROBERTS, A. J. & PEREGRINE, D. H. 1983 Notes on long-crested water waves. J. Fluid Mech. 135, 323–335. 





ROBERTS, A. J. & SCHWARTZ, L. W. 1983 The calculation of nonlinear short-crested gravity waves. Phys. Fluids 26 (9), 2388–2392. 





SHARMA, J. & DEAN, R. 1981 Second-order directional seas and associated wave forces. Soc. Petrol. Engng J. 21, 129–140. 





SMITH, D. H. & ROBERTS, A. J. 1999 Branching behaviour of standing waves – the signatures of resonance. Phys. Fluids 11 (5), 1051–1065. 





STIASSNIE, M. & GRAMSTAD, O. 2009 On Zakharov’s kernel and the interaction of non-collinear wavetrains in finite water depth. J. Fluid Mech. 639, 433–442. 





STOKES, G. G. 1847 On the theory of of oscillatory waves. Trans. Camb. Phil. Soc. 8, 441–455. 





ZAKHAROV, V. E. 1968 Stability of periodic waves of finite amplitude on the surface of a deep fluid. J. Appl. Mech. Tech. Phys. 9, 190–194. 





ZAKHAROV, V. 1999 Statistical theory of gravity and capillary waves on the surface of a finite depth fluid. Eur. J. Mech. (B/Fluids) 18, 327–344. 





ZHANG, J. & CHEN, L. 1999 General third-order solutions for irregular waves in deep water. J. Engng Mech. 125 (7), 768–779. 

