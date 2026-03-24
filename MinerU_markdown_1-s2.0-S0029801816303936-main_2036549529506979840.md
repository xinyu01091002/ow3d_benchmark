# Numerical wave tank study of extreme waves and wave-structure interaction using OpenFoam®

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/c9e0142f0c11407cce7453309081b963917588f401e09e3e12b38fecdf5afd3b.jpg)


sm 

Zheng Zheng $\mathrm { H u } ^ { * }$ , Deborah Greaves, Alison Raby 

School of Marine Science and Engineering, Faculty of Science and Environment, University of Plymouth, Plymouth, Devon PL4 8AA, UK 

# A R T I C L E I N F O

Keywords: 

OpenFoam® 

Wave2Foam 

Focused waves 

NewWave 

VOF phase-fraction 

Wave-structure interaction 

# A B S T R A C T

In the present work, the open source Computational Fluid Dynamics (CFD) package-Open Field Operation and Manipulation (OpenFoam®) is used to simulate wave-structure interactions and a new wave boundary condition is developed for extreme waves. The new wave boundary condition is implemented for simulation of interaction with a fixed/floating truncated cylinder and a simplified Floating Production Storage and Offloading platform (FPSO) and results are compared with physical experiment data obtained in the COAST laboratory at Plymouth University. Different approaches to mesh generation (i.e. block and split-hexahedra) are investigated and found to be suitable for cases considered here; grid and time convergence is also demonstrated. The validation work includes comparison with theoretical and experimental data. The cases performed demonstrate that OpenFoam® is capable of predicting these cases of wave-structure interaction with good accuracy (e.g. the value of maximum pressure on the FPSO is predicted within $2 . 4 \%$ of the experiment) and efficiency. The code is run in parallel using high performance computing and the simulations presented have shown that OpenFoam® is a suitable tool for coastal and offshore engineering applications, is able to simulate two-phase flow in 3D domains and to predict wave-structure interaction well. 

# 1. Introduction

As is well known, coastal and offshore structures whether a sea wall, oil and gas platform, wave energy device or other fixed or floating structure, must be designed to survive in a very hostile environment, including heavy storms. For example, an extremely high or steep wave impact on the hull of a moored FPSO (Floating Production Storage and Offloading) platform may result in damage due to the impact or to water on the deck. Known as green water, this may cause severe damage to the deck house or other deckside equipment. Thus, there is a need for simulation tools to predict impact loading and to provide more information of non-linear wave-structure interaction. 

The phenomenon of freak or rogue waves and their interaction with structures has been an active research area in recent years (see Walker et al. 2004). The first measurement of such a wave in location was the New Year wave, which was observed at the Draupner platform in the North Sea off the coast of Norway on 1st January (1995) (see Adcock and Taylor (2009) and Adcock et al., 2011). The high order components of such nonlinear waves can excite offshore structures and cause high frequency responses. This was first encountered during the tow out of large platforms and so-called ringing of offshore structure columns was then demonstrated in laboratory experiments by Chaplin et al. (1997). Gaps in understanding of extreme wave-structure 

interaction was further evidenced by unexpected damage to Schiehallion FPSO vessel in 1998 during the build-up to a modest storm. Much research has been carried out in this area to consider local forces on offshore structures, such as the SAFE-FLOW European project (duration: 2000–2003) to consider wave loads and the CresT joint industry project (duration: 2008–2011) to consider realistic extreme wave profiles. A recent research project considering extreme wave loading on offshore wave energy devices using CFD (see Westphalen et al. 2014) presented four different CFD approaches, including AMAZON-SC 3D code (see Hu et al., 2009, 2010), the Control Volume Finite Element (CV-FE) and the Finite Volume (FV) (see Westphalen et al. 2012) and the Smoothed Particle Hydrodynamics (SPH) (see Omidvar et al. 2013) methods. The project found that CFD modelling offers predictions of high accuracy and a high density of flow field data within a reasonable timescale to assist with design changes and improvements to coastal or offshore structures. 

Wave loading on cylinders has been widely studied in coastal and offshore engineering because many offshore platforms and jetties are pile-supported. Morison et al. (1950) proposed a semi-empirical and semi-theoretical equation to calculate the wave force for small-scale problems (where the ratio of diameter to wave length $D / L { \le } 0 . 1 5$ , in which $D$ is the cylinder diameter and $L$ is the wave length). MacCamy 

and Fuchs (1954) proposed a linear diffraction theory for vertical surface-piercing cylinder, which was extended to nonlinear wave loads by the Faltinsen (1993) and Newman (1996) theory. Semi-analytical solutions have been proposed for particular geometries such as bottommounted vertical cylinders (Chau and Eatock Taylor, 1992, Malenica and Molin 1995, Newman 1996). More recently, second-order theories (Kriebel, 1998, Rahman and Heaps, 1983 and Taylor and Hung 1987), which consider the nonlinearity and can give more accurate wave force predictions, were proposed. In addition, some numerical methods (Wu and Hu, 2004, Ma et al. 2001a, 2001b) based on the potential flow formulation have been developed for this problem. Furthermore, Zhang et al. (2009, 2012) simulated an oscillating water column device in waves by solving the Navier Stokes equations with the level set method for the free surface, and Boccotti et al. (2012), Boo (2002, 2006), Dixon et al. (1979), Chaplin et al. (1997) and Kriebel (1998) have also investigated wave forces on structures. 

In this paper, OpenFoam® (hereafter OF) has been adopted and modified to study highly nonlinear free surface waves and wave interaction with structures. A major advantage is that OF is an open source code, it is possible to gain control over the exact implementations of different features and to develop and implement new models and fit them into the overall code structure. At present, a number of applications based on OF have been published for the numerical simulation of the wave-structure interactions with various offshore structures and ranges of wave conditions. For example, Gerald et al. (2010) used the RasInterFoam, part of the OF library, for non-linear wave interaction with a cylinder; Chen et al. (2014) focused on the assessment of how OpenFOAM performs for wave interactions with a vertical surface cylinder, which compared with physical experiments under three regular waves and first-order focused wave groups; Higuera et al. (2014) used IHFOAM, a new implementation that can reduce drastically CPU time, for simulating a porous coastal structure under oblique incidence of irregular waves on a high mound breakwater; Jiang et al. (2015) investigated tsunami-like solitary waves propagating over a row of vertical -slotted-cylindrical piles; Li et al. (2012a, 2012b) simulated the interaction between wave and ship body in a tank sloshing and compared the numerical and experimental results of the agreement well with each other, including impact load on tank wall, wave forces acting on ship hull and the ship motion. It can be seen that OpenFOAM is very capable of accurate modelling of nonlinear wave interaction with offshore structures. 

Regarding coastal engineering applications, and those cases in which extreme free surface motion is important, OF supports two phase flow simulation; for example, the OF solver InterFoam/ interDyMFoam with the waves2Foam library has been used in this study to investigate nonlinear extreme free surface flows. In this paper, our effort is added an extreme wave boundary condition with second order components for implementation in OF within the waves2Foam library in order to simulate both challenging wave climates and also average sea states. The extreme wave formulation prescribed as an inlet condition follows Dalzell (1999) and Hu et al. (2011, 2014), which is based on a focused wave group generated using the second order Stokes wave theory. Extreme wave boundary conditions are described by Ning et al. (2009) for their fully nonlinear potential flow simulation and compared with experiments in a wave flume. Ning et al. (2008, 2009) studied the propagation of NewWave groups (see Tromans et.al. 1991) for 4 different wave heights up to the point where the waves almost break, and results from physical experiments and the numerical study are found to be in good agreement. NewWave simulations and comparison with the experiment data (see Ning et al. 2009) were also considered by Hu et al. (2009 and, 2011) using the Finite Volume Method (FVM) and by Westphalen et al. (2012) using the FVM and the CV-FE, for which good agreement was achieved. The numerical wave flume is based on the Madsen and Sørensen (1992) set of Boussinesq equations with focused wave groups and the non-linear shallow water equations. The numerical wave tank is based on the Boussinesq 

equations with second order focused wave groups and compared experiment for non-breaking waves, and the nonlinear shallow water equations for broken waves. Good agreement is achieved between the numerical predictions and laboratory measurements of free surface elevation, run-up distances and overtopping volumes for the test case with a plane beach and a seawall in a wave basin (see Orszaghova et al. 2014). In the OF model used here, a Volume of Fluid (VOF) interface capturing approach is taken for the free surface. Its principal advantages are that it is very simple, allowing very complex free surface configurations to be represented easily and that it involves no mesh motion. The incompressible Navier-Stokes equations are solved using the FVM on unstructured meshes (block mesh or split hexahedra mesh). The laminar flow assumption is selected in this study because for the cases considered here, which do not include wave breaking, turbulence is unlikely to have a significant effect. 

The method of parallel computing used by OF is known as domain decomposition, in which the geometry and associated domain are broken into sub-domains and allocated to separate processors for solution. The process of parallel computation involves: decomposition of mesh and domain, running the application in parallel and postprocessing the decomposed case as described in the user guide of the Open source CFD toolbox (http://www.openfoam.org). The parallel runs use the public domain openMPI implementation of the standard Message Passing Interface (MPI) to manage parallel distributed memory computers within acceptable CPUs. In this paper, first the case of a regular wave generated by a Numerical Wave Tank (NWT) has been investigated by using a code parallelisation implemented on a desktop and a cluster of high performance computing (hereafter HPC), which shows remarkable improvements in performance of the method. The parallel computing used is then applied to the remaining test cases. 

The aims of this work are to validate the OF models using the new extreme wave boundary condition and to provide new data to improve understanding of the propagation pattern using extreme waves induced by a transient wave groups. This paper is organised as follows. In Section 2, the numerical method used is presented and the extreme wave condition and solution procedure are described in Section 3. Validations and analysis are carried out in Section 4, which include: a transient wave groups first and second order Stokes waves; a fixed horizontal cylinder and a fixed vertical cylinder under regular waves, an extreme wave generated using first and second order focused wave groups; a fixed/floating truncated cylinder and a simplified FPSO under second order focused wave groups. Conclusions are provided in Section 5. 

# 2. Numerical method in OpenFoam

# 2.1. Governing equations

Consider two fluids (air and water) in a computational domain. The fluids are separated by an interface (free surface). Both air and water are assumed to be incompressible in the present study, thus the governing equations for the incompressible laminar fluid are as follows: 

$$
\frac {\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf {U}) = 0, \quad \text {a n d} \tag {1}
$$

$$
\frac {\partial \rho \mathbf {U}}{\partial t} + \nabla \cdot (\rho \mathbf {U} \mathbf {U}) - \nabla \cdot (\mu (\nabla \mathbf {U} + (\nabla \mathbf {U}) ^ {T}) = \rho g - \nabla p \tag {2}
$$

where $\pmb { U }$ denotes the velocity vector, $\rho$ the density, μthe dynamic viscosity, and $g$ the acceleration due to gravity. 

An additional equation must also be solved to describe the movement of the phases. The indicator phase function $\alpha$ is defined as the quantity of water per unit of volume in each cell. This means that if $\scriptstyle { a = 1 }$ the cell is full of water, if $\scriptstyle a = 0$ the cell is full of air, and in any other case it belongs to the air-water interface. It is straightforward to calculate 

any of the properties of the fluid at each cell, just by weighting them by the VOF function. For example, the fluid density and the dynamic viscosity of the cell are computed as follows: 

$$
\rho = \alpha \rho_ {\text {w a t e r}} + (1 - \alpha) \rho_ {\text {a i r}}
$$

$$
\mu = \alpha \mu_ {\text {w a t e r}} + (1 - \alpha) \mu_ {\text {a i r}} \tag {3}
$$

The starting point for the equation which tracks the fluid movement is an equation for advance the phase fraction field: 

$$
\frac {\partial \alpha}{\partial t} + \nabla \cdot (\mathbf {U} \alpha) = 0 \tag {4}
$$

OF makes use of an artifical compression term $\boldsymbol { \nabla } { \cdot } \mathbf { U } _ { \mathbf { c } } ( 1 - \alpha )$ (see Weller, 2002), which added in the phase equation (4). $\mathbf { U _ { c } }$ is the artifical compressive velocity and equal to $\mathbf { \left| U _ { c } \right| } = \operatorname* { m i n } { \left[ c _ { \alpha } \mathbf { | U | } \right. }$ U , max( ) ] , where the factor $c _ { \alpha }$ usually takes value 1. Therefore, the governing differential equation for $\alpha$ from Rusche (2002) is: 

$$
\frac {\partial \alpha}{\partial t} + \nabla \cdot (\mathbf {U} \alpha) + \nabla \cdot (\mathbf {U} _ {\mathrm {c}} \alpha (1 - \alpha)) = 0. \tag {5}
$$

The function $\alpha$ is calculated using this equation by means of a specially designed solver called MULES (Multidimensional Universal Limiter for Explicit Solution). It makes use of a limiter factor on the fluxes of the discretised divergence term to ensure a final value between 0 and 1. 

# 2.2. InterFoam/InterDyMFoam solver and waves2Foam library

Versions V.2.2.0/2.3.0 of OF have been used in this work. The interFoam solver is prepared for static meshes only. The interDyMFoam solves the same as interFoam equations but it can handle dynamic meshes (‘DyM’ stand for Dynamic Mesh) for simulating floating body. InterFoam/InterDyMFoam is one of the solvers included in OF and it can solve the three-dimensional Navier-Stokes equations for two incompressible phases using a finite volume discretisation and the VOF method. In the VOF method, each phase is described by a fraction $\alpha _ { i }$ occupied by the volume of fluid of ith material in the cell. Its principal advantages are that it is very simple, allowing very complex free surface configurations to be represented easily, and that it involves no mesh motion. 

The solver algorithm used by InterFoam/InterDyMFoam is called PIMPLE, and is a combination of PISO (Pressure Implicit with Splitting of Operators) and SIMPLE (Semi-Implicit Method for pressure-Linked Equations) algorithms. Its main structure is inherited from the original PISO, but it allows equation under-relaxation to ensure the convergence of all the equations at each time step. Both algorithms are thoroughly explained in applications with VOF by Jasak (1996). 

The library waves2Foam is a toolbox used to generate and absorb free surface water waves. The relaxation zone may be implemented to work simultaneously with wave generation at the wave inlet or to absorb waves only at the outlet. This feature is a key point for coastal engineering as it allows for a shorter computational domain to be used. The passive wave absorption method is used as it appears in Schäffer and Klopman (2000) and recently presented by Lara et al. (2011) and Jacobsen et al. (2012). Relaxation zones at both inlet and outlet (see Fig. 1) have been used for the NWT in this paper. The present relaxation technique is an extension to that of Mayer et al. (1998) and the relaxation function 

$$
\alpha_ {R} \left(\chi_ {R}\right) = 1 - \frac {\exp \left(\chi_ {R} ^ {3 . 5}\right) - 1}{\exp (1) - 1} \quad \text {f o r} \chi_ {R} \in [ 0: 1 ], \tag {6}
$$

is applied inside the relaxation zone in the following way, 

$$
\phi = \alpha_ {R} \phi_ {\text {c o m p u t e d}} + (1 - \alpha_ {R}) \phi_ {\text {t a r g e t}} \tag {7}
$$

where $\phi$ is either Uor $\alpha$ . The variation of $a _ { R }$ is the same as given by 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/090853ed96c9e38be224b815d794655cae005660a42a678801b4169b656d4536.jpg)



Fig. 1. The relaxation zone at both inlet and outlet.


Fuhrman et al. (2006). The definition of $\chi _ { R }$ is such that $a _ { R }$ is always 1 at the interface between the non-relaxed part of the computational domain and the relaxation zone, as illustrated in Fig. 1. 

A large number of wave theories are supported in Waves2Foam, including three regular wave generation types defined by Stokes first, second and fifth order, as well as cnoidal, stream function, Boussinesq solitary and irregular waves. Full details are provided in Higuera et al. (2013) and the OpenFoam website (http://openfoamwiki.net/index.php/Contrib/waves2Foam). In this work, the NewWave theory with second order components has been developed and implemented in OF to account for extreme waves and full details are given in next section. 

# 3. Implementation of new boundary conditions

As is well known, the exact velocity profile for a true physically realisable nonlinear wave under given conditions is not known a priori. Thus, a viable approach is to input reasonable approximate wave conditions along the input boundary to simulate the real phenomenon. This leads to the notion of the extreme wave formulation as a focused wave group in which many wave components in a spectrum are focused simultaneously at a position in space in order to model the average shape of an extreme wave profile consistent with the random process in a specified wave energy spectrum (see Tromans et al. 1991). The formulation here refers to the work by Dalzell ( 1999) and Hu et al. (2011, 2014) in which a first or second-order Stokes focused wave can be imposed in such a manner. 

# 3.1. NewWave formulation

A Cartesian coordinate system $\boldsymbol { O } - \boldsymbol { x } \boldsymbol { y } \boldsymbol { z }$ is defined with the origin located at the undisturbed equilibrium free surface, with the zcoordinate vertical and positive upwards. The $x$ -coordinate is zero at the wave-maker located at $x = 0 . 0 m$ , $x _ { 0 }$ is the focus point, t is the focus time and the water depth h. $A _ { i }$ is the input wave amplitude of the focused wave. 

The corresponding wave elevation $\eta$ , and horizontal and vertical velocities $u$ and $w$ are expressed as follows: 

$$
\eta = \eta^ {(1)} + \eta^ {(2)} \tag {8}
$$

$$
u = u ^ {(1)} + u ^ {(2)} \tag {9}
$$

$$
w = w ^ {(1)} + w ^ {(2)} \tag {10}
$$

where $\eta ^ { ( 1 ) } , u ^ { ( 1 ) }$ and $\boldsymbol { w } ^ { ( 1 ) }$ are the linear wave elevation and velocities, $\eta ^ { ( 2 ) } , u ^ { ( 2 ) }$ and $w ^ { ( 2 ) }$ correspond to the second-order wave elevation and velocities, respectively. Both velocity and wave elevation can be decomposed into $N$ components with different frequencies following Hu et al. (2011, 2014) and are included in Appendix A for completeness. 

For the simulations presented here, the incoming wave entering the computational domain is fluxed through the inlet boundary. This flux is defined in terms of either first order theory, or first order theory plus second order theory. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/82de5b86149872e711abeeed76c20daa8512487af70696a04907deb187a2886a.jpg)



Fig. 2. Flow chart of OpenFoam sequences.


# 3.2. Solving procedure

Implementation of the boundary condition is dependent on the user interface of the software package. OF uses its own script languages in C $^ { + + }$ to express mathematical equations and logical operations. As described earlier, the approach in OF is valid for regular waves and some irregular waves (i.e. bichromatic) with the waves2Foam library, however, the NewWave-based equations (see Appendix A) have to be expressed as a new boundary condition for extreme waves. The OF flow chart incorporating both solver and the new extreme boundary condition is shown in Fig. 2. The OF solver starts with the preprocessor, which is needed to set up wave properties at the initial conditions. When no other wave data is available, the sea climate may be described by using a design spectrum. The JONSWAP and Pierson-Moscovitz spectra are commonly used in offshore engineering. In this paper, input data included the wave amplitude of $A _ { i } ( f )$ and frequency of $f$ (i) are used in the experiment only. The wave number $k ( i ) = { \omega _ { i } } ^ { 2 } / g \operatorname { t a n h } ( k ( i ) h )$ , and phase-lag, $p h i ( i ) = k ( i ) \times x _ { 0 } - 2 \pi f ( i ) \times t _ { 0 } .$ , are calculated as individual parameters and together are input to the pre-processor. Subsequently, OF is called once, before the time step and updates are made during the time step loop. For the NewWave inlet boundary condition, the corresponding wave elevation $\eta$ , horizontal $u$ and vertical velocities $w$ are calculated (see Eqs. (9), (10) and (11)) based on the solution at time t and the $z$ coordinate at every time step. It is an incoming NewWave condition entering the NWT through a transparent inlet boundary. These boundary values are calculated as described before. 

# 4. Numerical simulation and results

In the following simulations, the water is assigned density $\rho { = } 1 0 0 0 \mathrm { k g / m ^ { 3 } }$ and kinematic viscosity $\nu { = } 1 0 ^ { - 5 } \mathrm { P a } / \mathrm { s }$ , while the air assigned density $\rho { = } 1 ~ \mathrm { k g } / \mathrm { m } ^ { 3 }$ and kinematic viscosity $\nu { = } 1 . 4 8 { \times } 1 0 ^ { - 4 } \mathrm { P a } /$ s. A low Reynolds number is specified for all cases (e.g. $R e = d | U | / \nu = 9 7 5$ in the case of a vertical cylinder in regular waves). The following boundary conditions were applied. At the inlet, the velocities and surface elevation $\eta$ are specified; the velocity is specified for the water component only and the velocity of the air at the inlet boundary is set to zero. The pressure is set to zero normal gradients at all boundaries and at the outlet zero gradients condition on velocity is applied. The top boundary and right far boundary are specified with a non-reflecting boundary condition allowing air to leave or enter the domain. The remaining boundaries and structure are set as rigid walls with no-slip boundary conditions. The force calculation is obtained by integration of the pressure and viscous force components around the body contour as follows: 

$$
F _ {\text {p r e s s u r e}} = \int_ {S _ {b}} p \mathbf {n} d A, \tag {11}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/c94150decc2df88c49dfe64b1ea010b7e4f78b180b620e226f877f803f0ad344.jpg)



Fig. 3. Comparison of surface elevation for the regular wave condition.


$$
F _ {\text {v i s c o u s}} = \int_ {S _ {b}} \left(\tau_ {x x} + \tau_ {y x} + \tau_ {z x}\right) d x + \left(\tau_ {x y} + \tau_ {y y} + \tau_ {z y}\right) d y + \left(\tau_ {x z} + \tau_ {y z} + \tau_ {z z}\right) d z, \tag {12}
$$

where $S _ { b }$ is the body surface and $\tau _ { x x } , \tau _ { y x } , \ldots \tau _ { \mathrm { z z } }$ are nine viscous stress components (see Versteeg and Malalaskera, 1995). 

# 4.1. Regular wave NWT

In the first test case, the 3D NWT is verified by simulating regular waves of amplitude $a { = } 0 . 1 \mathrm { m }$ , regular wave number of $k { = } 1 . 2 8 \ \mathrm { m } ^ { - 1 }$ , the steepness $k a { = } 0 . 1 2 8$ , the relative water depth of $k h { = } 3 . 5 7$ , the wave length of $4 . 9 \mathrm { m }$ and wave period of 1.77 s. The NWT has outer dimensions $2 0 \mathrm { m } { \times } 0 . 0 8 \mathrm { m } { \times } 3 . 5 \mathrm { m }$ and water depth $\times .$ . A relaxation zone of $5 \mathrm { m }$ length is defined both at the inlet and outlet. A uniform mesh was applied in this case with 700,000 cells and mesh spacing of $0 . 0 2 \mathrm { m }$ . Fig. 3 shows a comparison between the surface elevation time history predicted by the OF simulation and theoretical data. Fig. 4 shows the comparison of the surface elevation over one wave period. The maximum crest elevation predicted in the second order Stokes theory is $1 . 8 \%$ higher than the first order theory. It can be seen that the wave elevation predicted by OF and the theory agree well. 

The OF simulations were run in parallel on a desktop PC and a HPC (High Performance Computing) cluster. The desktop PC is a Quad core $3 . 4 \ : \mathrm { G H z }$ with 16 Gb RAM with a maximum of 8 cores and the HPC cluster is a Quad core 2.56 GHz with 16 Gb RAM with a maximum number 372 cores. In this case, the total simulation time is 20 s. The time step of the numerical simulation is set to $0 . 0 0 1 s$ and the adjustable time step is set to limit the maximum Courant number to 0.25 for stability. Details of the CPU computational time for each of the simulations are summarised in Table 1. It can be seen that the parallelisation of openMPI works well and increasing the number of cores for the parallel computation can save execution time. 

# 4.2. A fixed horizontal cylinder in regular waves with Stokes’ first order theory

The first wave-structure interaction case considered is the interaction between regular waves and a horizontal cylinder in a tank, where the cylinder is positioned such that the axis is at $z { = } { \mathrm { - } } 0 . 0 7 5 ~ \mathrm { m }$ and the still water level at $z { = } 0 . 0 \ \mathrm { m }$ . The purpose of the test case is again to provide validation, this time of the wave forces on the cylinder compared with the theory based on the modified Morison's equation (see Dixon et al. 1979 and Morison et al. 1950) and experimental results (see Dixon et al. 1979). According to the physical experiments by Dixon et al. (1979) the wave signal is accurate to first order. Therefore, a regular wave of Stokes first order is generated in the NWT to interact with the cylinder. 

Test parameters, including the cylinder diameter $D$ , wave number $k$ , wave steepness $k a$ , relative water depth $k h$ , Keulegan-Carpenter number $N _ { k c } ,$ wave amplitude a and wave period $_ T$ are shown in Table 2. To compare the numerical result in OF with those obtained by Dixon et al. (1979) the vertical relative force $F$ ′on the cylinder is exported. Dimensionless parameters are defined as: relative force 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/42ec2b3bf961ca1e37d0398691335b7e627f367318d7db1631619e733b90a6be.jpg)



Fig. 4. Comparison of surface elevation over one wave period.



Table 1 Details of CPU consumption for each HPC.


<table><tr><td>Computer</td><td>Cores</td><td>Execution time</td></tr><tr><td rowspan="4">Desktop: Quad core 3.4 GHz with 16 Gb RAM</td><td>1</td><td>14 h 15 min</td></tr><tr><td>2</td><td>10 h 38 min</td></tr><tr><td>4</td><td>8 h 36 min</td></tr><tr><td>6</td><td>7 h 34 min</td></tr><tr><td rowspan="4">HPC cluster: Quad core 2.56 GHz with 16 Gb RAM</td><td>4</td><td>15 h 27 min</td></tr><tr><td>16</td><td>5 h 8 min</td></tr><tr><td>32</td><td>2 h 52 min</td></tr><tr><td>64</td><td>2 h 28 min</td></tr></table>


Table 2 Parameter set up: horizontal cylinder simulation.


<table><tr><td>D (m)</td><td>k</td><td>ka</td><td>kh</td><td>Nkc</td><td>a (m)</td><td>T (s)</td></tr><tr><td>0.25</td><td>1.61</td><td>0.01</td><td>1.61</td><td>1.3</td><td>0.05</td><td>1.65</td></tr></table>

$F ^ { ' } = F _ { z } / [ g \rho ( 1 / 4 \pi D ^ { 2 } L ) ]$ , relative amplitude $A ^ { ' } = a / D$ , relative wave elevation $\eta ^ { ' } = \eta / D$ , relative wavelength $L ^ { ' } = L / D$ $\scriptstyle { \mathrm { ~ \mathcal { L } ~ } }$ is the wavelength) and relative axis depth $d ^ { ' } = d / D$ and the vertical force $F _ { z }$ on the cylinder. The NWT geometry has outer dimensions $1 2 \mathrm { m } \times 0 . 2 \mathrm { m } \times 1 . 5 \mathrm { m }$ with a water depth $h = 0 . 4 5 m$ . The relaxation zone has a length of 2 m at both inlet and outlet boundaries. The cylinder sits one wavelength downstream of the inlet and is defined as a wall. 

A non-uniform (split-hexahedra) mesh was used in the NWT. The regions close to the cylinder have been refined with spacings $\Delta z = 0 . 0 1 \mathrm { ~ m ~ }$ , 0.015 m and $0 . 0 2 \mathrm { m }$ to check the accuracy of the NWT. The grid convergence index (GCI) was examined for the root mean square (RMS) of the relative vertical force over one wave period T. The value of $G C I _ { 3 2 }$ (where 3 indicates a coarse mesh and 2 the corresponding intermediate mesh in Fig. 5), is approximately $4 . 3 \%$ and $G C I _ { 2 1 }$ (where 2 again indicates an intermediate mesh and 1 the corresponding fine mesh 1), is $2 . 6 \%$ (see Roache et al. 1986). This confirms that the calculations are mesh convergent. 

The intermediate mesh, 2, with a spacing of $\Delta z = 0 . 0 1 5 m$ (giving 521,300 cells) was selected for the cylinder as shown in Fig. 6. The total vertical force $F _ { z }$ on the cylinder is shown in Fig. 7 for the entire 13 s simulation. It shows that the force time history reaches a steady state after the third wave period. In Fig. 8, the relative vertical force over one period is compared with the theoretical (see Dixon et al. 1979) and experimental data and shows that there is general agreement. In the experiment, the force was measured over one wave cycle once steady state was reached. The force measurements were accurate to within $1 \%$ of the largest force measured and the initial force on the cylinder in still water was subtracted before the force measurements were taken. The numerical result is taken from one wave period (see Fig. 7) once a steady state has been reached. In Fig. 8 where the axis of $t / T$ ranges from 0.0 to 0.5, the numerical and experimental curves show an asymmetry not predicted by the theoretical equation. This can be 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/13c50ab26af9b835203b7d21eb89595d8721a1d8ff6b379f4393cda2df092b53.jpg)



Fig. 5. Grid convergence study for relative vertical force prediction.


explained by the following. Initially the waterline is at $0 . 2 5 D$ below the axis of the cylinder. However, the theory is based on the cylinder being totally submerged and predictions are based on first order buoyancy effects only. When the cylinder becomes further immersed, the prediction shows better agreement with the theory as $t / T$ varies from 0.5 to 1.0. Fig. 9 shows snapshots of the surface elevation over one period around the horizontal cylinder. These snapshots show the free surface deformations corresponding to the relative vertical force presented in Fig. 8. The features of deformation at the free surface around the cylinder can be clearly seen, which generates a depressed wave region in front of the cylinder at $t / T { = } 0 . 1 2$ and $t / T { = } 0 . 3 6$ , increased wave elevation in front of the cylinder at $t / T { = } 0 . 6$ and the cylinder fully immersed in water at $t / T { = } 0 . 7 3$ . The waterline at $t / T { = } 1 . 0$ after one period is the same as the initial waterline at $t / T { = } 0 . 0$ . The simulations were carried out using four cores running on a desktop with a CPU times of $1 7  { \mathrm { h } } \ 2 7  { \mathrm { m i n } }$ with a maximum Courant number of 0.25 and simulation time of ${ 1 3 \mathrm { s } }$ . 

# 4.3. A fixed vertical cylinder in regular waves with Stokes’ second order theory

This case describes the numerical simulation of a fixed bottommounted cylinder in regular waves. The main purpose of this case is to provide a comparison between the OF prediction and theoretical models which are based on linear and second-order diffraction theory (see Kriebel, 1998), as well as reported experimental work (see Kriebel, 1998). To correspond with 2nd order diffraction theory, a regular wave of Stokes’ second order is generated in the NWT to interact with the cylinder. 

Test parameters including the cylinder diameter $D$ , the wave number $k _ { : }$ , the scattering parameter kr (r is the cylinder radius), the steepness ka $( a { = } 0 . 0 5 3 5 \mathrm { m }$ is the wave amplitude), the Keulegan-Carpenter number $N _ { k c } ,$ , wave amplitude a and wave period $T$ are shown in Table 3. The NWT has outer dimensions of $1 2 \ : m \times 1 . 2 \ : m \times$ $_ { 0 . 9 m }$ and a water depth equal to $h = 0 . 4 5 m$ . The cylinder it positioned at one wavelength $\phantom { - } ( 3 . 7 7 m )$ from the wave maker. The relaxation zone is defined to be $2 \mathrm { m }$ in length, both at the inlet and outlet. Two simulations were run based on the block mesh with 1,487,040 cells and the split-hexahedra mesh with 966,692 cells, for each case, the mesh close to the cylinder is shown in Figs. 10 and 11. Results for the two different meshes are shown in Fig. 12 where the total horizontal force $F _ { x }$ on the cylinder is plotted as a function of time. The agreement is very good showing that the solution is independent of the mesh type used. In Figs. 13 and 14, the contribution to the force due to pressure and viscous forces on the cylinder are plotted separately. It is noticeable that the viscosity has little contribution to the total force in comparison with pressure force. This is expected because the case has a small Keulegan-Carpenter number $( N _ { k c } { = } 0 . 3 2 )$ ), and so lies in the range of wavelength-to-characteristic body length ratio where diffraction effects dominate the loads and flow separation is not significant. Nevertheless it cannot be neglected as the viscous force is related to the vortex-flow around the cylinder and is of special interest in wave 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/67b85c6973191512d3e3e06955bfcfcf57ecbc42d2adf1866446c6eec78d6d77.jpg)



Fig. 6. Split-hexahedra mesh around cylinder.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/13a8f8920957443291b7beb8f25423efade331093dd0687fcbdc4fb0c676dc25.jpg)



Fig. 7. Total vertical force.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/963e381da500aa9494acab490ffc3db3ef4915c2d7d9912c81d1c0ed65f082f3.jpg)



Fig. 8. Relative vertical forces.


load estimation. It is necessary to incorporate all significant effects 


Table 3 Parameter set up: vertical cylinder simulation.


<table><tr><td>D (m)</td><td>k</td><td>kr</td><td>ka</td><td>Nkc</td><td>a (m)</td><td>T (s)</td></tr><tr><td>0.325</td><td>1.667</td><td>0.271</td><td>0.09</td><td>0.32</td><td>0.0535</td><td>1.95</td></tr></table>

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/b0ac85f5b317173682ec9a143fcdb2b97be6ac7f417269f2bd99e29ad548acf0.jpg)



Fig. 10. Block mesh around vertical cylinder (blockMesh).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/41130198acfbcf906eda7e05776fa3a56ed3d88ed9d635fc32c8995b7cf3ea5e.jpg)



Fig. 11. Split-hexahedra mesh around vertical cylinder (snappyHexMesh).


associated with the nonlinear free surface flow. The force is normalized in a standard format as $F _ { x } / F _ { 0 }$ F kh ρgaHh kD , ( = / tanh( ), where $H$ is the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/82fa4c881f89d669a42688b1e4bf83d2ef139e60144f80dab799179eebe9e2b0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/dd00a7d7252caff504a59d680f0e06419b4f0238bb952757838bb2f39984d792.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/7532f9bc992e74d1bf74885436e5a7201f64c1f8de396b99f3cf640b5eb31b01.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/148b65cd0a3da7498ef6feee1f4b07c5061545e484410f59e8d80f5b432b9d81.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/e1085a1efa20f09d11c2275c78b08bb6825ff2f5593440af6ee128fa459f8dc8.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/982f6796d7cdeee193f2af32ab8fea7cc142d2b63a5321a87b4cb649047a7785.jpg)



(f)



Fig. 9. Wave profile over one period (a) $\mathrm { { t } / \mathrm { { T } = 0 . 0 } }$ , (b) $\mathrm { t / T } { = } 0 . 1 2$ , (c) t/T=0.36, (d) $\mathrm { t } / \mathrm { T } { = } 0 . 6$ , (e) $\mathrm { t } / \mathrm { T } { = } 0 . 7 3$ , (f) $\mathrm { { t } / \mathrm { { T } = 1 . 0 } }$ .


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/ef54728bc61e6b4a8ea0907ac554c4608bb11fcf2d73945df466f37cda53299b.jpg)



Fig. 12. Comparison of total horizontal force on cylinder with different mesh types.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/983c9e0d2556bd83d9c3cc03966e1dc95c7071a4ac5e5b7362baa33d837933a8.jpg)



Fig. 13. Pressure forces on the cylinder in the horizontal direction.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/fe8208d1f12912fc920e4d10997d9166091bc157e0aa4762bd80f699151a2c28.jpg)



Fig. 14. Viscous forces on the cylinder in the horizontal direction.


wave height on the cylinder over one period as shown in Fig. 15. The significance of including second order terms in the theory is demonstrated in Fig. 15. Comparison of maximum magnitude of normalized force in the negative direction (opposing the direction of wave motion) is shown in Fig. 15, the experimental force exceeds the first order theory by $14 \%$ while the second order theory provides much better agreement and exceeds the experiment by $1 . 7 \%$ , as described by Kriebel (1998). Fig. 15 also illustrates the OF prediction with second order Stokes wave boundary condition, the OF prediction agrees with the experiment at the maximum force is predicted within $0 . 5 \%$ and the mean force is predicted within $1 5 \%$ . In comparison with second-order theory, OF agrees well; the maximum force is predicted within $0 . 4 \%$ and the mean force is predicted within $0 . 5 \%$ . The nonlinearities are the result of dominance of the second order plane wave components and are common features of nonlinear inertial forces experienced at small $k r < 0 . 4$ , as discussed by Kriebel (1998). A small difference in the result is expected because the physical data is averaged over 10 wave periods, whereas the numerical result represents one wave cycle only once a steady state has been reached. Fig. 16 shows snapshots of the surface elevation over one period around the vertical cylinder. These snapshots demonstrate the free surface deformations corresponding to relative the normalized force presented in Fig. 15. Features of the free surface flow around the cylinder include increased wave elevation at $t / T { = } 0 . 1 2$ and $t / T { = } 0 . 3 6$ , a depressed wave region at $t / T { = } 0 . 6$ and $t / T { = } 0 . 7 3$ and after one period the initial profile at $t / T { = } 0 . 0$ is recovered at $t / T { = } 1 . 0$ . Execution time is about 7 h 33 min on the block mesh and 1 d 4 h 54 min on the split-hexahedra mesh, with four cores running on a desktop with a maximum Courant number of 0.25 for a simulation time 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/e422b1ee6b9e5cc7577a47068eda1425da395fd17fe31a5b85a6a48a57b1617d.jpg)



Fig. 15. Comparison of the normalized horizontal force on cylinder over one wave period.


of ${ \boldsymbol { 1 2 s } }$ . It shows the CPUs consumption under the block mesh with regular hexahedra cells is less than split-hexahedra mesh, in which the face of some prism cells is not parallel to the surface of the cylinder. The split-hexahedra mesh needs more time to deal with additional terms to account for the non-orthogonality. 

# 4.4. Extreme NWT

Fig. 17 provides a schematic of the physical ocean basin and the NWT set up in OF for this case. The physical experiments were carried out in the COAST laboratory at Plymouth University. The experimental set-up, measurement devices used and the paddle signal used to generate focused waves in the ocean basin were described by Mai et al. (2015). The NewWave technique (see Tromans et al. 1991) was used to define the input to the physical model. The wave basin is $3 5 \mathrm { m } { \times } 1 5 . 5 \mathrm { m } { \times } 3 \mathrm { m }$ and the water depth in the experiment was $2 . 9 3 \mathrm { m }$ . Waves are generated by a flap-type wave maker with force controlled wave absorption and there is a parabolic beach to dissipate energy at the downstream end of the tank. The wave characteristics are shown in Table 4. The number of wave components, $N _ { : }$ , used in the experiment is 244 and the corresponding wave amplitude with frequency is shown in Fig. 18. 

In the numerical simulation, the length of the NWT is defined as 5 times the characteristic wave length, i.e. $5 \lambda _ { p } { = } 1 6 . 5 5 m$ , in which the length of the relaxation zone at inlet and outlet are $\lambda _ { p }$ and $2 \lambda _ { p }$ respectively. The water depth in the NWT is $h = 2 . 9 3 m$ , matching the physical experiment, and the initial height of air above the water surface is $_ { 0 . 5 \mathrm { m } }$ . The simulation is approximately 2D $( x$ and $y )$ , although actually the tank is one cell wide with the z- direction for numerical purposes. The vertical distribution of grid cells is chosen to be relatively coarse near the sea bed, but then becomes finer towards the free surface, which is suitable in deep water. The region from $y { = } { - } 0 . 2 m$ to $y { = } 0 . 2 m$ , which contains the free surface, has a uniform mesh. In the horizontal direction, the cell size is uniform across the whole domain. The linear focus position and focus time are defined as $1 . 5 \lambda _ { p }$ and $8 T _ { p }$ as suggested by Ning et al. (2009) in their numerical simulation. In the simulation, the focus time and position are slightly different because of nonlinear wave-wave interactions in the flume. In order to compare with the physical experiment, the focus position and time are determined to be where there is a maximum crest and two symmetric deep troughs. At the start of the computation, a cosine ramp function is applied over a wave period to prevent the impulse-like behavior of the wavemaker and reduce the corresponding unnecessary transient waves. 

To investigate the potential for reducing CPU time, numerical convergence was carried out on the number of wave components. Figs. 19 and 20 show the comparison of wave spectra and wave elevation at the focus location, which are obtained with the number of wave components $N$ separately defined as 20, 25 and 30. From these results, it can be seen that the results obtained are completely identical for the latter two numbers, indicating that convergence was achieved 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/d0fc5f25222a0a12d5dc82073f7a74bb3ca391c8b39fb2ec840338427cb11ded.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/d78d5f9a104121f28adb89237f5688e058d11f8b32857b9ad9e639278a10b928.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/823d0a3ba5dabe2c0199544c8b52675c7eb5f1e0db42ca83132d6a10c9ab95b0.jpg)



(c）


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/4b849208edd9e0d9f1d7d687d01503f813729a1f36396e5f24e4be6ceaece4d9.jpg)



(@)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/fc20bf52bd11ca9934db48db7b51d2329d312319b2063d39bffc6586e53473a3.jpg)



(e)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/8677cc994d22cd3f54f7b6aaa7fd1d41543714da21140411de77268255c8f5c6.jpg)



Fig. 16. Wave profile over one period (a) $\mathrm { { t } / \mathrm { { T } = 0 . 0 } }$ , (b) $\mathrm { t / T } { = } 0 . 1 2$ , (c) $\mathrm { t } / \mathrm { T } { = } 0 . 3 6$ , (d) $\mathrm { { t } / \mathrm { { T } = 0 . 6 } }$ , (e) t/T=0.73, (f) t/T=1.0.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/d1b1b5820eff09a28a95425749d094ffd94b4c1fc25deb53884af7ef1674d061.jpg)



Fig. 17. Schematic sectional drawing of physical ocean basin and the NWT setup.



Table 4 Wave properties.


<table><tr><td>Wave components</td><td>Wave period Tp(s)</td><td>Wave length λp(m)</td><td>Wave height Hs (m)</td><td>Frequency Band (Hz)</td></tr><tr><td>244</td><td>1.456</td><td>3.31</td><td>0.103</td><td>0.1–2.0</td></tr></table>

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/3da50bc8abd1777f3cb732af44739c462dce325e612bd7c1c5df0d156d981596.jpg)



Fig. 18. Wave amplitude spectra.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/834d3bd2450c735ee7c36ed0bb03e282c557037b2e095cb136002eb18fa8df84.jpg)



Fig. 19. Comparison of wave spectra for N=244. from numerical simulation.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/361e8876b1625e725806542edf903a9c01fe9631f803017469d2942b07829e8b.jpg)



Fig. 20. Comparison of surface elevation at focus location for different wave components.


using a component number $N { = } 2 5$ . 

To check the accuracy of the numerical tank, mesh and time convergence tests have been carried out. Results were obtained using three different uniform meshes with spacings $\varDelta x = 0 . 0 2 m$ (mesh 1) $, \Delta x = 0 . 0 4 m$ (mesh 2) and $\Delta x = 0 . 0 6 m$ (mesh 3) (see Fig. 21). The value of grid converge index of $G C I _ { 3 2 }$ is $7 . 0 \%$ and $G C I _ { 2 1 }$ is $2 . 5 \%$ . Fig. 22 shows the results of three different uniform meshes with spacings $\Delta y = 0 . 0 1 m$ (mesh 1), $\Delta y = 0 . 0 2 m$ (mesh 2) and $\varDelta y = 0 . 0 4 m$ (mesh 3) in the region that contains the free surface and 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/7193cd6a663c5d28f13e6a2bfc3125f40e8cbde4282fc1d4a987bc632f9b02b4.jpg)



Fig. 21. Comparison of surface elevation at focus location for different grid sizes in $x$ direction.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/356cd2ffd97cfb2503a24be549e5674bcbfde7372c4119318c3485db8df8ae18.jpg)



Fig. 22. Comparison of surface elevation at focus location for different grid sizes in $y$ direction.


the value of grid convergence index of $G C I _ { 3 2 }$ is $1 . 8 \%$ and $G C I _ { 2 1 }$ is $0 . 9 \%$ . These confirm that the calculations are mesh convergent. An adjustable time step is set by defining a maximum Courant number of 0.5, 0.25 and 0.1 under a time-step of $\varDelta t = 0 . 0 5 8 \mathrm { ~ s ~ }$ $( = T _ { p } / 2 5 )$ ). Comparison of results with three different Courant numbers and with mesh $\varDelta x = 0 . 0 4$ mand $\Delta y = 0 . 0 2 m$ are shown in Fig. 23. It can be seen that the results have converged with respect to time step. The execution time with maximum Courant number of 0.25 and 24,308 cells is about 19 min on 4 cores running on the desktop and 16 min on 32 cores running on the HPC cluster with simulation time of 18 s. It should be noted that there is no obvious reduction in execution time after increasing the number of cores beyond 32 for this case with a small total number of cells. The reason is due to the time taken in transferring data for every core which is a higher proportion of the execution time for smaller cases. 

In order to assess the required length of the NWT domain, simulations with different domain lengths of $5 \lambda _ { p } ( = 1 6 . 5 5 \mathrm { m } )$ $7 \lambda _ { p } ( = 2 3 . 1 7 \mathrm { m } )$ and $1 0 \lambda _ { p } ( = 3 3 . 1 0 \mathrm { m } )$ have been carried out. Fig. 24 shows the comparison of wave elevation at the focus location, which the results obtained are completely identical for the three simulations. Therefore, the length of the domain with $5 \lambda _ { p } \mathrm { i s }$ used in following numerical simulation. 

Fig. 25 illustrates the surface elevation time history at the focus location with first and second order wave generation in the numerical model and the comparison with experimental data. For this particular wave steepness, there is very little difference between the first and second order solutions and they overlay one another for much of the time history. At the crest, a small difference is evident and the second order wave boundary condition gives a slight improvement on the linear case, although both underpredict the maximum crest elevation measured in the experiment. The benefit of using the second order accurate boundary condition is more evident for steeper waves (Hu et al. 2014). Both numerical results of the trough elevation after 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/f45a2b0d0b2aa47cb757c5df9b56af33da1849a24baf8a1a66a89b56a19784e2.jpg)



Fig. 23. Comparison of surface elevation at focus location for different time steps.


maximum crest are lower and do not exactly coincide with the result of the experimental data. However, the match between numerical and physical experiment is generally very good and the second order NewWave boundary condition will be implemented on the following simulations for extreme wave generation. 

# 4.5. A fixed/floating truncated vertical cylinder in extreme waves

In this case, the physical experiments were carried out in the ocean basin within the COAST laboratory at Plymouth University. The wave characteristics, relaxation zone, the length and height of the domain are the same as in the extreme wave only case, though with the width of the domain set to $3 \mathrm { m }$ in a 3D simulation. The model is a truncated cylinder, for which the radius is equal to $0 . 1 5 \mathrm { m }$ and the height is $0 . 3 \mathrm { m }$ . Its initial submergence in the water is equal to $0 . 1 5 \mathrm { m }$ and the distance between the inlet and the front of the cylinder is $4 . 1 5 \mathrm { m }$ , equal to the focus location $x _ { r } = 4 . 1 5 \mathrm { m }$ . Figs. 26 and 27 show the splithaxahedra mesh around the cylinder in horizontal and vertical sections with a total 1,512,800 cells. It is noticeable that there are some twisted cells in Fig. 27. This is an aberration caused by the ParaView visualization tool in producing a vertical section around cylinder in 3D, in reality the snappyHexMesh algorithm in OF is robust with a prespecified final mesh quality. The simulation was 3 day 16 h CPU running time using 32 cores of a HPC cluster with a maximum Courant number of 0.25 for a simulation time ${ 1 8 \mathrm { s } }$ . 

Figs. 28 and 29 show the total vertical and horizontal force around the cylinder. There is non-zero initial vertical force in Fig. 28, which is the hydrostatic force only. The hydrostatic force is defined by $F = \rho A d g$ , which $\rho$ is the water density, $d$ is the depth below the water surface of the cylinder and A is the surface area of the underside of the cylinder, in which means the area of the bottom cylinder. Comparing the numerical and physical maximum crest elevations at the position of the front of the cylinder (see Fig. 30), the numerical prediction underestimates the experiment by $4 . 7 \%$ with a value of $0 . 1 1 8 \mathrm { m }$ compared with $0 . 1 2 4 \mathrm { m }$ measured in the experiment. It can be seen that this shows good agreement between the experimental data and numerical result. Fig. 31 presents the pressure measured at the front of the cylinder and $0 . 0 5 \mathrm { m }$ below the initial free surface; it agrees well with the experiment. The maximum pressure in the experiment is $1 . 3 6 2 \mathrm { k P a }$ , which the OF result underestimates by $2 . 4 \%$ with a value of $1 . 3 3 0 \mathrm { k P a }$ . Fig. 30 shows the wave run up on the front of the cylinder, and Fig. 31 shows the pressure time history for the pressure transducer at the front centerline of the cylinder, positioned $0 . 0 5 \mathrm { m }$ below the initial free surface. Emergence of the pressure transducer in the wave troughs can be seen from the wave elevation plot, and in the pressure time history, this corresponds with the gauge pressure falling to zero. In the numerical prediction, the pressure appears to rise above zero before the pressure transducer is the remerged. This may be due to the numerical errors. Fig. 32 shows the wave profile around the cylinder at $t { = } 1 0 . 4 8 \ s$ . It can be seen that the numerical simulation predicts the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/586b8190a27d501aa614348a56c3ec8266f71f1dfbe4ca33551b062e57e06c02.jpg)



Fig. 24. Comparison of surface elevation at focus location for different length of NWT domain.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/aeac9bb851d116bae746a82b5f06b7b4eafe25a0d37aeddb711e06bca83182c7.jpg)



Fig. 25. Comparison of surface elevation at focus location.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/c9bb913dfc52d7c301d30e81d71599e58d740e33ea5924db0beb7bf03a5cb562.jpg)



Fig. 26. Split-hexahedra mesh (horizontal Section) around cylinder.


emergence of the pressure transducer as the wave troughs pass. The result includes the incoming NewWave and the diffracted wave created by the wave interaction with the cylinder. 

Next, a floating cylinder in heave has been investigated, in which the cylinder is allowed to move in the z -direction only and other degrees of freedom are restricted. The wave characteristics and the geometry of the domain are the same as for the fixed cylinder described above. The boundary condition for the surface of the floating cylinder is defined as a moving wall, in which the cylinder is free to respond to the fluid excitation and its motion is calculated as part of the solution. The mass of cylinder is set to $m { = } 1 1 \mathrm { k g }$ , which matches the experiment. Fig. 33 shows the time history of $z$ displacement for the cylinder predicted by the numerical model and measured in the experiment. It is noticeable that the numerical prediction has a larger amplitude of motion than measured in the experiment and this may be due to 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/7e6f9e10861c1524eb0a460191683d695f9b7364a914e74a2bc6f8a2b07a42de.jpg)



Fig. 27. Split-hexahedra mesh (vertical Section) around cylinder.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/51234abea4671844aebb9c2e2deefff2d9f91f22b546a36b56a0d7690a13901c.jpg)



Fig. 28. Total vertical force on the cylinder.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/8e4f103b6c5e14758466d84d5847be00f57a84620c049717b556b0d8cb5368ef.jpg)



Fig. 29. Total horizontal force on the cylinder.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/fed27978a20be5fa019bf495a9c640c7cd5ce2359a494341c552b65f61f9ea80.jpg)



Fig. 30. Comparison of surface elevation at front of the cylinder.


friction in the heave-only support system used in the experiment that is not included in the numerical simulation. Fig. 34 shows the pressure predicted by the numerical model and measured at the front of the cylinder at $0 . 0 5 \mathrm { m }$ below the initial free surface. The numerical prediction is in reasonable agreement with the experiment in that the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/73a78d69b351df32082e5c12441d2185d45f7075db989ff235f35391b2509b9d.jpg)



Fig. 31. Comparison of pressure at the front of the cylinder along the centre-line at $0 . 0 5 \mathrm { m }$ below the initial free surface.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/c6ea7bae2ec3fa2b9195b698e1ec79957d678111aab90215672acb23cd85dbd5.jpg)



Fig. 32. Wave profile around the cylinder at $t { = } 1 0 . 4 8 \ s$ .


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/82d8b0518ae08a7f0cb5b5ced9600ef72fd27d41e8c61fb69043ecae3afb7d81.jpg)



Fig. 33. Comparison of z displacement of the cylinder.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/6843a9e4eea16f8b02daf6f40f9f38a34d1cde6111592a782b49e584b92a4276.jpg)



Fig. 34. Comparison of pressure at the front of the cylinder along the centre-line at $0 . 0 5 \mathrm { m }$ below the initial free surface.


frequency and phase of the pressure oscillation are predicted well. Emergence of the pressure transducer in the wave troughs can be seen from the wave elevation plot, and in the pressure time history, this corresponds with the gauge pressure falling to zero. 

# 4.6. A fixed simplified FPSO in extreme waves

The physical experiments were carried out in the COAST laboratory 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/2a6eced197d044d18e9f951184115208cf7be1338b397d5918c669ed4e6c6e66.jpg)



Fig. 35. Split-hexahedra mesh (horizontal section) around FPSO.


at Plymouth University. The experimental set-up, measurement devices used and the paddle signal used to generate focused waves in the ocean basin were described by Mai et al. (2015). The wave characteristics, relaxation zone, the length and height of the domain are the same as in extreme wave only case and the width of domain is set to $4 \mathrm { m }$ in a 3D simulation. The geometry of a simplified FPSO, which is formed by a rectangle and a semicircle at each end, is as follows: the radius of the semicircle is taken as $r { = } 0 . 1 5 \ \mathrm { m }$ , the total length of the FPSO is taken as $1 . 2 \mathrm { m }$ and height as $0 . 3 \mathrm { m }$ . Its submergence in the water is equal to $0 . 1 5 \mathrm { m }$ and the distance between the inlet and the bow of the FPSO is $4 . 1 5 \mathrm { m }$ , equal to the focus location of $4 . 1 5 \mathrm { m }$ . Figs. 35 and 36 show the split-haxahedra mesh around FPSO as horizontal and vertical sections with a total of 1,446,638 cells. The CPU time was 3 days $2 \mathrm { { h } }$ using 32 cores on the HPC cluster, with a maximum Courant number of 0.25 for a simulation time 18 s. 

Figs. 37 and 38 give total vertical and horizontal forces around the FPSO. In Fig. 37, the hydrostatic force at initial time is calculated by $F = \rho A d g$ , in which A is surface area of the underside of the FPSO. Comparing the numerical and physical maximum crest elevation at the bow of the FPSO (see Fig. 39), the numerical prediction underestimates the experiment by $6 . 8 \%$ with a value of $0 . 1 2 3 \mathrm { m }$ compared with $0 . 1 3 2 \mathrm { m }$ measured in the experiment. It can be seen that this shows reasonable agreement between the experimental data and numerical result. Fig. 40 presents the pressure measured at the front of the FPSO and $0 . 0 5 \mathrm { m }$ below the initial free surface; it agrees well with the experiment. The maximum pressure in the experiment is $1 . 2 9 2 \mathrm { k P a }$ , which the OF result overestimates by $2 . 6 \%$ with a value of $1 . 3 2 6 \mathrm { k P a }$ . Fig. 41 shows a wave profile around the FPSO at 10.48 s. The result includes the incoming NewWave and the diffracted wave created by the wave interaction with the FPSO. 

# 5. Conclusions

The goal of this paper is to develop numerical simulation of extreme waves and wave-structure interaction using OF. A new wave boundary condition is presented, which is based on NewWave for representation of the extreme wave event together with first or second-order Stokes wave theories for the individual wave components. The new boundary condition has been integrated in waves2Foam within OF as the inlet boundary condition. In the case of extreme wave NWT, mesh and time convergence tests to check the accuracy of the numerical tank are presented, and the number of wave components needed to represent properly the extreme wave has been investigated. After that, the second order NewWave boundary condition is implemented on a fixed/floating cylinder and a FPSO in extreme waves. The results are compared with physical experiments, which include the surface elevation and the pressure at front of the structure, and are in good agreement. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/161734626e53e8042cd38402142175530f3266d0243484f4fb6638e8810f187a.jpg)



Fig. 36. Split-hexahedra mesh (vertical section) around FPSO.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/c8fafe2f736785dbcdff22b643752934c93de8d9110c3f63ed88efea3f9ff662.jpg)



Fig. 37. Total vertical force on the FPSO.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/383f85a1f2a30d15a66075acf4448f33482d9f1ba98351f2a1c3e38b37469f13.jpg)



Fig. 38. Total horizontal force on the FPSO.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/420ebae8fc6552947a47330da28a14b1a500b929d156303d1886e6fe9b5a19a2.jpg)



Fig. 39. Comparison of surface elevation at the bow of the FPSO.


Regular wave interaction with a fixed horizontal and vertical cylinder have also been simulated in OF. The validation work includes comparison of the predicted wave force on a fixed cylinder with theoretical and experimental data and are in good agreement. In each case, the results allow details of the free surface deformation as the incident wave interacts with the cylinder. For the partially submerged horizontal cylinder, the free surface is seen to engulf the cylinder completely during the wave cycle. For the vertical cylinder, the maximum run up occurs as the wave crest passes the cylinder and the maximum horizontal wave force occurs at the same time. Two kinds of mesh generation (block and split-hexahedra) have been used, grid 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/0a66bd1ddfb96789f665e8c0c1f273dbfd6cd50ed8a390bdcaae22e7154cdb60.jpg)



Fig. 40. Comparison of pressure at the front of the FPSO along the centre-line at $0 . 0 5 \mathrm { m }$ below the initial free surface.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-25/0a7ae3ba-6e41-4f1a-a97b-be5e5eb4b535/f5d5caa380cb270c4f9ce83e76f56996b199e670092b3c01b489c19a4cf6c701.jpg)



Fig. 41. Wave profile around the FPSO at t=10.48 s.


and time convergence are demonstrated and show that OF is capable simulating structure interaction well and efficiently when code parallelisation is used. 

In general, the results confirm that the model in OF is well placed for extension to many coastal engineering applications to simulate a wide range of nonlinear wave conditions. Future work will include extension to wave interaction with floating and elastic structures under extreme wave conditions. 

# Acknowledgments

This work is supported by the UK Engineering and Physical Sciences Research Council (EPSRC), FROTH (Fundamentals and Reliability of Offshore sTructure Hydrodynamics) project of reference number: EP/J012866/1 for which the authors are most grateful. The authors would also like to acknowledge use of the experimental data supplied by the COAST laboratory at Plymouth University. 

# Appendix A

$$
\eta^ {(1)} = \sum_ {i = 1} ^ {N} A _ {i} \cos \left[ \left(k _ {i} \left(x - x _ {0}\right) - \omega_ {i} \left(t - t _ {0}\right) + \varepsilon_ {i} \right], \right. \tag {A1}
$$

$$
\begin{array}{l} \eta^ {(2)} = \sum_ {i = 1} ^ {N} \sum_ {j > i} ^ {N} \{A _ {i} A _ {j} B ^ {+} \cos [ (k _ {i} + k _ {j}) (x - x _ {0}) - (\omega_ {i} + \omega_ {j}) (t - t _ {0}) + (\varepsilon_ {i} + \varepsilon_ {j}) ] \\ + A _ {i} A _ {i} B ^ {-} \cos \left[ \left(k _ {i} - k _ {j}\right) \left(x - x _ {0}\right) - \left(\omega_ {i} - \omega_ {j}\right) \left(t - t _ {0}\right) + \left(\varepsilon_ {i} - \varepsilon_ {j}\right) \right] \rbrace \\ + \sum_ {i = 1} ^ {N} \frac {A _ {i} ^ {2} k _ {i}}{4 \tanh  (k _ {i} h)} (2 + \frac {3}{\sinh^ {2} (k _ {i} h)}) \cos [ 2 (k _ {i} (x - x _ {0}) - \omega_ {i} (t - t _ {0}) + \varepsilon_ {i}) ] \\ - \frac {A _ {i} ^ {2} k _ {i}}{2 \sinh (2 k _ {i} h)}, \tag {A2} \\ \end{array}
$$

$$
u ^ {(1)} = \sum_ {i = 1} ^ {N} \frac {g A _ {i} k _ {i}}{\omega_ {i}} \frac {\cosh \left(k _ {i} z\right)}{\cosh \left(k _ {i} h\right)} \cos \left[ k _ {i} \left(x - x _ {0}\right) - \omega_ {i} \left(t - t _ {0}\right) + \varepsilon_ {i} \right], \tag {A3}
$$

$$
\begin{array}{l} u ^ {(2)} = \sum_ {i = 1} ^ {N} \sum_ {j > i} ^ {N} \left\{A _ {i} A _ {j} A ^ {+} \left(k _ {i} + k _ {j}\right) \frac {\cosh \left(k _ {i} + k _ {j}\right) z}{\cosh \left(k _ {i} + k _ {j}\right) h} \cos \left[ \left(k _ {i} + k _ {j}\right) \left(x - x _ {0}\right) \right. \right. \\ - (\omega_ {i} + \omega_ {j}) (t - t _ {0}) + (\varepsilon_ {i} + \varepsilon_ {j}) ] + A _ {i} \mathrm {A _ {j}} \mathrm {A ^ {-}} (k _ {i} - k _ {j}) \frac {\cosh (k _ {i} - k _ {j}) z}{\cosh (k _ {i} - k _ {j}) h} \\ \cos \left[ (k _ {i} - k _ {j}) (x - x _ {0}) - (\omega_ {i} - \omega_ {j}) (t - t _ {0}) + (\varepsilon_ {i} - \varepsilon_ {j}) \right] \\ + \sum_ {i = 1} ^ {N} \frac {3 k _ {i} A _ {i} ^ {2} \omega_ {i}}{4} \frac {\cosh (2 k _ {i} z)}{\sinh^ {4} (k _ {i} h)} \cos [ 2 \left(k _ {i} \left(x - x _ {0}\right) - \omega_ {i} \left(t - t _ {0}\right) + \varepsilon_ {i}\right) ], \tag {A4} \\ \end{array}
$$

$$
\begin{array}{l} w ^ {(1)} = \sum_ {i = 1} ^ {N} \frac {g A _ {i} k _ {i}}{\omega_ {i}} \frac {\sinh (k _ {i} z)}{\cosh (k _ {i} h)} \sin \left[ k _ {i} \left(x - x _ {0}\right) - \omega_ {i} \left(t - t _ {0}\right) + \varepsilon_ {i} \right], \quad \text {a n d} (A5) \\ w ^ {(2)} = \sum_ {i = 1} ^ {N} \sum_ {j > i} ^ {N} \{A _ {i} A _ {j} A ^ {+} (k _ {i} + k _ {j}) \frac {\sinh (k _ {i} + k _ {j}) z}{\cosh (k _ {i} + k _ {j}) h} \sin [ (k _ {i} + k _ {j}) (x - x _ {0}) \\ - (\omega_ {i} + \omega_ {j}) (t - t _ {0}) + (\varepsilon_ {i} + \varepsilon_ {j}) ] + A _ {i} \mathrm {A _ {j}} \mathrm {A ^ {-}} (k _ {i} - k _ {j}) \frac {\sinh (k _ {i} - k _ {j}) z}{\cosh (k _ {i} - k _ {j}) h} \\ \sin \left[ \left(k _ {i} - k _ {j}\right) \left(x - x _ {0}\right) - \left(\omega_ {i} - \omega_ {j}\right) \left(t - t _ {0}\right) + \left(\varepsilon_ {i} - \varepsilon_ {j}\right) \right] \rbrace \\ + \sum_ {i = 1} ^ {N} \frac {3 k _ {i} A _ {i} ^ {2} \omega_ {i}}{4} \frac {\sinh (2 k _ {i} z)}{\sinh^ {4} (k _ {i} h)} \sin \left[ 2 \left(k _ {i} \left(x - x _ {0}\right) - \omega_ {i} \left(t - t _ {0}\right) + \varepsilon_ {i}\right) \right], (A6) \\ \end{array}
$$

where $g$ is the gravitational acceleration, $h$ is the water depth, $k$ is the wave number given by $k _ { i } = \omega _ { i } { } ^ { 2 } / g$ tanh(  k h) and the frequency is given by $\omega _ { i } = 2 \pi f _ { i }$ . The phase angle $\varepsilon$ is set to zero for the calculations in this work. The coefficients are given as follows: 

$$
\begin{array}{l} A ^ {+} = - \frac {\omega_ {i} \omega_ {j} (\omega_ {i} + \omega_ {j})}{D ^ {+}} \left(1 - \frac {1}{\tanh (k _ {i} h) \tanh (k _ {j} h)}\right) + \frac {1}{2 D ^ {+}} \left(\frac {\omega_ {i} ^ {3}}{\sinh^ {2} (k _ {i} h)}\right) \\ \left. + \frac {\omega_ {j} ^ {3}}{\sinh^ {2} (k _ {j} h)}\right), \\ \end{array}
$$

$$
\begin{array}{l} A ^ {-} = \frac {\omega_ {i} \omega_ {j} (\omega_ {i} - \omega_ {j})}{D ^ {-}} \left(1 + \frac {1}{\tanh (k _ {i} h) \tanh (k _ {j} h)}\right) + \frac {1}{2 D ^ {-}} \left(\frac {\omega_ {i} ^ {3}}{\sinh^ {2} (k _ {i} h)} \right. \\ \left. - \frac {\omega_ {j} ^ {3}}{\sinh^ {2} (k _ {j} h)}\right), \\ \end{array}
$$

$$
\begin{array}{l} B ^ {+} = \frac {\omega_ {i} ^ {2} + \omega_ {j} ^ {2}}{2 g} - \frac {\omega_ {i} \omega_ {j}}{2 g} \left(1 - \frac {1}{\tanh  (k _ {i} h) \tanh  (k _ {j} h)}\right) \\ \left(\frac {(\omega_ {i} + \omega_ {j}) ^ {2} + g (k _ {i} + k _ {j}) \tanh [ (k _ {i} + k _ {j}) h ]}{D ^ {+}}\right) + \frac {\omega_ {i} + \omega_ {j}}{2 g D ^ {+}} \left(\frac {\omega_ {i} ^ {3}}{\sinh^ {2} (k _ {i} h)} \right. \\ \left. + \frac {\omega_ {j} ^ {3}}{\sinh^ {2} (k _ {j} h)}\right), \\ \end{array}
$$

and 

$$
\begin{array}{l} B ^ {-} = \frac {\omega_ {i} ^ {2} + \omega_ {j} ^ {2}}{2 g} + \frac {\omega_ {i} \omega_ {j}}{2 g} \left(1 + \frac {1}{\tanh  (k _ {i} h) \tanh  (k _ {j} h)}\right) \\ \left(\frac {\left(\omega_ {i} - \omega_ {j}\right) ^ {2} + g \left(k _ {i} - k _ {j}\right) \tanh  \left[ \left(k _ {i} - k _ {j}\right) h \right]}{D ^ {-}}\right) + \frac {\omega_ {i} - \omega_ {j}}{2 g D ^ {-}} \left(\frac {\omega_ {i} ^ {3}}{\sinh^ {2} \left(k _ {i} h\right)} \right. \\ \left. - \frac {\omega_ {j} ^ {3}}{\sinh^ {2} (k _ {j} h)}\right), \\ \end{array}
$$

where 

$$
D ^ {\pm} = \left(\omega_ {i} \pm \omega_ {j}\right) ^ {2} - g \left(k _ {i} \pm k _ {j}\right) \tanh  \left[ \left(k _ {i} \pm k _ {j}\right) h \right].
$$

# References



Adcock, T.A.A., Taylor, P.H., 2009. Estimating ocean wave directional spreading from an Eulerian surface elevation time history. Proc. R. Soc. A: Math., Phys. Eng. Sci. 465 (2111), 3361–3381. 





Adcock, T.A.A., Taylor, P.H., Yan, S., Ma, Q.W., Janssen, P.A.E.M., 2011. Did the Draupner wave occur in a crossing sea? Proc. R. Soc. A: Math. Phys. Eng. Sci. 467 (2134), 3004–3021. 





Boccotti, P., Arena, F., Fiamma, V., Barbaro, G., 2012. Field experiment on random Wave forces acting on vertical cylinders. Probabilistic Engineering Mechanics. Elsevier Ltd, Langford Lane, Kidlington, Oxford, 39–51. 





Boo, S.Y., 2002. Linear and nonlinear irregular waves and forces in a numerical wave tank. Ocean Eng. 29 (5), 475–493. 





Boo, S.Y., 2006. Measurements of higher harmonic wave forces on a vertical truncated circular cylinder. Ocean Eng. 33 (2), 219–233. 





Chaplin, J.R., Rainey, R.C.T., Yemm, R.W., 1997. Ringing of a vertical cylinder in waves. J. Fluid Mech. 350, 119–147. 





Chau, F.P., Eatock Taylor, R., 1992. Second order wave diffraction by a vertical cylinder. J. Fluid Mech. 240, 571–599. 





Chen, L.F., Zang, J., Hillis, A.J., Morgen, G.C.J., Plummer, A.R., 2014. Numerical investigation of wave-structure interaction using OpenFOAM. Ocean Eng. 88, 91–109. 





CresT (Cooperative research on extreme seas and their impact) duration: 2008–2011. Joint Industry project (MARIN project 21003). 





Dalzell, J.F., 1999. A note on finite depth second-order wave-wave interactions. Appl. Ocean Res. 21, 105–111. 





Dixon, A.G., Greated, C.A., Salter, S.H., 1979. Wave forces on partially sub-merged cylinders. J. Waterw. Port. Coast. Ocean Div. 105 (4), 421–438. 





Faltinsen O., 1993. Sea loads on ship and offshore structures, Cambridge Ocean Technology Series, 〈isbn: 9780521458702〉. 





Fuhrman, D.R., Madsen, P.A., Bingham, H.B., 2006. Numerical simulation of lowestorder short-crested wave instabilities. J. Fluid Mech. 563, 415–441. 





Gerald C., Morgan J., Zang, J., 2010. Using the RasInterFoam CFD model for non-linear wave interaction with a cylinder. In: Proceedings of 20th Internation Socity of Offshore and Polar Engineers (ISOPE) Conference. Beijing, China. 





Higuera, P., Lara, J.L., Losada, I.J., 2013. Realistic wave generation and active wave absorbtion for Navier-Stokes models. Coast. Eng. 71, 102–118. 





Higuera, P., Lara, J.L., Losada, I.J., 2014. Three-dimensional interaction of waves and porous coastal structures using OpenFOAM®. Part 2: Application. Coast. Eng. 88, 259–270. 





Hu Z.Z., Greaves D., Raby A., 2014. Simulation of Extreme Free Surface Waves Using OpenFOAM. In: Proceedings of the 5th Conference on the Application of Physical Modelling to Port and Coastal Protection. Varna, Bulgaria. Vol. 2, pp. 243–252. 





Hu Z.Z., Causon D.M., Mingham C.M., Qian L., 2009. Numerical wave tank study of a wave energy converter in heave. In: Proceedings 19th Internation Socity of Offshore and Polar Engineers (ISOPE) Conference. Osaka, Japan. 383-388. 





Hu Z.Z., Causon D.M., Mingham C.M., Qian L., 2010. Numerical simulation of nonlinear wave interactions with a wave energy converter. In: Proceedings of 20th Internation Socity of Offshore and Polar Engineers (ISOPE) Conference. Beijing, China. 871-876 in Wave Energy Issue. 





Hu, Z.Z., Causon, D.M., Mingham, C.M., Qian, L., 2011. Numerical simulation of floating bodies in extreme free surface waves. J. Nat. Hazards Earth Syst. Sci. 11 (2), 519–527. http://dx.doi.org/10.5194/nhess-11-519-2011. 





Jacobsen, N.G., Fuhrman, D.R., Fredsøe, J., 2012. A wave generation toolbox for the opensource CFD library: OpenFoam. Int. J. Numer. Methods Fluid 70 (9), 1073–1088. 





Jasak, H., 1996. Error analysis and estimation for the finite volume method with applications to fluid flows (PhD thesis). Imperial College of Science, Technology and Medicine. 





Jiang C.B., Yao Y., Deng Y., and Deng B. 2015. Numerical Investigation of Solitary Wave Interaction with a Row of Vertical Slotted Piles. Journal of Coastal Research In-Press. 





Kriebel, D.L., 1998. Nonlinear’wave interaction with a vertical circular cylinder: Wave forces. Ocean Eng. 25 (7), 597–605. 





Lara, J.L., Ruju, A., Losada, I.J., 2011. Reynolds averaged Navier-Stokes modelling of long waves induced by a transient wave group on a beach. Proc. R. Soc. A467, 1215–1242. 





Li, Y.L., Zhu, R.C., MiAO, G.P., Fan, J., 2012a. Simulation of tank sloshing based on OpenFoam and coupling with ship motions in time domain. J. Hydrodyn. 24 (3), 450–457. 





Li, Y.L., Zhu, R.C., Miao, G.P., Fan, J., 2012b. Simulation of tank sloshing based on OpenFOAM and coupling with ship motions in time domain. J. Hydrodyn. Ser. B 24 (3), 450–457. 





Ma, Q.W., Wu, G.X., Eatock Taylor, R., 2001a. Finite element simulation of fully nonlinear interaction between vertical cylinders and steep waves. Part 1: methodology and numerical procedure. Int. J. Numer. Methods Fluids 36 (3), 265–285. 





Ma, Q.W., Wu, G.X., Eatock Taylor, R., 2001b. Finite element simulations of fully nonlinear interaction between vertical cylinders and steep waves. Part 2: numerical results and validation. Int. J. Numer. Methods Fluids 36 (3), 287–308. 





MacCamy, R.C., Fuchs, R.A., 1954. Wave forces on piles: a diffraction theory. Beach Erosion Board, Tech. Mem., 69. U.S. Army Corps of Engineers. 





Madsen, P.A., Sørensen, O.R., 1992. A new form of the Boussinesq equations with improved linear dispersion characteristics. Part 2. A slowly-varying bathymetry. Coast. Eng. 18, 183–204. 





Mai T., Hu Z.Z., Greaves D., Raby A. 2015. Investigation of hydroelasticity:wave impact on a truncated vertical wall. In: 25th Internation Socity of Offshore and Polar Engineers (ISOPE) conference. Hawaii, USA. pp. 647–654. 〈isbn:978-1-880653-89- 0〉. 





Malenica, S., Molin, B., 1995. Third harmonic wave diffraction by a vertical cylinder. J. Fluid Mech. 302, 203–229. 





Mayer, S., Garapon, A., Sørensen, L.S., 1998. A fractional step method for unsteady freesurface flow with applications to non-linear wave dynamics. Int. J. Numer. Methods Fluids 28 (2), 293–315. 





Morison, J.R., O’Brien, M.P., Johnson, J.W., Shaaf, S.A., 1950. The force exerted by surface waves on piles. Petrol. Trans. Am. Inst. Min. Pet. Eng. 189, 149–154. 





Newman, J.N., 1996. The second-order wave force on a vertical cylinder. J. Fluid Mech. 320, 417–443. 





Ning, D.Z., Teng, B., Eatock Taylor, R., Zang, J., 2008. Numerical simulation of nonlinear regular and focused waves in an infinite water-depth. J. Ocean Eng. 35 (8–9), 887–899. 





Ning, D.Z., Zang, J., Liu, S.X., Eatock Taylor, R.E., Teng, B., Taylor, P.H., 2009. Free surface and wave kinematics for nonlinear focused wave groups. J. Ocean Eng. 36 (15–16), 1226–1243. 





Omidvar, P., Stansby, P.K., Rogers, B.D., 2013. SPH for 3D floating bodies using variable mass particle distribution. Int. J. Numer. Methods Fluids 72 (4), 427–452. http://dx.doi.org/10.1002/fld.3749. 





Orszaghova, J., Taylor, P.H., Bothwick, A.G.L., Raby, A., 2014. Importance of secondorder wave generation for focused wave group run-up and overtopping. J. Coast. Eng. 93, 63–79. 





Rahman, M., Heaps, H.S., 1983. Wave forces on offshore structures: nonlinear wave diffraction by large cylinders. J. Phys. Oceanogr. 13 (12), 2225–2235. 





Roache, P.J., Ghia, K., White, F., 1986. Editorial policy statement on the control of numerical accuracy. ASME J. Fluids Eng. 108 (1), 2. 





Rusche, H., 2002. Computational fluid dynamics of dispersed two-phase flows at high phase fractions (PhD thesis). Imperial College, London. 





SAFE-FLOW. SAFE-FLLOating offshore structures under impact loading of shipped 





green water and waves. duration: 2000 -2003, EU contract No: G3RD-CT-2000- 00271. 





Schäffer, H.A., Klopman, G., 2000. Review of multidirectional active wave absorption methods. J. Waterw. Port. Coast. Ocean Eng., 88–97, (March/April). 





Taylor, R.E., Hung, S.M., 1987. Second order diffraction forces on a vertical cylinder in regular waves. Appl. Ocean Res. 9 (1), 19–30. 





Tromans P.S., Anaturk A.R., Hagemeijer P., 1991. A new model for the kinematics of large ocean waves-application as a design wave. In: Proceedings 1st Inter, Offshore and Polar Engineering Conference. Edinburgh, U.K. 64-71. 





Versteeg H.K., Malalaskera W.,1995. An introduction to computational fluid dynamicsthe finite volume method. 〈isbn:0-582-21884-5〉. 





Walker, dA.G., Taylor, P.H., Eatock Taylor, R.,, 2004. The shape of large surface waves on the open sea and the Draupner New Year. Appl. Ocean Res. 26 (34), 73–83. 





Weller H.G., 2002. Derivation, modelling and solution of the conditionally averaged twophase flow equations. Technical report TR/HGW/02. Nabla Ltd. 





Westphalen, J., Greaves, D.M., Williams, C.J.K., Hunt-Raby, A.C., Zang, J., 2012. Focused waves and wave-structure interaction in a numerical wave tank. Ocean Eng. 45, 9–21. 





Westphalen, J., Greaves, D.M., Raby, A., Hu, Z.Z., Causon, D.M., Mingham, C.G., Omidvar, P., Stansby, P.K., Rogers, B.D., 2014. Investigation of wave-structure interaction using state of the art CFD techniques. Open J. Fluid Dyn. 4 (1), 18–43. http://dx.doi.org/10.4236/ojfd.2014.41003. 





Wu G.X., Hu Z.Z., 2004. Simulation of nonlinear interactions between waves and floating bodies through a finite element based numerical tank. Proceedings of the Royal Society A, 460, 2797-2817. 





Zhang, Y., Zou, Q., Greaves, D.M., 2009. Numerical simulation of two phase flow using the level set method with global mass correction. Int. J. Numer. Methods Fluids 63 (6), 651–680. http://dx.doi.org/10.1002/fld.2090 (10). 





Zhang, Y.L., Zou, Q.P., Greaves, D., 2012. Air–water two-phase flow modelling of hydrodynamic performance of an oscillating water column device using a two-phase flow model. Renew. Energy 41, 159–170. http://dx.doi.org/10.1016/j.renene.2011.10.011). 

