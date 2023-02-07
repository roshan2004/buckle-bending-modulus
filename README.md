# Calculation of Bending Modulus (Kc) via Buckled Membrane
This repository provides a simple workflow for the calculation of the Bending Modulus (Kc) of a buckled membrane, from start to finish.

Requirements:
- GROMACS
- insane.py

The method to calculate bending modulus comes from:
> Hu, M., Diggins, P., & Deserno, M. (2013). Determining the bending modulus of a lipid membrane by simulating buckling. The Journal of chemical physics, 138(21), 214110. https://doi.org/10.1063/1.4808077

Tests on optimal system size and overall workflow/method development originated in:
> Eid, J., Razmazma, H., Jraij, A., Ebrahimi, A., & Monticelli, L. (2020). On Calculating the Bending Modulus of Lipid Bilayer Membranes from Buckling Simulations. The Journal of Physical Chemistry B, 124 (29), 6299-6311. https://doi.org/10.1021/acs.jpcb.0c04253 

### Simulation Explanations:
Two production runs are required for the calculation of Kc: 
- Unbuckled bilayer unrestrained in x, but with pressure coupling turned off in y, allowing the bilayer only to move in x. This allows for the calculation of the average length of x.
- Buckle bilayer restrained in both x and y, fixing the buckle shape. This allows for the calculation of the force the buckle exerts on the x side of the simulation box due to being fixed in a buckled state.

Further explanation of all simulation parameters may be found in the corresponding .mdp files, provided in /mdp.

Following the creation of a bilayer, which is at least 32nm x 8nm x 20nm (Eid, 2020), the system is energy minimised and equilibrated for 500ns. 

From this equilibration, (3.eq), the first production run (4.lx) is performed for 8 microseconds.

Also from this equilibration, the bilayer is deformed over 240ns in the x-dimension (5.deform), and a 2 microsecond equilibration (6.buckle_eq) is run on the resulting system.

The second production run (7.restrain) can then be conducted on the buckled equilibration, again for 8 microseconds.


### Simulation hierarchy
                            ┌────────────────┐
                            │                │
                            │     1.init     │
                            │                │
                            └───────┬────────┘
                                    │
                            ┌───────▼────────┐
                            │                │
                            │      2.em      │
                            │                │
                            └───────┬────────┘
                                    │
                            ┌───────▼────────┐
                            │                │                                             
                            │      3.eq      │
                            │                │
                        ┌───┴────────────────┴────┐
                        │                         │
                        │                         │
                        │                         │
               ┌────────▼───────┐        ┌────────▼───────┐
               │                │        │                │
               │      4.lx      │        │    5.deform    │
               │                │        │                │
               └────────────────┘        └────────┬───────┘
                                                  │
                                         ┌────────▼───────┐
                                         │                │
                                         │   6.buckle_eq  │
                                         │                │
                                         └────────┬───────┘
                                                  │
                                         ┌────────▼───────┐
                                         │                │
                                         │   7.restrain   │
                                         │                │
                                         └────────────────┘
