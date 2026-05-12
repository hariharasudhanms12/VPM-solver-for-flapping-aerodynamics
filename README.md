# VPM-solver-for-flapping-aerodynamics
This is a Vortex-Panel-Method solver coded in Python for capturing the wake characteristics and classifying if von Karman or reverse von Karman wake. However, the values of the lift and drag coefficients are different from the expectation although the graphs are similar. 

This is applicable only for very thin airfoils, as the vortices are placed along the chord of the airfoil. The airfoil can have both pitching and heaving motion (plunging motion).

Parameters to be specified:
1. Number of panels
2. Chord length
3. Time-step
4. Total time of simulation
5. Freestream velocity
