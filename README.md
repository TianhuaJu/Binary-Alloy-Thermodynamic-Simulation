# Binary-Alloy-Thermodynamic-Simulation

Using this simulation software:
1. Install Python3 and these non-standard packages: matplotlib, plotly, numpy
2. To run simulation with one set of parameters only, use main.py file

3. The most interesting results come from running functions in the simulationAnalysis.py file

Summary description:
This simulation models the arrangement of atoms in a binary alloy as a lattice gas, considering the laws of thermodynamics,
in both 2 and 3 dimensions. It swaps adjacent atoms if this move lowers the system energy or if there is sufficient
thermal energy for such an event. The atoms to be swapped are determined randomly and potential swaps are considered
many times (on the order of millions). The input parameters varied were system temperature, composition and alloy-matrix
bond energy (E_am) and, by quantifying the order of the system, one is able to examine the relationship between these
parameters and the atomic arrangement, within the confines of this model.
For more detail please see: Simulation Details.pdf


Explanation of simulation parameters:
>T = system temperature, in Kelvin, 300 K is a good starting point
>alloyFraction = percentage of atoms in the system that will be foreign/alloying atoms, integer values between
1 and 50 are recommended
>E_am = interaction energy (in eV) between unlike (host & alloying) atoms; it is recommended to use values
of -0.1, 0.0 or 0.1 eV
>gridLength = length of numpy array. E.g. for 2D, gridLength=5 will run simulation on a 5x5 array; values of
70 and 40 are recommended for 2D and 3D simulations respectively
nIterations = no. of iterations that simulation will run for, with a single set of parameters; on one iteration
a random pair of nearest neighbour atoms are swapped if it is energetically favourable to do so; It is
recommended to use the plotEnergyConvergence function in simulationAnalysis.py to determine how many iterations
are required before system energy stops changing by large amounts with more iterations, as the optimum no.
of iterations varies depending on gridLength & no. of dimensions; e.g. 250,000 is suitable for 70x70 grid in 2D
>dimensions = no. of dimensions to run simulation in; valid values are 2 or 3, for 2D or 3D