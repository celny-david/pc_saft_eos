# Info
Implementation of PC-SAFT equation of state in Fortran 95.
A part of research into Equation of states, specifically Perturbed Chain Statistical Associating Fluids Theory (for more info see [wiki](https://en.wikipedia.org/wiki/PC-SAFT)). This type of EoS was first introduced in 2001 by J. Gross & G. Sadowski in [this article](https://doi.org/10.1021/ie0003887).

# Capabilities
This implementation is using a state machine approach to effective calculation of residual properties from the given temperature, density input. This is intentional as this way no aditional iteration is required (such as in temperature, pressure input case). The state machine allows for recalculation of only the necessary contributions when the input values change making repated calls with similar inputs cheaper. Another important point is the implementation uses analytical approach to contribution differenciation (making it espetially complex for composition and higher order differentiation).

## fluid characteristics:
 - standard behavior accounted for with Hard sphere, Hard Chain and Dipersion contributions
 - Polarity including Diploar-Dipolar, Dipolar-Quadrupolar, Quadrupolar-Quadrupolar contributions
 	- beware that dipolar quadrupolar sunbstance testing was limited and results produced should be put under higher scrutiny
 - Association contribution including 2A,2B,3A,3B,3C,4A,4B,4C
 	- beware that assotiation is not fully tested and should be considered unreliable

## properties
calculated **residual** properties include:
 - pressure and its first and second density derivative
 - chemical potentials
 - fugacity coefficients
 - enthalpy
 - entropy
 - Gibbs free energy

## call options
 - simple executable accepting CLI arguments from console
 	- alternatively one can supply a file containing the arguments with example shown in `res/default.txt`
 - dynamic library with the example call in `test/test_fortran_library_dynamic.f90`
 	- this includes python interface example call in `test/test_FORTRAN_python_interface.py`
 	- also showcase MATLAB interface `test/test_FORTRAN_MATLAB_interface.m`
 		- beware that this was last tested for MATLAB 2017

## substances
Substances included (`res/parameter_list.txt`) in this showcase version are intentionally limited to basic alkanes, ax the substance fiting is part of ongoing research and is a coslty endeavour. It is still possible to fill in own parameters into the file to remedy this situation. Format is provided within the file itself.

## call diagram
For better understanding of the design a call diagram is provided showcasing the structure of the program. The diagram is available in graph format and png figure in `res/` folder. In same folder is also example for standard forms used for functions,subroutines and modules within the project to preserve the code structure.

# Conclusions
If you find any issue with the EoS plese do not hesitate to report it. Also be aware that this implementation was primarily developed on/for Linux. Windows is supported within the `CMakeLists.txt` but not regularly verified so be careful.
