This program solves the 2D steady state heat diffusion problem provided boudany conditions and a single heat point source.

1. The coordinate of the heat source can be given in the function forcing_source defined at the end in geometry.cxx file.
The heat source value (watts/m^2) is also provided in the same function.
2. This program solves for two boundary conditions(upper and lower). Any integral values of flags can be choosen and the flags value and the actual temperature value at those boundaries is given in the function named boundary_value. 
3. This program used five different regions. These regions are defined in .poly file in examples folder. the mattype value is the conductivity value which can be directly given in .poly file and program will honour for those values. (conductivity value should be between 0 and 1.


