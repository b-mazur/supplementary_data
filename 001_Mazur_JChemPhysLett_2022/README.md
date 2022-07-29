Supporting information for ["Quasicontinuous Cooperative Adsorption Mechanism in Crystalline Nanoporous Materials"](https://doi.org/10.1021/acs.jpclett.2c01752), B. Mazur, F. Formalik, K. Roztocki, V. Bon, S. Kaskel, A. V. Neimark, L. Firlej, and B. Kuchta, **2022**, _J. Chem. Phys. Lett._, DOI: [10.1021/acs.jpclett.2c01752](https://doi.org/10.1021/acs.jpclett.2c01752)

**Contents**

- [simulation_input](simulation_input): RASPA2 input files with python files for running simulation / results analysis
  - [blocked_pore](simulation_input/blocked_pore): cif files with large (lp) and small (sp) pore in the middle of cell and simulation.input file for running simulation with moves restriction to pore in the middle
  - [density_maps](simulation_input/density_maps): simulation.input file for running simulation with computation of methane density, vtk2numpy.py file for plotting density maps from VTK files
  - [energy_surface](simulation_input/energy_surface): simulation.input file for running simulation with one MD move to calculate potential energy at given location and energy_map.py file for running simulations in loop to cover area
  - [GC-TMMC](simulation_input/GC-TMMC): simulation.input file to run GC-TMMC simulation with given number of molecules in the system and calcDOS.py file to posprocess data to obtain density of states **COMMENT:** to calculate DOS you have to run 0..Nmax simulations
  - [GCMC](simulation_input/GCMC): simulation.input file used to run standard GCMC simulations
