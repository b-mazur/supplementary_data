Supporting information for "Application of Free Energy Landscape Interpolation and Extrapolation for Large-Scale Computational Screening of MOFs for Use in Water Adsorption", B. Mazur, L. Firlej, and B. Kuchta, **2024**

**Contents**

- [isobars](isobars): numerical data of adsorption isobars of TIP4P water model adsorption on MOF-303, MOF-LA2-1 and NU-1000 at 298 - 343 K at the pressure given by the last value in the file name. The value in the columns refers to the number of macrostates simulated directly.
- [isotherms](isotherms): numerical data of adsorption isotherms of TIP4P water model adsorption on MOF-303, MOF-LA2-1 and NU-1000 at 298 and 343 K.
  - [extrapolated](isotherms/extrapolated): numerical data of isotherms extrapolated from data collected at 298 K. The last value in the file name refers to the number of macrostates simulated directly.
  - [interpolated](isotherms/interpolated): numerical data of isotherms constructed using interpolation scheme. The last value in the file name refers to the number of macrostates simulated directly.
  - [reference](isotherms/reference): numerical data of isotherms calculated using full data and directly at given temperature.
- [lnp_files](lnp_files): numercial data of natural logarithm of macrostate probability $\ln\Pi$ at given macrostate calculated using WL/TMMC with [FEASST](https://doi.org/10.6028/jres.123.004).
- [prob_files](prob_files): numerical data of transtion probabilites ($P(N->N+1)$ and $P(N->N-1)$) at given macrostate calculated using *NVT + ghost swap* approach with [RASPA2](https://github.com/b-mazur/RASPA2_GC-TMMC.git).
- [simulation_input](simulation_input): RASPA2 and FEASST simulation input files.
  - [feasst](simulation_input/feasst): .fstprt files with information about TIP4P water molecule and frameworks (MOF-303, MOF-LA2-1 and NU-1000) and force field used, and .py example file to run WL/TMMC simulation using FEASST (for NU-1000 in this case).
  - [raspa](simulation_input/raspa): all RASPA2 input files needed to run simulations. Files with force field information are grouped in corresponding folders with the name of the framework. 