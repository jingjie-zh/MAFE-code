# MAFE-code
This contains the code used in _Countercyclical Unemployment Benefits A General Equilibrium Analysis of Transition Dynamics_. , which will reproduce the tables and figures therein. 

- `parameters.m`: MATLAB file contains most parameters needed.
- `library`: folder contains 4 MATLAB files `uniformjobloss.m`,`StationaryE.m`,`HJB_transition.m` and `dynamics.m`, where are the 4 basic functions used to compute the stationary and time-dependent equilibrium in our model.
- `main.m`: MATLAB file to generate equilibrium results under different settings.
-  `pre.mat`, `post_cc.mat`, `post_ac.mat`: MAT-files contains results in the baseline model. 
- `pre_lesspatient.mat`, `post_cc_lesspatient.mat`, `post_ac_lesspatient.mat`: MAT-files contains results in the "less-patient" model where the time preference parameter `\rh0` is smaller than that in the baseline model.
-  `CompShockSize.mat` MAT-file contains the data required for studying the comparative static exercise in response to changes in shock size.
-  `table2.m`,  `table3.m`,  `table4.m`,  `figure1.m`, `figure2.m`, `figure3.m`:  MATLAB files to generate the corresponding tables and figures in the paper.

If you wish to replicate a specific table or figure, ensure that the required data is successfully loaded and then execute the corresponding MATLAB file.

For reproducing the equilibrium results, running 'main.m' will regenerate the MAT-files. Please note that this process might take a couple of hours.
