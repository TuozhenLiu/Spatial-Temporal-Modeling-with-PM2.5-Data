

This is the R code for the paper "Spatially Dependent Autoregression Model".

To use the code, please set the working path to be "../simu/".


### 1. Main folder
 - simulator.R: functions which simulate the dataset
 - est_infer_func.R: functions which estimate the model and make inference
 - predictFunc2.R: functions which conduct predictive kriging
 
### 2. simu folder

Please run the following R files to obtain estimation results:
 - simu_main_exp.R: Simulation results for exponential model with 1,000 replications. The results are summarized in Table 1.
 - simu_main_quad.R: Simulation results for rational quadratic model with 1,000 replications. The results are summarized in Table 2.
 - simu_main_mixed.R: Simulation results for mixed model with 1000 replications. The results are summarized in Table 3.

