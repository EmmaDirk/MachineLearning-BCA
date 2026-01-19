In this folder are scripts that are outdated, which were used in earlier stages in the project. They are kept in case useful parts can be integrated at a later stage. Currently, the following scripts are available:

In `01_data_simulation`: 
* `01_equal_effects_over_time`: contains the first versions of data simulation scripts for simulating data with equal effects of confounders over time. Data simulation is now done another way, but the ways of including non-linearity may be useful. 
* `02_varying_effects_over_time`: contains the first versions of data simulation scripts for simulating data with varying effects of confounders over time. Data simulation is now done another way, but the ways of including non-linearity may be useful.

In `02_research_report`:

* `00_models_T=5`: some lavaan synatx for models with T=5 measurement occasions. Deprecated since models are now built using string builders that adapt to T. 
* `01_results_v1`: plotting of results from simulation studies (v1). Now integrated in the main simulation engine.
* `02_results_v2`: plotting of results from simulation study (v2). Now integrated in the main simulation engine.
* `031_sim_study_1`: earlier version of simulation study 1. Now integrated in the main simulation engine.
* `041_sim_study_2`: earlier version of simulation study 2. Now integrated in the main simulation engine.
* `032_sim_study_1_no_sampling`: earlier version of simulation study 1 without sampling each beta-matrix for a run each time. This was an improved verison of 031, which now is included in the main simulation engine.
* `042_sim_study_2_no_sampling`: earlier version of simulation study 2 without sampling each beta-matrix for a run each time. This was an improved verison of 041, which now is included in the main simulation engine.