# Overview scripts for additional simulation results

These scripts provide a reusable way to inspect combined simulation results produced by the simulation runner. This README summarizes the main workflow: identify which models are present, select the scenario or conditions of interest, and generate a standard set of summary tables and plots.

The scripts are intended for browsing many different result files, not only the example shown in `03_example_usage.R`.

## Scripts

### `01_overview_models.R`

This script defines `overview_models()`, a compact helper for checking which fitted models are present in a combined results data frame.

Use it when you first load a results file and want to see, for example:

- which model families are included: CLPM, RI-CLPM, DPM, or their free-loading variants;
- whether models use no adjustment, SEM-based confounder adjustment, or residualization;
- which residualizer was used, such as linear model, elastic net, or BCA/XGBoost;
- which scenarios, sample sizes, and occasions are represented;
- how many rows and replications are available per model specification.

This is mainly an inspection step. It helps confirm that the results file contains the models and conditions you think it contains before making plots.

### `02_plot_overview_suite.R`

This is the main overview script. It defines `plot_overview_suite()`, which takes a combined results data frame and returns one omnibus object containing both tables and plots.

The returned object is useful for comparing fitted models across occasions, scenarios, and sample-size conditions. It includes summaries for:

- SEM performance: bias, relative bias, RMSE, mean standard error, detection probability, and standard-error calibration;
- residualizer diagnostics: prediction MSE and prediction \(R^2\), where available;
- true confounder-effect trajectories, if `true_delta_t_vector` is present;
- convergence and improper-solution diagnostics, based on the available flag columns;
- family-specific plots for CLPM, RI-CLPM, and DPM methods.

The function is deliberately broad. It can be used for a single scenario and sample size, or for a larger collection of scenarios and sample sizes. When multiple scenarios or sample sizes are retained, the summaries keep them separate and the plots facet where appropriate.

### `03_example_usage.R`

This script gives a minimal working example. It loads one combined `.rds` file, sources the overview functions, prints the model overview, creates an overview object for one scenario and one sample size, and then displays selected tables and plots.

Use this file as a template. For a different result file or condition, the main changes are usually the path to the `.rds` file and the filters supplied to `overview_models()` or `plot_overview_suite()`.




