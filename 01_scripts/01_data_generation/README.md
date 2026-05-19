# Alternative specifications of the simulation runner

The current `12_simulation_runner.R` is set up to reproduce the final thesis scenarios. The same runner can also be used to define other simulation conditions. This README gives the main substantive switches: how to vary the data-generating mechanism and how to vary the models fitted to each generated data set.

This is not a full user manual. It does not document logistical settings such as the number of replications, bootstrap resamples, cores, output paths, thread controls, or tuning-grid details. Those settings can be changed in the runner.

## Basic structure

Each simulation condition is defined by two choices:

1. **The data-generating mechanism**: how the panel data are generated, how strong the confounding is, whether the confounding is linear or nonlinear, and whether the confounding pattern changes over time.
2. **The fitted model set**: which longitudinal model is fitted and what information about the confounders the analyst is allowed to use.

In the runner, these choices are combined inside `scenario_defs` and then crossed with the selected sample sizes.

## The data-generating mechanism

Data are generated from a cross-lagged causal model for repeated `X` and `Y` measurements. The lagged dynamics are controlled by `Phi`, the contemporaneous variance-covariance structure of `X` and `Y` is controlled by `Sigma`, and the direct contribution of baseline confounders is controlled by a time-specific coefficient matrix, `Delta_t`.

The main DGM settings in the runner are:

```r
k        <- 5      # number of baseline confounders
T_waves  <- 5      # number of observed waves
burn_in  <- 20     # generated waves discarded before analysis

R2_total <- 0.15   # direct confounder contribution before any step
R2_nl    <- 0.05   # absolute nonlinear contribution, when interactions are used
R2_new   <- 0.35   # direct confounder contribution after a step
step_at  <- 3L     # first observed wave after which the step applies
```

The other key DGM objects are:

```r
Omega11  # covariance matrix of the baseline confounders C1, ..., Ck
Phi      # autoregressive and cross-lagged effects between X and Y
Sigma    # variance-covariance matrix of X_t and Y_t innovations
```

### Choosing the initial confounding pattern

The function `sample_delta_t()` samples the initial direct confounder effects. The important distinction is whether the confounding is linear only or whether it also contains interactions among the baseline confounders.

For **linear confounding only**, use the main effects of the baseline confounders:

```r
delta_linear <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = 0,
  include_2way = FALSE,
  include_3way = FALSE
)
```

For **nonlinear confounding**, include two-way and/or three-way interaction terms. In the thesis runner, the nonlinear contribution is set through `R2_nl / R2_total`:

```r
delta_nonlinear <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = R2_nl / R2_total,
  include_2way = TRUE,
  include_3way = TRUE
)
```

Here, `R2_total` controls the overall direct contribution of the confounder term at the scaling point. `rho_int` controls what share of that direct contribution is assigned to the interaction block.

## Three ways to generate confounding over time

After the initial `Delta_t` has been sampled, the runner builds a full `Delta_list`: one confounder-effect matrix for each generated wave. The helper scripts support three main time patterns.

### 1. Constant confounding effects

Use `generate_Delta_constant()` when the confounding structure should stay the same across all observed waves.

```r
Delta_list <- generate_Delta_constant(
  Delta_initial = delta_linear$Delta,
  n_waves = T_waves,
  burn_in = burn_in
)
```

This can be used with either `delta_linear` or `delta_nonlinear`. With `delta_linear`, the DGM has constant linear confounding. With `delta_nonlinear`, it has constant nonlinear confounding.

### 2. Stepwise, rank-stable confounding effects

Use `generate_Delta_stepwise()` when the strength of confounding should change at a chosen wave, but the coefficient pattern should remain the same.

```r
Delta_list <- generate_Delta_stepwise(
  Delta_initial = delta_nonlinear$Delta,
  n_waves = T_waves,
  burn_in = burn_in,
  step_at = step_at,
  R2_old = R2_total,
  R2_new = R2_new
)
```

This keeps the signs and relative ordering of the coefficients stable. From `step_at` onward, all coefficients are rescaled so that the direct confounder contribution changes from `R2_old` to `R2_new`.

### 3. Stepwise mixture confounding effects

Use `generate_Delta_stepwise_mixture()` when the confounding structure itself should change after the step.

```r
Delta_list <- generate_Delta_stepwise_mixture(
  Delta_initial = delta_nonlinear$Delta,
  n_waves = T_waves,
  Omega = delta_nonlinear$Omega,
  burn_in = burn_in,
  step_at = step_at,
  R2_old = R2_total,
  R2_new = R2_new,
  lambda_L = 0.50,
  lambda_int = 0.50,
  seed = 123
)
```

This creates a post-step coefficient pattern by mixing the original coefficients with a new random direction. `lambda_L` controls how much the main-effect coefficients are allowed to change. `lambda_int` does the same for the interaction coefficients. Values closer to `0` keep the post-step pattern closer to the original pattern; values closer to `1` allow more rank-order change.

## Defining scenarios in the runner

The simulation runner defines conditions through `scenario_defs`. A scenario mainly specifies:

- a scenario id and label;
- the analyst-side confounder information;
- the function used to generate `Delta_list`.

A minimal scenario definition looks like this:

```r
list(
  id = 1L,
  label = "constant_linear",
  analyst_order = 1L,
  exclude_general = NULL,
  make_Delta = function() {
    generate_Delta_constant(
      Delta_initial = delta_linear$Delta,
      n_waves = T_waves,
      burn_in = burn_in
    )
  }
)
```

The thesis runner uses four scenarios:

1. constant linear confounding;
2. constant nonlinear confounding;
3. stepwise nonlinear confounding with stable coefficient ranks;
4. the same stepwise nonlinear DGM as scenario 3, but with `C1` and `C2` omitted from every analyst-side adjustment.

Additional scenarios can be added by adding new entries to `scenario_defs`.

## Varying the fitted models

The fitted models are built in `make_model_set()`. Each model specification combines two choices:

1. the longitudinal SEM to fit;
2. the adjustment strategy used before or inside that SEM.

### Longitudinal model type

The main model choice is controlled by `sem_model`:

```r
sem_model = "clpm"    # cross-lagged panel model
sem_model = "riclpm"  # random-intercept CLPM
sem_model = "dpm"     # dynamic panel model
```

The thesis runner fits all three model types in each scenario. Results are reported only for the RI-CLPM and CLPM. 

### Adjustment strategy

Adjustment is controlled by `residualizer` and by the analyst-side confounder order.

```r
residualizer = "none"    # no residualisation; SEM may include observed confounders directly
residualizer = "linear"  # residualise X and Y using a linear model
residualizer = "xgb"     # residualise X and Y using XGBoost
residualizer = "enet"    # residualise X and Y using elastic net
```

When `residualizer = "none"`, confounders can be included directly in the SEM through `sem_c_order`. When `residualizer` is not `"none"`, the SEM is fit to residualised `X` and `Y`, and the confounders used during residualisation are controlled by `residualizer_c_order`.

The confounder order has the same meaning in both places:

```r
0L  # no observed confounders
1L  # baseline confounder main effects: C1, ..., Ck
2L  # main effects plus all two-way interactions
3L  # main effects plus all two-way and three-way interactions
```

Specific confounder terms can be removed with `sem_exclude` or `residualizer_exclude`. In the runner, this is passed through `exclude_general`, for example:

```r
omit_vars <- c("C1", "C2")
```

This is useful for scenarios where the DGM contains confounding terms that are not fully observed by the analyst.

## Current thesis model set

The current runner fits nine models in every scenario:

1. CLPM without adjustment;
2. RI-CLPM without adjustment;
3. DPM without adjustment;
4. CLPM with observed linear confounder control;
5. RI-CLPM with observed linear confounder control;
6. DPM with observed linear confounder control;
7. CLPM after XGB residualisation;
8. RI-CLPM after XGB residualisation;
9. DPM after XGB residualisation.

Models 1, 3, 6, and 9 are dropped from the output figures and tables. To add models, edit `make_model_set()` and add further `make_model_spec()` entries. The most important fields to change are usually `name`, `sem_model`, `residualizer`, `sem_c_order`, `sem_exclude`, `residualizer_c_order`, and `residualizer_exclude`.


