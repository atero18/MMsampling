# Estimation of t0

Estimation of the total \\t_0\\ under an almost sure equality between
the potential outcomes of each unit. The total can be estimated using
the mode selection probabilities or the global answer probabilities.

## Usage

``` r
estim_tot_0_no_p1(Yobs, modes, pi, p2)

estim_tot_a(Yobs, responses, pi, pa)
```

## Arguments

- Yobs:

  vector of the observed outcomes. (numeric vector of size N the size of
  the population).

- modes:

  vector of the selected mode of each unit. The first mode of the
  protocol is defined as "m1", the second as "m2" (character vector or
  factor of size N).

- pi:

  vector containing the first-order inclusion probabilities (numeric
  vector of size N).

- p2:

  vector containing the true or estimated selection probabilities of the
  second mode (numeric vector of size N).

- responses:

  vector indicating which unit answered. If a unit answered then the
  corresponding value is "r" (character vector or factor of size N).

- pa:

  vector containing the true or estimated global answer probability
  (numeric vector of size N).

## Value

the estimation of the total \\t_0\\ (scalar).

## Functions

- `estim_tot_0_no_p1()`: Sum an estimator per mode. Do not use the first
  mode selection probabilities.

- `estim_tot_a()`: Use the probability to answer by any mode instead of
  by each mode.

## See also

Other total estimators:
[`estim_tot_m1()`](https://atero18.github.io/MMsampling/reference/estim_tot_m1.md),
[`estim_tot_m2()`](https://atero18.github.io/MMsampling/reference/estim_tot_m2.md),
[`estim_tot_phi_y()`](https://atero18.github.io/MMsampling/reference/estim_tot_phi_y.md)
