# Estimation of t2

Estimation of the total of the second mode of response potential
outcomes \\t_2\\.

## Usage

``` r
estim_tot_m2(Yobs, modes, pi, p1, p2, phi = rep(1, length(Yobs)))
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

- p1:

  vector containing the true or estimated selection probabilities of the
  first mode (numeric vector of size N).

- p2:

  vector containing the true or estimated selection probabilities of the
  second mode (numeric vector of size N).

- phi:

  optional vector containing weights for the total. The values must be
  between 0 and 1 (numeric vector size N).

## Value

the estimation of the total \\t_2\\, possibly weighted (scalar).

## See also

Other total estimators:
[`estim_t0`](https://atero18.github.io/MMsampling/reference/estim_t0.md),
[`estim_tot_m1()`](https://atero18.github.io/MMsampling/reference/estim_tot_m1.md),
[`estim_tot_phi_y()`](https://atero18.github.io/MMsampling/reference/estim_tot_phi_y.md)
