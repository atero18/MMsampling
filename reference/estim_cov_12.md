# Estimation of the covariance between the potential outcomes of two modes of response under homoscedasticity assumptions and linearity

We suppose the homoscedasticity within the first and second and between
the modes of response potential outcomes. We will have an unbiased
estimator if the non-informativeness, unconfoundedness and conditional
mutual independence assumptions are verified ; the measure effect y_1k -
y_2k follow a linear model with x_k and for some units both potential
outcomes are known (which can't be under a sequential protocol). We must
have at least p + 1 (p the number of covariates) units such that their
measure effect are known or well estimated. The standard deviations
`sd1` and `sd2` must be unbiased as well.

## Usage

``` r
estim_cov_12(
  Y1,
  Y2,
  sd1,
  sd2,
  X,
  Yobs = NULL,
  modes = NULL,
  clamp = FALSE,
  warnClamp = TRUE,
  ...
)
```

## Arguments

- Y1:

  vector of the first mode outcomes (numeric vector of size N the size
  of the population).

- Y2:

  vector of the second mode outcomes (numeric vector of size N).

- sd1:

  true or estimated value of the standard deviation of the first mode of
  response potential outcomes (positive scalar).

- sd2:

  true or estimated value of the standard deviation of the second mode
  of response potential outcomes (positive scalar).

- X:

  design matrix (numeric matrix with N rows and p columns).

- Yobs:

  vector of the observed outcomes. Optional. Useful if `Y1` and `Y2` are
  not given. In that case the counterfactuals of the respondents are
  estimated with the function `estim_counterfactuals`. the value is not
  considered and therefore can be equal to NA (numeric vector of size N
  the size of the population).

- modes:

  vector of the selected mode of each unit. Optional. Used if the
  counterfactuals must be estimated (character vector or factor of size
  N).

- clamp:

  TRUE if the estimation of the covariance must be clamped if its
  absolute value is superior to `sd1` \* `sd2` (boolean).

- warnClamp:

  TRUE if a warning must be sent when a clamp is made (boolean).

- ...:

  arguments for the function
  [`MatchIt::matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).

## Value

an estimator of the covariance between the m1 and m2 potential outcomes,
unbiased under the assumptions, if `sd1`^2 and `sd2`^2 are unbiased and
the counterfactuals are the true values (scalar).
