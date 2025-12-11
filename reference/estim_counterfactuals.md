# Estimation of the counterfactuals of each respondent

In this version considers only the case with two modes. Could be
generalised. Requires the package
[MatchIt](https://cran.r-project.org/web/packages/MatchIt/).

## Usage

``` r
estim_counterfactuals(Yobs, modes, X, ...)
```

## Arguments

- Yobs:

  vector of the observed outcomes. (numeric vector of size N the size of
  the population).

- modes:

  vector of the selected mode of each unit. The first mode of the
  protocol is defined as "m1", the second as "m2" (character vector or
  factor of size N).

- X:

  design matrix. The rows corresponding to the non-respondents can
  contain NA (numeric matrix with N rows and p columns).

- ...:

  arguments for the function
  [`MatchIt::matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).

## Value

for each mode a vector of size N containing for each respondent its
outcome or its estimated counterfactual depending on if it answered by
the mode or the other. Equal to NA for the non-respondents (numeric
matrix of dimension (N,2)).

## See also

[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
