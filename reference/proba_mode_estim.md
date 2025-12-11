# Estimation of mode choising probabilities

Estimation of mode choising probabilities

## Usage

``` r
estim_response_prob_global(
  I,
  modes,
  Z,
  RGH = NULL,
  constRGH = TRUE,
  chosenOnly = FALSE
)

estim_response_prob_sequential(
  I,
  Z,
  modes,
  orderModes,
  RGH = NULL,
  constRGH = FALSE,
  link = "logit",
  chosenOnly = TRUE
)
```

## Arguments

- orderModes:

  the order modes are considered. The first mode will have an absolute
  probability of selection, the second a probability a selection
  conditionally to the non-selection of the first one, etc.

## Functions

- `estim_response_prob_global()`: Estimation by a multimodal point of
  view

- `estim_response_prob_sequential()`: Estimation sequentially with
  conditional probabilities (interesting particularly with sequential
  mixed-mode protocols)
