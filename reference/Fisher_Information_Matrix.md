# Fisher Information Matrix (FIM) for a logistic model and weights equal to 1

Fisher Information Matrix (FIM) for a logistic model and weights equal
to 1

## Usage

``` r
Fisher_Information_Matrix(prob, Z, maskSubset = !logical(nrow(Z)))
```

## Arguments

- prob:

  probability to answer. Can be true or estimated. Can be equal to NA
  for the units that are not in the subset defined by `maskSubset`
  (numeric vector of size N the size of the population).

- Z:

  design matrix. The rows corresponding to the units that are not in the
  subset can contain NA (numeric matrix with N rows and q columns).

- maskSubset:

  mask indicating which unit should be used in the calculation of the
  FIM. Default to the entire set (logical vector of size N).

## Value

the Fisher Information Matrix or its estimation (numeric matrix of order
q).
