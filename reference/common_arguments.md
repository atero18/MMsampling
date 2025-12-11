# Arguments commonly used

Arguments commonly used

## Arguments

- Yobs:

  vector of the observed outcomes. (numeric vector of size N the size of
  the population).

- modes:

  vector of the selected mode of each unit. The first mode of the
  protocol is defined as "m1", the second as "m2" (character vector or
  factor of size N).

- I:

  vector containing the membership indicators equal to TRUE if unit k
  has been selected (logical vector of size N).

- pi:

  vector containing the first-order inclusion probabilities (numeric
  vector of size N).

- p1:

  vector containing the true or estimated selection probabilities of the
  first mode (numeric vector of size N).

- pi_mat:

  matrix containing the second-order inclusion probabilities (symmetric
  numeric matrix of order N).

- p2:

  vector containing the true or estimated selection probabilities of the
  second mode (numeric vector of size N).

- phi:

  optional vector containing weights for the total. The values must be
  between 0 and 1 (numeric vector size N).
