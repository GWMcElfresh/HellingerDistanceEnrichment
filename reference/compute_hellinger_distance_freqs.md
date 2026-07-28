# Hellinger distance between two frequency vectors on a simplex.

Part of Paul Edlefsen's original enrichment procedure (subject-level
Hellinger geometry on softened multinomial compositions).

## Usage

``` r
compute_hellinger_distance_freqs(freqs1, freqs2)
```

## Arguments

- freqs1:

  Numeric frequency vector (sums to 1).

- freqs2:

  Numeric frequency vector (sums to 1).

## Value

Scalar Hellinger distance in \eqn\[0, 1\].
