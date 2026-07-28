# Hellinger enrichment on synthetic compositions

Understanding whether categorical compositions differ between
experimental groups is critical when cell-level clustering is summarized
to subject-level frequencies. However, treating those frequencies as
fixed points in a simplex ignores sampling variation in the estimated
multinomial distribution. Effect size in this package is the mean
between-group pairwise Hellinger distance divided by the mean
within-group pairwise distance. Here, we walk through the
`HellingerDistanceEnrichment` workflow on entirely synthetic data so
that the inference path remains visible without cohort-specific noise.

## Build a synthetic long counts table

We first construct a long count table with subject identifiers, category
labels, group membership, and counts. A modest enrichment is then
planted in the Treatment group for a single category. That planted
structure is the known signal against which the permutation and Bayes
procedures can be compared.

``` r

library(HellingerDistanceEnrichment)

set.seed(2026)
long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
  n_subjects_per_group = 6,
  n_categories = 5,
  groups = c("Control", "Treatment"),
  seed = 2026
)

# Plant a modest enrichment in Treatment for one category.
long_table <- HellingerDistanceEnrichment:::plant_group_structure(
  long_table,
  target_group = "Treatment",
  target_category = "Category0",
  boost = 4
)

head(long_table)
#>     subjectId  category   group  n
#> 1 Control_S01 Category0 Control 22
#> 2 Control_S01 Category1 Control 20
#> 3 Control_S01 Category2 Control 15
#> 4 Control_S01 Category3 Control 20
#> 5 Control_S01 Category4 Control 23
#> 6 Control_S02 Category0 Control 20
```

## Extract and compare compositions

[`ExtractClusterComposition()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/ExtractClusterComposition.md)
aggregates the long table into a subject-by-category count matrix with a
group factor.
[`CompareGroupCompositions()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/CompareGroupCompositions.md)
then returns an omnibus contrast together with pairwise group contrasts.
The permutation method shuffles group labels while holding the estimated
compositions fixed, which tests whether the observed between-to-within
distance ratio is large relative to a label-exchange null.

``` r

composition <- ExtractClusterComposition(long_table)

perm_result <- CompareGroupCompositions(
  composition,
  method = "permutation",
  nPermutations = 200,
  seed = 42
)

perm_result
#> HellingerEnrichmentResult (permutation): omnibus effectSize=2.1878, pValue=0.0050
#>   1 contrast(s) tested
```

The Bayes path samples subject-level compositions from a conjugate
Dirichlet posterior and, at each draw, runs a full permutation null. In
this vignette the posterior and permutation counts are kept small for
runtime; production analyses should increase both.

``` r

bayes_result <- CompareGroupCompositions(
  composition,
  method = "bayes",
  nPosterior = 20,
  nPermutations = 50,
  seed = 42
)

bayes_result$omnibus
#> $effectSize
#> [1] 1.843022
#> 
#> $pValue
#> [1] 0.02156863
#> 
#> $effectCiLow
#> [1] 1.554202
#> 
#> $effectCiHigh
#> [1] 2.152769
```

## Plot contrasts

The permutation panel marks significance against the label-shuffle null.
The Bayes panel shows posterior-mean effect sizes with 95% intervals, so
composition uncertainty and group-label uncertainty appear on the same
display.

``` r

PlotCompositionContrasts(perm_result)
```

![Permutation contrast panel with omnibus
reference.](hellinger-enrichment_files/figure-html/plot-perm-1.png)

Permutation contrast panel with omnibus reference.

``` r

PlotCompositionContrasts(bayes_result)
```

![Bayesian posterior-mean contrasts with 95%
intervals.](hellinger-enrichment_files/figure-html/plot-bayes-1.png)

Bayesian posterior-mean contrasts with 95% intervals.

The two methods therefore answer related but non-equivalent questions.
Plain permutation asks whether group labels are exchangeable given fixed
composition estimates. The posterior-predictive Bayes path nests
composition uncertainty inside that same label null, which the plain
permutation test does not account for. When category counts are sparse,
that nesting can change which contrasts appear stable. The synthetic
enrichment used here is strong enough that both paths recover the
planted structure; real cohorts will often require the Bayes path to
decide whether an apparent contrast survives sampling variation in the
compositions themselves.
