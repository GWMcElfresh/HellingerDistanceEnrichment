# Compare categorical compositions across groups.

Tests whether between-group Hellinger distances exceed within-group
distances using a permutation null or conjugate Dirichlet
posterior-predictive Bayes. Always reports an omnibus statistic. When
`contrasts` is omitted, also tests every pairwise group comparison; when
`contrasts` is supplied, those specifications are the contrast family
(auto pairwises are omitted unless `includePairwise = TRUE`).

## Usage

``` r
CompareGroupCompositions(
  x,
  method = c("permutation", "bayes"),
  contrasts = NULL,
  nPermutations = 1000,
  nPosterior = 100,
  nCores = 1,
  seed = NULL,
  priorPseudocounts = 1/2,
  pAdjustMethod = "holm",
  credibleInterval = c("hpdi", "cpi"),
  includePairwise = NULL,
  nope = "shuffle",
  energyNope = "shuffle",
  nNopeShuffles = 200,
  ...
)
```

## Arguments

- x:

  A CategoryComposition object or long counts table.

- method:

  Inference method: `"permutation"` or `"bayes"`.

- contrasts:

  Planned comparisons, or `NULL` for every pairwise contrast. A length-2
  character vector is one pairwise comparison. A list may mix length-2
  character vectors and lists with `collapse` (named remap of group
  levels) and/or `groups` (subset to these levels before testing). Names
  become `contrastId` values; unnamed elements get `A_vs_B` identifiers.

- nPermutations:

  Number of label permutations (permutation method only).

- nPosterior:

  Number of posterior draws (bayes method).

- nCores:

  Number of parallel workers (default 1).

- seed:

  Random seed for reproducibility.

- priorPseudocounts:

  Jeffreys prior increment per category (default 1/2).

- pAdjustMethod:

  Multiple-testing adjustment for permutation contrasts (`"holm"`,
  `"BH"`, `"none"`). Ignored for Bayes.

- credibleInterval:

  Bayes interval type: `"hpdi"` (default, shortest 95% interval on the
  draws) or `"cpi"` (equal-tailed 2.5%/97.5% quantiles). Ignored for
  permutation.

- includePairwise:

  Include every pairwise comparison in addition to `contrasts`. Default
  `NULL` includes pairwises when `contrasts` is omitted (the usual
  all-pairs analysis) and skips them when `contrasts` is supplied so
  Holm is paid only on the planned family. `FALSE` skips auto pairwises
  even when `contrasts` is omitted (omnibus only). `TRUE` keeps all
  pairwises even when `contrasts` is supplied (the old additive
  overlay).

- nope:

  Bayes NOPE interval on \\R\\. `"shuffle"` (default) sets a
  permutation-calibrated null reference interval whose upper bound comes
  from label shuffles within each posterior draw. A single number above
  1 sets a fixed NOPE bound on \\R\\ that you choose in advance. Either
  adds `R_PPGTnope` and `R_nopeHigh` to the output. `"none"` drops them.
  Ignored for permutation.

- energyNope:

  Bayes NOPE interval on the energy distance \\E\\, in squared Hellinger
  units. `"shuffle"` (default) sets a NOPE interval from the same label
  shuffles as `nope`. A single number above 0 sets a fixed NOPE bound on
  \\E\\;
  [`ComputeEnergyNopeBound()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/ComputeEnergyNopeBound.md)
  turns a composition shift into one. Either adds `E_PPGTnope` and
  `E_nopeHigh` to the output. `"none"` drops them. Ignored for
  permutation.

- nNopeShuffles:

  Label shuffles per posterior draw used to set the NOPEs (default 200).
  Each shuffle builds one group summary per shuffled statistic, so cost
  scales as `nPosterior * nNopeShuffles` summaries per contrast for one
  NOPE and twice that for both. Used only when `method = "bayes"` and
  `nope` or `energyNope` is `"shuffle"`.

- ...:

  Ignored.

## Value

A HellingerEnrichmentResult with omnibus summary, contrast table, and
(Bayes only) long-format posterior `draws` including the omnibus draws.
Bayes summaries carry `R_PPGTnope` and `R_nopeHigh` unless
`nope = "none"`, `E_PPGTnope` and `E_nopeHigh` unless
`energyNope = "none"` (each edge column is the median across draws of
the per-draw upper edge), and `enoughShufflesForNope` when `nope` or
`energyNope` is `"shuffle"` (`FALSE` when groups are too small for the
NOPE, which is why the NOPE probabilities are then `NA`); `draws`
carries the per-draw `R_nopeHigh` and `E_nopeHigh`. Both methods carry
`energyDistance` beside the ratio `effectSize`; Bayes adds
`energyCiLow`/`energyCiHigh` and per-draw `energyDistance` in `draws`.

## Details

The between/within Hellinger ratio, Jeffreys softening,
observed-inclusive label-permutation p-value, and collapse/subset
contrasts follow Paul Edlefsen's original procedure. `method = "bayes"`
instead draws subject compositions from a conjugate Dirichlet posterior,
rebuilds the Hellinger ratio \\R\\ on each draw, and reports the
posterior probability that \\R \> 1\\ (`PPGT1`) with a 95% highest
posterior density interval on \\R\\ by default.

`PPGT1` alone overstates enrichment when counts are deep. The Dirichlet
posterior only describes count noise within each subject, so as depth
grows it concentrates on the sample \\R\\, which is rarely exactly 1
even when groups do not differ. With `nope = "shuffle"` (the default),
the package reports a Null-Oriented Permutation Equivalence (NOPE)
interval on \\R\\: a permutation-calibrated null reference interval for
conditional effects. Within each posterior draw, the package shuffles
group labels on that draw's own distance matrix, so the interval is
conditional on the observed design (group sizes, depths, within-group
spread) and on that draw's compositions. The goal parallels a
permutation test, but NOPE reports a region instead of a p-value. The
NOPE interval for \\R\\ runs from 1 to the 95th percentile of the
shuffled \\R\\; `R_nopeHigh` is the median of that upper bound across
draws, and `R_PPGTnope` is the posterior probability that \\R\\ clears
it. `energyNope = "shuffle"` (also the default) builds the matching
interval for the energy distance \\E\\ from the same shuffles, from 0 to
`E_nopeHigh`, and reports `E_PPGTnope`.

NOPE borrows the equivalence logic of Kruschke's region of practical
equivalence (ROPE) and changes two things. In Kruschke's procedure the
analyst fixes the region before seeing the data, from domain knowledge
of which effects are negligible, and the decision compares the 95%
highest density interval with that region. Few analysts have that
intuition for a ratio of Hellinger distances or an energy distance, so
NOPE takes its upper bound from the data: the largest value of the
statistic these subjects routinely produce when labels carry no
information. The bound also differs from draw to draw, and each NOPE
probability counts the draws that clear their own bound. NOPE therefore
acts as a null reference calibrated by permutation, closer to a
permutation test than to a statement of practical significance. When you
know the smallest effect that matters, pass it as a fixed NOPE bound,
which plays the role of Kruschke's preset region:
`nope = <number above 1>` on \\R\\, or `energyNope = <number above 0>`
on \\E\\ in squared Hellinger units. A fixed bound needs no shuffles, so
it also works for groups too small for a NOPE interval.

To set a fixed bound on \\E\\ from a composition shift, pass the
reference and shifted compositions to ComputeEnergyNopeBound and hand
its result to `energyNope`. It returns \\2(1 - \mathrm{BC})\\, with BC
the Bhattacharyya coefficient between the two compositions, which equals
\\E\\ when every subject sits at its group's mean and approximates it
otherwise.

Because shuffling keeps the observed group sizes, depths, and
within-group spread, the NOPE widens or narrows with the data rather
than relying on a fixed number. Each draw's edge never falls below the
no-difference value (1 for \\R\\, 0 for \\E\\), so `R_PPGTnope <= PPGT1`
always holds, and `E_PPGTnope` never exceeds the posterior probability
that \\E \> 0\\. Reading a NOPE probability above 0.95 asks that nearly
every draw clear its own 95th shuffle percentile, which runs well below
a 5% false-positive rate and costs power; a looser threshold such as 0.5
behaves more like the permutation test. The sampling-depth vignette
measures both. Very small groups limit the NOPE the same way they limit
the permutation p-value. With 3 versus 3 subjects there are only 10
distinct splits, so the true split alone is 10% of all shuffles and
neither statistic can clear its 95th percentile. Whenever the true split
is at least 5% of shuffles (e.g. 2 versus 2, 3 versus 3, 2 versus 4),
the function sets each NOPE probability to `NA` with a warning rather
than a misleading 0, and `enoughShufflesForNope` is `FALSE` so the
reason travels with the result; 4 versus 4 (35 splits) is the smallest
balanced design that supports it.

Both methods also report the energy distance between groups
(`energyDistance`): twice the mean between-group squared Hellinger
distance minus each group's mean within-group squared distance, averaged
over group pairs for the omnibus. Squaring the existing distance matrix
costs nothing extra. \\R\\ measures how far apart the groups are
relative to how spread out they are; \\E\\ measures the shift in mean
composition itself, as the squared distance between group centroids in
square-root composition space, with no division. Independent count noise
adds equally to between- and within-group squared distances and cancels
from \\E\\ in expectation, but Jeffreys softening and the curvature of
the square root still shrink \\E\\ at shallow depth, though less than
they shrink the excess of \\R\\ over 1 (see the sampling-depth
vignette). \\E\\ responds mainly to shifts in mean composition; a
difference in spread alone moves it only slightly, likely through that
same curvature. The sample estimate excludes each subject's distance to
itself, so it can dip slightly below 0 when groups do not differ. The
p-values and `PPGT1` stay on \\R\\; under Bayes the package reports a
posterior interval on \\E\\ (`energyCiLow`, `energyCiHigh`) and, with
`energyNope`, `E_PPGTnope`.

## References

Edlefsen, P. Hellinger between/within enrichment for subject-level
categorical compositions (Jeffreys softening, ratio statistic, label
permutation, collapse/subset contrasts).

Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. Advances in Methods and Practices in Psychological
Science, 1(2), 270-280.
[doi:10.1177/2515245918771304](https://doi.org/10.1177/2515245918771304)
. Article page:
<https://www.psychologicalscience.org/journals/ampps/2515245918771304/>;
open-access record:
<https://scholarworks.iu.edu/dspace/items/bada026d-fe1b-4f53-ad5b-db6e8bf99711>
