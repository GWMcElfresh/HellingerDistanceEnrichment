# Simulate a two-group cohort with controlled depth, heterogeneity, and shift.

Each subject draws a latent composition from a Dirichlet centered on its
group's mean composition, then draws `reads_per_subject` counts from a
multinomial on that composition. Three knobs separate the sources of
variation the sampling vignette studies: `reads_per_subject` sets count
noise within a subject, `concentration` sets how far subjects scatter
around their group mean (larger is more homogeneous), and `shift` moves
composition mass between two categories in the second group only.

## Usage

``` r
simulate_composition_cohort(
  n_subjects_per_group = 8,
  n_categories = 10,
  reads_per_subject = 2000,
  concentration = 50,
  shift = 0,
  groups = c("Control", "Treatment"),
  seed = NULL
)
```

## Arguments

- n_subjects_per_group:

  Subjects per group; length 1 or 2.

- n_categories:

  Number of categories (\>= 2).

- reads_per_subject:

  Multinomial count total per subject.

- concentration:

  Dirichlet concentration per group; length 1 or 2.

- shift:

  Composition mass moved in the second group, in \[0, first category's
  baseline share).

- groups:

  Two group labels.

- seed:

  Random seed, or `NULL` to use the current stream.

## Value

List with `long_table` (subjectId, category, group, n), `latent_counts`
(subjects by categories latent compositions scaled to
`reads_per_subject`, for noise-free reference statistics), and `group`
(named factor).

## Details

The baseline mean composition decays linearly over categories (e.g. 10
categories run from about 0.18 down to 0.02), so the shift moves `shift`
mass from the first, most abundant category to the last, rarest one.
