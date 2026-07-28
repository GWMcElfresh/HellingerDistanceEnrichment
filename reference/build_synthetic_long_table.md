# Build a long counts table with exchangeable multinomial compositions per group.

Each subject receives an independent symmetric multinomial draw over
categories within its group. Before planting, all groups share the same
generative process, so group labels are exchangeable under the null.

## Usage

``` r
build_synthetic_long_table(
  n_subjects_per_group = 5,
  n_categories = 4,
  groups = c("Control", "Treatment"),
  seed = 1
)
```

## Arguments

- n_subjects_per_group:

  Number of subjects per group level (\>= 1).

- n_categories:

  Number of category levels (Category0, Category1, ...).

- groups:

  Character vector of distinct group labels (\>= 2 levels).

- seed:

  Random seed for reproducibility.

## Value

data.frame with columns subjectId, category, group, n.
