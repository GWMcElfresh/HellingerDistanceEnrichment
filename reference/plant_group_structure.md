# Plant a composition shift in one group-category cell.

Multiplies counts at (target_group, target_category) by boost. Other
group levels are unchanged, so a multi-group design can plant enrichment
in one arm (for example Condition1) while leaving another arm
exchangeable with Control (for example Condition2).

## Usage

``` r
plant_group_structure(
  long_table,
  target_group = "Treatment",
  target_category = "Category0",
  boost = 5
)
```

## Arguments

- long_table:

  Output from build_synthetic_long_table.

- target_group:

  Group label that receives elevated counts.

- target_category:

  Category label to enrich.

- boost:

  Positive multiplier applied only to matching rows.

## Value

Modified long table with the planted shift.
