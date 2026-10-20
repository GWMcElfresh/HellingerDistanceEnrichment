# Stack the ratio and energy distance into one long table for faceting.

Each input row (a posterior draw or a contrast summary) appears once per
shown statistic, with the statistic's value in `value` and its panel
label in `statistic`. The ratio panel always comes first.

## Usage

``` r
stack_plot_statistics(source_data, show_energy, statistic_labels)
```

## Arguments

- source_data:

  Draws (Bayes) or contrast rows (permutation) carrying `effectSize`
  and, when `show_energy` is TRUE, `energyDistance`.

- show_energy:

  Whether to include the energy distance panel.

- statistic_labels:

  Named panel labels (`ratio`, `energy`).

## Value

data.frame with `value` and a factor `statistic`.
