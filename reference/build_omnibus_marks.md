# Omnibus line positions and labels for each shown panel.

The ratio panel's label carries the decision statistics (the permutation
p-value, or `PPGT1` and `R_PPGTnope` under Bayes); the energy panel's
label carries \\E\\ and, under Bayes, `E_PPGTnope`. Each label sits on
whichever side of its line has more room within the panel, counting the
reference line, so the panel edge does not clip it.

## Usage

``` r
build_omnibus_marks(result, long_data, statistic_labels, shown_labels)
```

## Arguments

- result:

  A HellingerEnrichmentResult.

- long_data:

  Output of stack_plot_statistics.

- statistic_labels:

  Named panel labels (`ratio`, `energy`).

- shown_labels:

  Panel labels present in the plot.

## Value

data.frame with `statistic`, `value`, `label`, and `hjust`.
