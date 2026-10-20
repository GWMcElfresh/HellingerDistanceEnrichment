# Resolve whether auto pairwises run given planned contrasts.

NULL follows the planned-comparison rule: all pairwises when `contrasts`
is omitted, none when a family is supplied. FALSE is omnibus-only if no
family is supplied; TRUE restores the old additive overlay.

## Usage

``` r
resolve_include_pairwise(include_pairwise, has_custom_contrasts)
```
