# Default ALR pivot pool for circulating-miRNA panels (v1)

Returns the curated v1 pivot pool used as the geometric-mean denominator
by
[`ws_alr_pivot`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md).
Six low-variance circulating miRNAs with no documented contamination
from haemolysis (miR-451a/16/486/144/223) or platelet activation
(miR-223 family). The pool is fixed at v1; any revision will be tagged
`circulating_v2`, etc., with a separate justification log under the
single-sample scoring bank's pivot-pool justification notes.

Notably absent from the v1 pool: miR-92a-3p (recent serum/plasma
literature treats it as erythrocyte-derived; codex Round 1 review of the
single-sample scoring bank plan v0.2 flagged its prior inclusion as a
concern).

## Usage

``` r
ws_default_pivot_pool()
```

## Value

Character vector of six miRNA IDs.
