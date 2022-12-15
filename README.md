
> <br/> Kristof Haneca<br/> 15 december, 2022<br/>
>
> [![](./figures/RG.png)](https://www.researchgate.net/profile/Kristof_Haneca)

# dendro_functions

A set of functions for tree-ring analysis.

The [:file_folder: R](/R) directory contains the R-script for each
function.

## `corrTable`

!! *Test version. Currently under development* !!

Computes correlations between a set of tree-ring series `x` (class
‘rwl’, see `dplR`: <https://github.com/AndyBunn/dplR>) and dated
reference chronologies `y`.

When only one set of tree-ring series is provided, internal correlations
are computed.

``` r
corrTable(x,
          y = NULL,
          min_overlap = 30,
          remove.duplicates = TRUE,
          output = "table",
          values = c("r_pearson", "t_St", "glk", "glk_p", "tBP", "tHo"), 
          sort_by = "tHo")
```

**Function parameters:**

-   x: tree-ring data as data.frame of class ‘rwl’
-   y: tree-ring data (or reference chronologies) as data.frame of class
    ’rwl
-   min_overlap: numeric, minimum overlap of series and chronologies
-   remove.duplicates: logical, if TRUE duplicate comparisons between
    pairs of series/chronologies are removed
-   output: should be one of c(“matrix”, “table”)
-   values: which correlation values should be computed: c(“r_pearson”,
    “t_St”, “glk”, “glk_p”, “tBP”, “tHo”)
-   sort_by: correlation value to sort resulting table by

**Output values:**

-   overlap: n° of shared years between two series
-   glk: the *Gleichläufigkeit* or percentage of parallel variation
    (%PV)
-   glk_p: the probability of exceedence for the *Gleichläufigkeit*
    value
-   tBP: *t*-values according to the Baillie-Pilcher (1973) algorithm
-   tHo: *t*-values according to the Hollstein (1980) algorithm
-   r_pearson: Pearson correlation between two series
-   t_St: Student’s *t*-value of the Pearson correlation
