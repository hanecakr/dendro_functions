
> <br/> Kristof Haneca<br/> 12 mei, 2021<br/>
> 
> [![](./figures/RG.png)](https://www.researchgate.net/profile/Kristof_Haneca)

# dendro\_functions

A set of functions for tree-ring analysis.

The [:file\_folder: R](/R) directory contains the R-script for each
function.

## `CorrTable`

Computes correlations between a set of tree-ring series `x` (class
‘rwl’, see `dplR`: <https://github.com/AndyBunn/dplR>) and dated
reference chronologies `y`.

When only one set of tree-ring series is provided, internal correlations
are computed.

``` r
CorrTable(x,
          y = NULL,
          min_overlap = 30,
          remove.duplicates = TRUE,
          output = "table", #c("matrix", "table") a list of matrices for each correlation variable, or a single table as output
          values = c("r_pearson", "t_St", "glk", "glk_p", "tBP", "tHo"), # for selecting a reduced set of correlation variables
          sort_by = "tHo") #c("glk", "t_St", "tBP", "tHo", "r_pearson")) sorting variable when "table" is chosen as output
```

**Output values:**

  - overlap: n° of shared years between two series
  - glk: the *Gleichlaüfigkeit* or percentage of parallel variation
    (%PV)
  - glk\_p: the probability of exceedence for the *Gleichläufigkeit*
    value
  - tBP: *t*-values according to the Baillie-Pilcher (1973) algorithm
  - tHo: *t*-values according to the Hollstein (1980) algorithm
  - r\_pearson: Pearson correlation between two series
  - t\_St: Student’s *t*-value of the Pearson correlation
