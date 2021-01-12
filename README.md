
> <br/> Kristof Haneca<br/> 12 januari, 2021<br/>
> 
> [![](./figures/RG.png)](https://www.researchgate.net/profile/Kristof_Haneca)

# dendro\_functions

A set of functions for tree-ring analysis.

## `CorrTable`

Computes correlations between two sets of tree-ring series (class ‘rwl’,
see `dplR`). When only one set of trs is provided, internal correlations
are computed.

Output values:

  - overlap: n° of shared years between two series
  - glk: the *Gleichlaüfigkeit* or percentage of parallel variation
    (%PV)
  - tBP: t-values according to the Baillie-Pilcher algorithm
  - tHo: t-values accoriding to the Hollstein 1980 algorithm
  - r\_pearson: Pearson correlation between two series
