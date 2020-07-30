
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The dvir (Disaster Victim Identification) library

DNA profiles are available from victims (post mortem, pm, data) and
reference families (ante mortem, am, data) with missing persons (mp-s).
The problem is to identify the mp-s. Some (or all) victims may not be
among the mp-s. Similarly, there may be mp-s not in the list of victims.
A search strategy is implemented. All victims are initially tried, one
at a time, in all mp positions. Results are sorted according to the
likelihood and moves with a LR (compared to the null likelihood) below a
user specified limit are omitted from further search. If mutations are
modelled all LR-s will typically be positive and the limit must be
specified to a negative number to include all possibilities in the
future search. Based on this initial screening, all possible moves of
victims are generated. Note that only a subset, possibly none, of
victims may be mapped to mp-s. The resulting list of moves may be
prohibitively large and for this reason it possible to restrict the
search by specifying that only the `nbest` moves for each victim be
considered.

In addition, a forward stepwise approach, conceptually similar to
variable selection in regression analysis, is implemented in the
function `forward`. This approach is generally faster, but may fail to
find a solution.

## Installation

To get the lastest version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install mut2 from GitHub
devtools::install_github("thoree/dvir")
```

The implementation relies heavily on the `pedtools` suite of
R-libraries, in particular the `forrel` and `pedmut` libraries which can
be installed by running

``` r
devtools::install_github("magnusdv/forrel")
devtools::install_github("magnusdv/pedmut")
```

## Example 1

The number of combinations

``` r
library(dvir)
ncomb(3,3,1,1)
#> [1] 68
```
