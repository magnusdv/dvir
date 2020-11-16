
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The dvir (Disaster Victim Identification) library

We assume DNA profiles are available from victim samples (post mortem,
pm data) and reference families (ante mortem, am, data) with missing
persons (MP-s). There may be several samples from the same victim,
potentially of low quality leading to *drop-outs*. The problem is to
identify the MP-s. Some (or all) victims may not be among the MP-s.
Similarly, there may be MP-s not in the list of victims. A search
strategy is implemented. All victims are initially tried, one at a time,
in all MP positions. Results are sorted according to the likelihood and
assignments with a LR (compared to the null likelihood) below a user
specified limit are omitted from further search. If mutations are
modelled all LR-s will typically be positive and the limit must be
specified to a negative number to include all possibilities in the
future search. Based on this initial screening, all possible assignments
of victims are generated. Note that only a subset, possibly none, of
victims may be mapped to MP-s. The resulting list of assignments may be
prohibitively large and for this reason it possible to restrict the
search by specifying that only the `nbest` assignments for each victim
be considered.

## Installation

To get the latest version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install dvir from GitHub
devtools::install_github("thoree/dvir")
```

The implementation relies heavily on the `ped suite` of R-libraries, in
particular `forrel` and `pedmut`. These are automatically installed by
the above command.

## Load libraries

We start by loading **dvir**. We also load **pedtools** for creating and
plotting pedigrees.

``` r
library(dvir)
library(pedtools)
```

## Example 1

We consider the following example

``` r
# Attributes of a single marker
locAttr = list(name = "m", alleles = 1:3, afreq = c(1, 1, 1)/3)

# PM data (victims)
n = 7
ids.from = paste0("V", 1:n)
sex = c(rep(1, n-1), 2)
df = data.frame(famid = ids.from, id = ids.from, fid = 0, mid = 0, sex = sex,
                m = c("1/1", "2/2", "1/1", "1/1", "2/2", "2/2", "2/2"))
from = as.ped(df, locusAttributes = locAttr)

# AM data (families)
MPs = c("MP1", "MP2", "MP3")
to = nuclearPed(3, father = "R1", mother = "R2", children = MPs)
m = marker(to, "R1" = "1/1", "R2" = "1/1", name = "m")
to = setMarkers(to, m, locusAttributes = locAttr)
```

``` r
# Plot both
plotPedList(list(from, to), marker = 1, 
            hatched = typedMembers, 
            col = list(red = MPs), 
            titles = c("PM data. 7 victims", "AM data. 3 MP-s"))
```

![](man/figures/README-ex1-ped-1.png)<!-- -->

We do not consider mutations or other artefacts. We assume that copies
of victim samples have been identified and merged so that without extra
information like (*todo:* add age), there are six symmetric solutions.
If the ages are known, for instance age(V1) \> age(V2) \> age(V3), the
solution may be unique. Recall the convention that individuals are
ordered left to right in the pedigree based on age.

### The number of assignments

The number of possible assignments is:

``` r
ncomb(nVfemales = 1, nMPfemales = 0, nVmales = 6, nMPmales = 3)
#> [1] 229
```

The complete list of these assignments is generated as follows (see also
the function `generateMoves()` below).

``` r
moves = list(V1 = c("V1", MPs ), V2 = c("V2", MPs ), V3 = c("V3", MPs ), V4 = c("V4", MPs),
             V5 = c("V5", MPs ), V6 = c("V6", MPs ), V7 = "V7" )
a = expand.grid.nodup(moves)
length(a)
#> [1] 229
```

### The search

The following search ranks all possible solutions by their likelihood.
The ten best solutions are shown.

``` r
res = global(from, to, MPs, moves = NULL, limit = -1, verbose = F)
res[1:10, ]
#>     V1 V2  V3  V4 V5 V6 V7     loglik  LR  posterior
#> 1  MP3 V2 MP2 MP1 V5 V6 V7  -8.788898 729 0.12326682
#> 2  MP2 V2 MP3 MP1 V5 V6 V7  -8.788898 729 0.12326682
#> 3  MP3 V2 MP1 MP2 V5 V6 V7  -8.788898 729 0.12326682
#> 4  MP1 V2 MP3 MP2 V5 V6 V7  -8.788898 729 0.12326682
#> 5  MP2 V2 MP1 MP3 V5 V6 V7  -8.788898 729 0.12326682
#> 6  MP1 V2 MP2 MP3 V5 V6 V7  -8.788898 729 0.12326682
#> 7  MP2 V2 MP1  V4 V5 V6 V7 -10.986123  81 0.01369631
#> 8  MP3 V2 MP1  V4 V5 V6 V7 -10.986123  81 0.01369631
#> 9  MP1 V2 MP2  V4 V5 V6 V7 -10.986123  81 0.01369631
#> 10 MP3 V2 MP2  V4 V5 V6 V7 -10.986123  81 0.01369631
```

Next, we exemplify how to limit the search, which may be necessary in
larger cases. First all sex-consistent marginal moves are generated.

``` r
moves = generateMoves(from, to, MPs)
```

Keep only the three best marginal candidates for each victim.

``` r
moves2 = marginal(from, to,  MPs, moves, limit = -1, sorter = T, nkeep = 3)
res = global(from, to, MPs, moves = moves2[[1]], limit = -1, verbose = F)
head(res)
#>    V1 V2  V3  V4 V5 V6 V7    loglik  LR posterior
#> 1 MP3 V2 MP2 MP1 V5 V6 V7 -8.788898 729 0.1666667
#> 2 MP2 V2 MP3 MP1 V5 V6 V7 -8.788898 729 0.1666667
#> 3 MP3 V2 MP1 MP2 V5 V6 V7 -8.788898 729 0.1666667
#> 4 MP1 V2 MP3 MP2 V5 V6 V7 -8.788898 729 0.1666667
#> 5 MP2 V2 MP1 MP3 V5 V6 V7 -8.788898 729 0.1666667
#> 6 MP1 V2 MP2 MP3 V5 V6 V7 -8.788898 729 0.1666667
```

## Example 2

We now consider a larger dataset, loaded as follows:

``` r
load(url("http://familias.name/BookKETP/Files/Grave.RData"))
```

The loaded dataset contains objects `from`, `to`, `ids.to` and `moves`.

``` r
summary(from)
#> List of 8 singletons.
#> Labels: V1 (female), V2 (male), V3 (female), V4 (female), V5 (female), V6 (female), V7 (male), V8 (male).
#> 23 attached markers.
summary(to)
#> Pedigree with 23 members.
#> 23 attached markers.
#> 5 typed members.
```

The family `to` has 8 missing persons, labelled MP1-MP8, and 5 genotyped
family members, labelled R1-R5. The pedigree is shown below.

``` r
refs = paste0("R", 1:5)
mps = paste0("MP", 1:8)
plot(to, title = "AM data", labs = c(refs, mps), hatched = c(refs, mps), 
     col = list(red = refs, blue = mps), deceased = mps)
```

![](man/figures/README-ex2-ped-1.png)<!-- -->

The list of singletons in `from` contains female victims V1, V3, V4, V5,
V6, and male victims V2, V7, V8. The *a priori* possible number of
assignments, ignoring symmetries, is

``` r
ncomb(nVfemales = 5, nMPfemales = 5, nVmales = 3, nMPmales = 3)
#> [1] 52564
```

We restrict the number of assignments by requiring \(LR > 0.99\) for
*marginal* moves. For instance, based on the below, the possibility `V1
= MP1` will be considered since the \(LR\) comparing this assignment
(and no further victims identified) to the null hypothesis (no victims
identified) exceeds 0.99.

``` r
moves = generateMoves(from, to, ids.to)
m = marginal(from, to, ids.to, limit = 0.99, moves = moves,  sorter = T)
m[[1]]
#> $V1
#> [1] "MP1" "V1" 
#> 
#> $V2
#> [1] "MP2" "V2" 
#> 
#> $V3
#> [1] "MP3" "V3" 
#> 
#> $V4
#> [1] "MP4" "MP5" "V4" 
#> 
#> $V5
#> [1] "MP4" "MP5" "V5" 
#> 
#> $V6
#> [1] "MP6" "V6" 
#> 
#> $V7
#> [1] "MP7" "MP8" "V7" 
#> 
#> $V8
#> [1] "V8"
```

The complete list of marginal LR values are contained in `m[[2]]`, shown
below.

``` r
m[[2]]
#> [[1]]
#>       MP1        V1       MP3       MP4       MP5       MP6 
#> 479971259         1         0         0         0         0 
#> 
#> [[2]]
#>          MP2           V2          MP8          MP7 
#> 6.776011e+10 1.000000e+00 5.512209e-01 0.000000e+00 
#> 
#> [[3]]
#>          MP3           V3          MP1          MP4          MP5          MP6 
#> 6.409841e+14 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#> 
#> [[4]]
#>        MP4        MP5         V4        MP1        MP3        MP6 
#> 1.8036e+12 1.8036e+12 1.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 
#> 
#> [[5]]
#>          MP4          MP5           V5          MP1          MP3          MP6 
#> 103006682220 103006682220            1            0            0            0 
#> 
#> [[6]]
#>          MP6           V6          MP1          MP3          MP4          MP5 
#> 8.817392e+12 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#> 
#> [[7]]
#>          MP7          MP8           V7          MP2 
#> 16946051.333      295.839        1.000        0.000 
#> 
#> [[8]]
#>       V8      MP8      MP2      MP7 
#> 1.000000 0.268489 0.000000 0.000000
```

From the above this we realise that `limit = 0` might be a better
option, since it would only add `V8 = MP8` and `V2 = MP8`. However,
changing to `limit = 0` would increase computation time considerably.

We perform the search using the assignments in `m[[1]]`:

``` r
res1 = global(from, to, ids.to, limit = 0.99, moves = m [[1]])
head(res1)
#>    V1  V2  V3  V4  V5  V6  V7 V8    loglik           LR    posterior
#> 1 MP1 MP2 MP3 MP4 MP5 MP6 MP7 V8 -737.0038 1.290829e+95 9.999998e-01
#> 2 MP1 MP2 MP3 MP4 MP5 MP6  V7 V8 -752.3418 2.816133e+88 2.181646e-07
#> 3  V1 MP2 MP3 MP4 MP5 MP6 MP7 V8 -768.7262 2.157791e+81 1.671632e-14
#> 4 MP1  V2 MP3 MP4 MP5 MP6 MP7 V8 -773.6762 1.528448e+79 1.184082e-16
#> 5 MP1 MP2  V3 MP4 MP5 MP6 MP7 V8 -774.0113 1.093148e+79 8.468571e-17
#> 6 MP1 MP2 MP3  V4 MP5 MP6 MP7 V8 -774.8047 4.944458e+78 3.830451e-17
```

We check the assignment with the identification `MP8 = V8` added.

``` r
res2 = global(from, to, ids.to, 
       moves = list(V1 = "MP1", V2 = "MP2", V3 = "MP3", V4 = "MP4",
                    V5 = "MP5", V6 = "MP6", V7 = "MP7", V8 = "MP8"))
res2
#>    V1  V2  V3  V4  V5  V6  V7  V8    loglik           LR posterior
#> 1 MP1 MP2 MP3 MP4 MP5 MP6 MP7 MP8 -737.8061 5.786551e+94         1
exp(res2$loglik - res1$loglik[1])
#> [1] 0.4482818
```

Finally, the code for performing an exhaustive search (i.e., `limit
= 0`) is shown below. With parallelisation (activated by setting the
argument `numCores` larger than 1) this computation should take around
15-20 minutes.

``` r
res3 = global(from, to, ids.to, moves = NULL, numCores = 4)
head(res3)
#>    V1  V2  V3  V4  V5  V6  V7  V8    loglik           LR    posterior
#> 1 MP1 MP2 MP3 MP4 MP5 MP6 MP7  V8 -737.0038 1.290829e+95 6.904732e-01
#> 2 MP1 MP2 MP3 MP4 MP5 MP6 MP7 MP8 -737.8061 5.786551e+94 3.095266e-01
#> 3 MP1 MP2 MP3 MP4 MP5 MP6  V7  V8 -752.3418 2.816133e+88 1.506369e-07
#> 4 MP1 MP2 MP3 MP4 MP5 MP6  V7 MP8 -753.3430 1.034770e+88 5.535057e-08
#> 5  V1 MP2 MP3 MP4 MP5 MP6 MP7  V8 -768.7262 2.157791e+81 1.154217e-14
#> 6  V1 MP2 MP3 MP4 MP5 MP6 MP7 MP8 -769.5285 9.672985e+80 5.174146e-15
```
