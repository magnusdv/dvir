
<!-- README.md is generated from README.Rmd. Please edit that file -->
The dvir (Disaster Victim Identification) library
=================================================

DNA profiles are available from victims (post mortem, pm, data) and reference families (ante mortem, am, data) with missing persons. A forward stepwise approach, conceptually similar to variable selection in regression analysis, is implemented in the function `forward`.

Installation
------------

To get the lastest version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install mut2 from GitHub
devtools::install_github("thoree/dvir")
```

The implementation relies heavily on the `pedtools` suite of R-libraries, in particular the `forrel` library which can be installed by running

``` r
devtools::install_github("magnusdv/forrel")
```

Example 1
---------

The data is simulated and summarised by the below figure

<img src="man/figures/dvi.data.png" >
-------------------------------------

One of the 13 CODIS, TPOX, is included. Five victims are shown to the left. The reference families are on the right hand side. Based on this marker, we see for instance that the only possibility for `MP3` is `V3` as mutations are disregarded in this example. The preparation of data, including modelling of mutations, is exemplified in the documentation of the main function `forward`. The below data has been prepared and is loaded and analysed below:

``` r
library(dvir)
library(forrel)
data(dvi.data)
am = dvi.data$am
pm = dvi.data$pm
vp = dvi.data$vict
mp = dvi.data$miss
res = forward(pm, am, vp, mp)
res[[3]]
#>   from   to           lik steps
#> 0 <NA> <NA> 2.625697e-132     0
#> 1   V2  MP2 1.525893e-129     1
#> 2   V3  MP3 1.660819e-120     2
#> 3   V1  MP1 4.418114e-117     1
```

The first line of the table only gives the null likelihood, the likelihood prior to any attempt of identification. This can be calculated directly using only the `forrel`library

``` r
prod(forrel::LR(list(pm, am), 1)$likelihoodsPerSystem)
#> [1] 2.625697e-132
```

In the next steps there are 5\*3 = 15 moves that can be made. A `move` is an assignment of one victim to one missing person, such as e.g. `M3 := V3`. The 15 likelihoods are calculated and the maximum is given in line two of the table. The LR compares `M3 := V3` to `V3-unrelated` and equals

``` r
2.857873e-123/2.625697e-132
#> [1] 1088424521
```

An identification, summarised in a line in the table, is included if this LR exceeds a prescribed level determined by the variable `LRlimit`in `dviForward`. The default is `1`. It is by no means obvious how specify the next column, the `prior`. Here we use a flat prior. The sample space consists of the 15 possible moves and one option corresponding to no indentification. The flat prior used is `1/(15+1) = 0.625`. Finally,the posterior is found using Bayes theorem. The complete table gives the result corresponding to the one from which data was simulated. Observe that the LR may increase. In the above table this is the case. The reason is that it is easier to identify `V1` given that a brother has already been identified.

It's relevant to study how the procedure works when the victims are unrelated to all missing persons. This can be studied by simulation ...

\`\`\`

Example 2
=========

``` r
library(dvir)
library(forrel)
data(dvi.nfi)
res = forward(from = pm, to = am, ids.from = vict, ids.to = miss, eliminate = FALSE, 
singleStep = TRUE, twoStep = FALSE)
res[[2]][[1]]
#>       from   to           lik           LR prior posterior
#> step0 <NA> <NA> 7.902327e-181           NA    NA        NA
#> step1   V1  MP1 3.740769e-170 4.733756e+10    NA 0.9999969
#> step2   V3  MP3 1.282852e-164 3.429381e+05    NA 0.9529405
#> step3   V2  MP2 1.608944e-157 1.254193e+07    NA 0.9833694
#> step4   V5  MP5 1.095323e-147 6.807710e+09    NA 0.9999840
#> step5   V4  MP4 2.170296e-144 1.981422e+03    NA 0.9991241
```
