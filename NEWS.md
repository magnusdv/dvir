# dvir 3.3.0

* `dviSolve()` and `amDrivenDVI()` now uses the *generalised likelihood ratio* (GLR) for identifications in families with more than one missing person.

# dvir 3.2.1

## New features
*	Add file `dvi-example.fam` as a system file. Hence avoid internet download during example (gave CRAN error).
*	New function `dviSolve()` implementing a complete pipeline for solving a DVI case.
*	New function `amDrivenDVI()` implementing AM-driven analysis.
* New functions `setPairing()` and `excludePairing()` for manually fixing or disallowing certain pairings.
* In `findUndisputed()` replace argument `relax` with its negation `strict`, and reverse the default to be `strict = FALSE`.
*	Harmonise the output summaries of the main functions.
*	New function `dviSim()` for simulating marker data onto a DVI dataset.
*	New function `findNonidentifiable()`.
*	Export (experimental) `checkDVI()`.
*	`mergePM()` gains argument `verbose`.

## Other
*	Enforce family names in `consolidateDVI()`.
* Reorganise and synchronise code in `findExcluded()` and `exclusionMatrix()`.
*	Use `cat()` instead of `message()` throughout.
*	Remove deprecated `summariseDVI()` (replaced by `print.dviData()`).

## Bug fixes
*	`findUndisputed()` miscounted steps if `verbose = F`.
*	`jointDVI()` sometimes dropped victims after `findUndisputed()`.


# dvir 3.1.0

* The __dvir__ package is now maintained by Magnus D Vigeland.

* In `print.dviData()`, add info on sex of victims/missing.

* Add `report` to output of `findExcluded()`.

* New function `plotUndisputed()`.

# dvir 3.0.1

* Fix bad URL reported by CRAN.


# dvir 3.0.0

Version 3.0.0 constitutes a major rewrite of **dvir**, with many new features reflecting a broader scope of the package. Furthermore,
the syntax has been greatly simplified, due to the new `dviData` container class for DVI datasets. 

It should be noted that these syntax changes are not backwards compatible.


### Breaking changes

* Most functions of **dvir** now expects a `dviData` object as input.

* All datasets have been regenerated as `dviData` objects. 

* Datasets from the book Mass Identifications (Kling, Egeland, Tillmar, Prieto) have been renamed with prefix KETP.

* `summariseDVI()` is deprecated in favour of the new `print()` method for `dviData` objects.

* As of version 3.0.0, **dvir** requires R >= 4.1.0 and recent versions of **pedtools**, **forrel** and **pedprobr**.


### New features

* New S3 class `dviData`, and constructor with syntax `dviData(pm, am, missing)`.

* New function `plotDVI()` for visualising DVI datasets.

* New function `relabelDVI()` greatly simplifying relabelling tasks.

* New function `findExcluded()` for identifying victim samples not matching any of the missing persons - and vice versa.

* New functions `directMatch()` and `mergePM()` for analysing and merging victim samples coming from the same individual.

* New function `plotSolution()` helps visualising the output of `jointDVI()`.

* New function `getFamily()` for extracting the family name (or index) of missing persons.

* New internal function `consolidateDVI()` ensures well-formed datasets. This is called in the beginning of all major functions.

* New function `compactJointRes()` simplifies the output of `jointRes()`.

* New function `familias2dvir()` for parsing `.fam` files written by the DVI module of Familias.

* `sequentialDVI()` has more informative output, in sync with `findUndisputed()`. While these functions do almost the same thing, the latter is generally preferred in practice.

* `jointDVI()` gains arguments `nkeep` and `maxAssign`. The latter triggers a gracefully exit if the number of assignments is too large.

* `jointDVI()` now has a progress bar (but only when `numCores = 1`).


### Other

* Revamped docs and examples.
* Several functions are more efficient due to better code organisation.


# dvir 2.2.0

* This is intended to be the last version in the 2.x series. Version 3 will include some breaking changes in syntax.


# dvir 2.1.0

* Various improvements of README and documentation.


# dvir 2.0.0

* Initial CRAN release.
