# CRAN submission comments — attrition 1.0.0

## Test environments

* macOS 26.5 (local), R 4.6.0
* GitHub Actions, all passing:
  * macOS-latest (release)
  * Windows-latest (release)
  * Ubuntu-latest (R-devel, release, oldrel-1)
* win-builder, R Under development (unstable) (2026-08-08 r90381 ucrt),
  x86_64-w64-mingw32

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the expected one for a first submission:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Alexander Coppock <acoppock@gmail.com>'

    New submission

    Possibly misspelled words in DESCRIPTION:
      Coppock (11:65)
      Imbens (15:5)
      Manski (14:6, 15:12)
      Nonignorable (2:19)
      nonignorable (13:22)
      poststratification (17:20)

All six flagged words are spelled correctly. Coppock, Imbens, and Manski are
surnames of cited authors. "Nonignorable" and "poststratification" are standard
terms in the missing-data and survey-sampling literatures, and both appear in
the titles of the works cited in the Description.

## Reverse dependencies

None (first release).
