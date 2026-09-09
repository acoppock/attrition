# CRAN submission comments — attrition 1.0.0

## Test environments

* macOS 26.6.2 (local), R 4.6.0
* GitHub Actions, all passing, tests FAIL 0 | WARN 0 | SKIP 0 | PASS 278 on each:
  * macOS-latest (release)
  * Windows-latest (release)
  * Ubuntu-latest (R-devel, release, oldrel-1)
* win-builder, x86_64-w64-mingw32, both Status: 1 NOTE (the note below):
  * R Under development (unstable) (2026-09-08 r90509 ucrt)
  * R 4.6.1 (2026-06-24 ucrt)

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

## Suggested packages

The second vignette, "Drawing the extreme value bounds", uses vayr and estimatr,
both in Suggests. Its chunks are guarded on both packages being available and on
vayr being at least 1.1.0, the version that added the two functions it calls, so
the vignette builds either way. On both win-builder rounds the guard passed and
the figure rendered, which was confirmed by finding the embedded PNG in
doc/drawing-the-bounds.html inside the built binary rather than by assuming it.

## Reverse dependencies

None (first release).
