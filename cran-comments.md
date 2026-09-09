# CRAN submission comments — attrition 1.0.0

> **STALE as of 2026-09-08. Do not submit on these results.** Everything below
> describes a tarball built from `e9a6677`. The estimator output names changed
> after that commit (see NEWS), so the win-builder rounds recorded here checked
> code that no longer exists. Re-run step 3 of the release protocol and rewrite
> this file before submitting.

## Test environments

* macOS 26.6.2 (local), R 4.6.0
* GitHub Actions, all passing, tests FAIL 0 | WARN 0 | SKIP 0 | PASS 278 on each:
  * macOS-latest (release)
  * Windows-latest (release)
  * Ubuntu-latest (R-devel, release, oldrel-1)
* win-builder, x86_64-w64-mingw32, both Status: 1 NOTE (the note below):
  * R Under development (unstable) (2026-08-27 r90452 ucrt)
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

## Reverse dependencies

None (first release).
