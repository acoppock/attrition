# Perceived polarization under double sampling

A two-wave survey experiment on Amazon Mechanical Turk, replicating
Levendusky and Malhotra (2016) and reported as the application in
Coppock, Gerber, Green, and Kern (2017). Subjects read a news article
describing the electorate either as sharply divided (the polarized
condition) or as focused on common ground (the moderate condition).
Outcomes were measured immediately in Wave 1 and again in a Wave 2
survey ten days later.

## Usage

``` r
levendusky_replication
```

## Format

A tibble with 1,980 rows and 10 columns:

- X_party_id:

  Party identification: `Dem`, `Ind`, or `Rep`. The poststratification
  variable used in Table 3.

- Z_condition:

  The condition as assigned: `Moderate` or `Polarized`.

- Z:

  Polarized (1) versus moderate (0). The contrast analyzed throughout.

- R1:

  Responded in the Wave 2 initial sample.

- Attempt:

  Selected for the follow-up sample and offered the higher incentive.

- R2:

  Responded to the follow-up attempt.

- Y_polarization_w1:

  Perceived polarization at Wave 1, from 0 to 6.

- Y_polarization_w2:

  Perceived polarization at Wave 2, from 0 to 6. The dependent variable
  throughout, and the one with missing values.

- Y_extremity_w1:

  Perceived extremity at Wave 1, from 0 to 3.

- Y_extremity_w2:

  Perceived extremity at Wave 2, from 0 to 3.

Perceived polarization is built from a battery of policy questions.
Subjects gave their own view and then guessed how a typical Democratic
voter and a typical Republican voter would answer. The outcome is the
average absolute difference between the two guesses.

The response indicators rather than the missing values define who
responded. Sixteen subjects have a Wave 2 outcome recorded in the
archive despite `R1 == 0` and `Attempt == 0`, which is why 448 outcomes
are `NA` where 536 subjects did not respond. Every estimator here keys
on `R1`, `Attempt`, and `R2`, as the paper does, so those sixteen
outcomes go unused.

## Source

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2016). Replication Data for: Combining double sampling and bounds to
address non-ignorable missing outcomes in randomized experiments.
Harvard Dataverse.
[doi:10.7910/DVN/AQB4MP](https://doi.org/10.7910/DVN/AQB4MP) , file
2887314 (`levendusky_mturk_clean.csv`).

## Details

Wave 2 is where the attrition happens, and it is what makes the data
useful here. Of the 1,980 subjects, 1,444 responded in Wave 2. Exactly
50 nonrespondents were then drawn at random from each condition and
offered \$4.00 rather than the original \$1.00 to participate. Of those
100 subjects, 72 completed the survey. The follow-up sample is small,
and that is the point: because it is a random sample of the
nonrespondents, the outcomes it recovers stand in for the outcomes of
every nonrespondent, and the worst-case bounds narrow sharply.

The experiment ran a third condition, a placebo group that read nothing
on the topic. It is not analyzed in the paper and is not shipped here,
so the data are the polarized-versus-moderate contrast the paper reports
and nothing else. `data-raw/levendusky_replication.R` builds them from
the archive.

## References

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206.
[doi:10.1017/pan.2016.6](https://doi.org/10.1017/pan.2016.6)

Levendusky, Matthew, and Neil Malhotra (2016). Does Media Coverage of
Partisan Polarization Affect Political Attitudes? *Political
Communication* 33(2):283-301.

## Examples

``` r
# Table 1: attrition by condition
with(levendusky_replication, table(Z_condition, R1))
#>            R1
#> Z_condition   0   1
#>   Moderate  264 731
#>   Polarized 272 713

# Table 3, column 2
estimator_ds(Y_polarization_w2, Z, R1, Attempt, R2,
             minY = 0, maxY = 6, data = levendusky_replication)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      -0.3417454       0.5718157       0.1134230       0.1053948      -0.5283097 
#>       conf.high 
#>       0.7451748 
```
