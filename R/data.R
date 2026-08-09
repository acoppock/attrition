#' Perceived polarization under double sampling
#'
#' A two-wave survey experiment on Amazon Mechanical Turk, replicating
#' Levendusky and Malhotra (2016) and reported as the application in Coppock,
#' Gerber, Green, and Kern (2017). Subjects read a news article describing the
#' electorate either as sharply divided (the polarized condition) or as focused
#' on common ground (the moderate condition), or read nothing on the topic (the
#' placebo condition). Outcomes were measured immediately in Wave 1 and again in
#' a Wave 2 survey ten days later.
#'
#' Wave 2 is where the attrition happens, and it is what makes the data useful
#' here. Of the 1,980 subjects in the polarized and moderate conditions, 1,444
#' responded in Wave 2. Exactly 50 nonrespondents were then drawn at random from
#' each condition and offered $4.00 rather than the original $1.00 to
#' participate. Of those 100 subjects, 72 completed the survey. The follow-up
#' sample is small, and that is the point: because it is a random sample of the
#' nonrespondents, the outcomes it recovers stand in for the outcomes of every
#' nonrespondent, and the worst-case bounds narrow sharply.
#'
#' The placebo condition is not analyzed in the paper. Every published result
#' uses \code{subset(levendusky, !is.na(Z1))}, the polarized-versus-moderate
#' contrast.
#'
#' @format A data frame with 2,955 rows and 12 columns:
#' \describe{
#'   \item{Z_lev}{Treatment assignment: \code{Placebo}, \code{Moderate}, or \code{Polarized}.}
#'   \item{Z1}{Polarized (1) versus moderate (0); \code{NA} in the placebo condition.
#'     The contrast analyzed in the paper.}
#'   \item{Z2}{Polarized (1) versus placebo (0); \code{NA} in the moderate condition.}
#'   \item{Z3}{Moderate (1) versus placebo (0); \code{NA} in the polarized condition.}
#'   \item{R1}{Responded in the Wave 2 initial sample.}
#'   \item{Attempt}{Selected for the follow-up sample and offered the higher incentive.}
#'   \item{R2}{Responded to the follow-up attempt.}
#'   \item{pid_3_recoded}{Party identification: \code{Dem}, \code{Ind}, or \code{Rep}.
#'     The poststratification variable used in Table 3.}
#'   \item{L_dif}{Perceived polarization at Wave 1, from 0 to 6.}
#'   \item{L_dif_w2}{Perceived polarization at Wave 2, from 0 to 6. The dependent
#'     variable throughout, and the one with missing values.}
#'   \item{L_ex}{Perceived extremity at Wave 1, from 0 to 3.}
#'   \item{L_ex_w2}{Perceived extremity at Wave 2, from 0 to 3.}
#' }
#'
#' Perceived polarization is built from a battery of policy questions. Subjects
#' gave their own view and then guessed how a typical Democratic voter and a
#' typical Republican voter would answer. The outcome is the average absolute
#' difference between the two guesses.
#'
#' @source Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L.
#'   Kern (2016). Replication Data for: Combining double sampling and bounds to
#'   address non-ignorable missing outcomes in randomized experiments. Harvard
#'   Dataverse. \doi{10.7910/DVN/AQB4MP}
#'
#' @references
#' Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern (2017).
#' Combining Double Sampling and Bounds to Address Nonignorable Missing Outcomes
#' in Randomized Experiments. \emph{Political Analysis} 25(2):188-206.
#' \doi{10.1017/pan.2016.6}
#'
#' Levendusky, Matthew, and Neil Malhotra (2016). Does Media Coverage of Partisan
#' Polarization Affect Political Attitudes? \emph{Political Communication}
#' 33(2):283-301.
#'
#' @examples
#' dat <- subset(levendusky, !is.na(Z1))
#'
#' # Table 1: attrition by condition
#' with(dat, table(Z_lev, R1))
#'
#' # Table 3, column 2
#' estimator_ds(L_dif_w2, Z1, R1, Attempt, R2, minY = 0, maxY = 6, data = dat)
"levendusky"
