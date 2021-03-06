% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_facing_functions.R
\name{estimator_ds}
\alias{estimator_ds}
\title{Extreme Value Bounds with Double Sampling}
\usage{
estimator_ds(Y, Z, R1, Attempt, R2, minY, maxY, strata = NULL, alpha = 0.05,
  data)
}
\arguments{
\item{Y}{The (unquoted) outcome variable. Must be numeric.}

\item{Z}{The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.}

\item{R1}{The (unquoted) initial sample respose indicator variable. Must be numeric and take values 0 or 1.}

\item{Attempt}{The (unquoted) follow-up sample attempt indicator variable. Must be numeric and take values 0 or 1.}

\item{R2}{The (unquoted) follow-up sample respose indicator variable. Must be numeric and take values 0 or 1.}

\item{minY}{The minimum possible value of the outcome (Y) variable.}

\item{maxY}{The maximum possible value of the outcome (Y) variable.}

\item{strata}{A single (unquoted) variable that indicates which strata units are in.}

\item{alpha}{The desired significance level. 0.05 by default.}

\item{data}{A dataframe}
}
\value{
A results matrix
}
\description{
Extreme Value Bounds with Double Sampling
}
\examples{
set.seed(343) # For reproducibility
N <- 1000

# Potential Outcomes
Y_0 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
Y_1 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))

R1_0 <- rbinom(N, 1, prob = 0.7)
R1_1 <- rbinom(N, 1, prob = 0.8)

R2_0 <- rbinom(N, 1, prob = 0.9)
R2_1 <- rbinom(N, 1, prob = 0.95)

# Covariate
strata <- as.numeric(Y_0 > 2)

# Random Assignment
Z <- rbinom(N, 1, .5)

# Reveal Initial Sample Outcomes
R1 <- Z*R1_1 + (1-Z)*R1_0 # Initial sample response
Y_star <- Z*Y_1 + (1-Z)*Y_0 # True outcomes
Y <- Y_star
Y[R1==0] <- NA # Mask outcome of non-responders

# Conduct Double Sampling
Attempt <- rep(0, N)
Attempt[is.na(Y)] <- rbinom(sum(is.na(Y)), 1, .5)

R2 <- rep(0, N)
R2[Attempt==1] <-  (Z*R2_1 + (1-Z)*R2_0)[Attempt==1]

Y[R2==1 & Attempt==1] <- Y_star[R2==1 & Attempt==1]

df <- data.frame(Y, Z, R1, Attempt, R2, strata)

# Without post-stratification
estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, data=df)

# With post-stratification
estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, strata=strata, data=df)

}

