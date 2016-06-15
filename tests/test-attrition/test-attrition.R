library(attrition)
set.seed(343) # For reproducibility
N <- 1000

# Potential Outcomes
Y_0 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
Y_1 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))

R1_0 <- rbinom(N, 1, prob = 0.7)
R1_1 <- rbinom(N, 1, prob = 0.8)

R2_0 <- rbinom(N, 1, prob = 0.7)
R2_1 <- rbinom(N, 1, prob = 0.75)

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

df <- within(df,{
  Z_rev <- 1-Z
})

# Sensitivity
# Pos
sens <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY=0, maxY=5, data=df)
sens$sensitivity_plot
# Neg
sens <- sensitivity_ds(Y, Z_rev, R1, Attempt, R2, minY=0, maxY=5, data=df)


sens <- sensitivity_ds(Y, Z_rev, R1, Attempt, R2, minY=0, maxY=2, data=df)


debugonce(sensitivity_ds)
