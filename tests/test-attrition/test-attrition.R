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




# test extreme value bounds -----------------------------------------------

rm(list=ls())
N <- 100

Y_0 <- rep(c(0,1), c(50, 50))
Y_1 <- rep(c(0,1), c(50, 50))

R1_0 <- rbinom(N, 1, prob = 0.5)
R1_1 <- rbinom(N, 1, prob = 0.5)


Z <- rbinom(N, 1, .5)

# Reveal Initial Sample Outcomes
R1 <- Z*R1_1 + (1-Z)*R1_0 # Initial sample response
Y_star <- Z*Y_1 + (1-Z)*Y_0 # True outcomes
Y <- Y_star
Y[R1==0] <- NA # Mask outcome of non-responders

test_df <- data.frame(Y, Z, R1)

estimator_ev(Y, Z, R1, 0, 1, data = test_df)


# Manual

Y_upp <- Y
Y_upp[is.na(Y)] <- 1
Y_low <- Y
Y_low[is.na(Y)] <- 0

mean(Y_upp[Z==1]) - mean(Y_low[Z==0])
mean(Y_low[Z==1]) - mean(Y_upp[Z==0])





RT1 <- mean(!is.na(Y[Z==1]))
YT1 <- mean((Y[Z==1]), na.rm = TRUE)

RC1 <- mean(!is.na(Y[Z==0]))
YC1 <- mean((Y[Z==0]), na.rm = TRUE)

Y_L <- 0
Y_U <- 1


a_hat_prime <- RT1*YT1 + (1-RT1)*Y_L
b_hat_prime <- RT1*YT1 + (1-RT1)*Y_U
c_hat_prime <- RC1*YC1 + (1-RC1)*Y_L
d_hat_prime <- RC1*YC1 + (1-RC1)*Y_U

upp_est_prime <- b_hat_prime - c_hat_prime
low_est_prime <- a_hat_prime - d_hat_prime





