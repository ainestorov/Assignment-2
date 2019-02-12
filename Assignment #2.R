# Alex Nestorov
# Econ 613
# Assignment #2

# Exercise 1
# First set a seed, and then construct the four vectors of 10K draws from uniform (range = 1:3), 
# gamma (shape = 3, scale = 2), binomial (probability = 0.3), and normal distribution (mean = 2, 
# sd = 1).
set.seed(15)
x1 <- c(runif(10000, 1, 3))
x2 <- c(rgamma(10000, 3, 0.5))
x3 <- c(rbinom(10000, 1, 0.3))
eps <- c(rnorm(10000, 2 ,1))

# Create the variable Y using the given function.
y <- 0.5 + 1.2*x1 + (-0.9)*x2 + 0.1*x3 + eps

# Create the variable ydum where ydum = 1 if Y > mean(Y) and 0 otherwise.
ydum <- as.numeric(y > mean(y))

# Exercise 2
# Calculate the correlation between Y and X1. 
cor(y, x1)
## The correlation between Y and X1 is 0.214395.

# How different is it from 1.2?
cor(y, x1)/1.2
## The correlation between Y and X1 is only about 18% of 1.2. The correlation should definitely be
## between -1 and 1, so the value of 0.214395 makes sense. 1.2 would not make sense. 

# Calculate the coefficients on the regression of Y on X where X = (1, X1, X2, X3)
## First, I have created a column vector of 1's for the first column of the X matrix. 
ones <- rep(1, 10000)
## Next, I create the X matrix by binding the vector of 1's with each individual x vector.
X = as.matrix(cbind(ones, x1, x2, x3))
## Finally, I create the y matrices, one for y and one for ydum.
Y = as.matrix(y)
YDUM <- as.matrix(ydum)
## In order to calculate the coefficients of the regression of Y on X, I use the standard econometric
## formula, beta = ((X'X)^(-1))X'Y.
betas <- solve(t(X)%*%X)%*%t(X)%*%Y
print(betas)

## I will now calculate the standard errors of the regression by using standard econometric formulas.
## First, I calculate the sample variance using the estimated coefficients from above.
var_OLS <- sum(((Y - X%*%betas)^2)/(nrow(X)-ncol(X)))
## Next, I multiply this variance by (X'X)^-1.
std_err_OLS <- sqrt(diag(var_OLS*chol2inv(chol(t(X)%*%X))))
## Finally, I print the betas and standard errors so that I can view them together.
print(cbind(betas, std_err_OLS))
## I then double-check that my answers are correct by using the built-in linear regression formula.
summary(lm(y~x1+x2+x3))

# Calculate the standard errors using bootstrap with 49 and 499 replications
## First, I bind the Y and X matrices into a data matrix so that I can more easily run the bootstrap
## loop.
data <- as.matrix(cbind(Y, X))
## Then, I set up the bootstrap loop. First, I set B = the number of replications.
B = 49
## Next, I create an empty matrix which the loop will populate with the different values that I get
## from the bootstrap. This matrix has 4 rows and B columns based on the number of replications.
beta_boot_49 = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
## Finally, I create the for loop that will populate the matrix above.
for(i in 1:B) {
  ## In the loop, I tell R to sample the data from the data matrix 10000 times with replacement, which
  ## is the definition of the bootstrap. This means that some data points will repeat in each column.
  boot_sample = data[sample(nrow(data), 10000, replace = TRUE),]
  ## I then separate back out the X and Y variables into "sample" matrices so that I can calculate the
  ## betas for each of the 49 bootstrap replications.
  X_sample <- boot_sample[, 2:5]
  Y_sample <- boot_sample[, 1]
  ## Now, I use the same ((X'X)^(-1))X'Y formula as above to calculate the betas for each replication.
  beta_sample = solve(t(X_sample)%*%X_sample)%*%t(X_sample)%*%Y_sample
  ## These betas then populate the empty matrix that I created before the for loop.
  beta_boot_49[i, ] <- c(beta_sample)
}
## Finally, I calculate the standard errors for each X column and print them to compare. Note that the
## values are fairly similar to the values from the standard error calculation from the OLS formulas 
## that I did above. I see that the standard errors from the bootstrap with 49 replications are:
## 0.0404, 0.01715, 0.00287, 0.02166 vs. 0.04362, 0.01765, 0.00274, 0.02326 for beta_0 through beta_3
## respectively. These are close but not exact, which I expact given that the bootstrap is an
## approximation.
std_err_B49 <- apply(beta_boot_49, 2, sd)
print(std_err_B49)

## I then run the bootstrap for 499 replications. I use the exact same methodology as described above,
## updating my B value to reflect 499 replications.
B = 499
beta_boot_499 = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
for(i in 1:B) {
  boot_sample = data[sample(nrow(data), 10000, replace = TRUE),]
  X_sample <- boot_sample[, 2:5]
  Y_sample <- boot_sample[, 1]
  beta_sample = solve(t(X_sample)%*%X_sample)%*%t(X_sample)%*%Y_sample
  beta_boot_499[i, ] <- c(beta_sample)
}
## As before, the bootstrap gets values that are fairly close but not exactly equal to the values from
## the standard OLS calculation, given it is an approximation.
std_err_B499 <- apply(beta_boot_499, 2, sd)
print(std_err_B499)

# Exercise 3
# Write a function that returns the likelihood of the probit estimation of ydum on X.
## First, I set up the function that will be calculating the likelihood of the probit estimation of
## ydum on X.
probit_mll <- function(probit_betas) {
  ## First, I create the linear predictor X*beta that will be used in the log-likelihood function.
  epsilon <- X %*% probit_betas
  ## Since the probit model corresponds to the case where F(x) is the CDF of a standard normal
  ## distribution function, the probability p = F(X*beta) is distributed in the form of a standard
  ## normal distribution.
  p <- pnorm(epsilon)
  ## To calculate the likelihood, I must then take the product of the probabilities raised to the 
  ## dependant variables, as seen in equation (13) in the slides. This would be equivalent to taking
  ## the sum of the log likelihood, as in equation (14). Since I are maximizing the log likelihood, 
  ## I multiply the sum by negative 1 as the negative minimization is equivalent to the maximization.
  -sum(YDUM*log(p)+(1-YDUM)*log(1-p))
}
probit_mll(betas)
## After I create and run this function, I get 2418.051 for the likelihood of the probit.

# Implement the steepest ascent optimization algorithm to maximize the above likelihood.
## This algorithm is used as an approximation of the derivative by first calculating (f(x+e)-f(x))/e.
## Thus, I first create a small value for e and a matrix that has the e along the diagonals, so I
## can easily add this e to my matrix of betas (see more information below).
e = 0.0001
e_mtx = as.matrix(cbind(c(e,0,0,0), c(0,e,0,0), c(0,0,e,0), c(0,0,0,e)))
## I also create a value for alpha. This will be used to calculate each new set of betas using the
## formula beta_0 - alpha*derivatives = beta_1, derivatives are calculating as described above.
alpha = 0.00001
## I also create a value for a threshold, so that the loop gets close enough to the maximum likelihood
## value (i.e. the difference between the maximum likelihood of the old betas and new betas is very 
## tiny) it will stop running and approximate the maximum likelihood that way.
threshold = 0.00000001
## Finally, I create a few new values so that the loop does not get stuck when it runs for the very 
## first time. These represent just starting points for the loop and are not important in any way.
betas_new <- rbind(0, 0, 0, 0)
mll_new = 0
mll_old = 1
## As mentioned above, the loop will keep running until the difference between the maximum likelihoods
## is smaller than the threshold. This is how I set up the first line of the loop.
while (abs(mll_new - mll_old) > threshold) {
  ## I then move the new maximum likelihood created through the formulaic approach below into an old
  ## value, so that I can run through the loop and compare the new value created through the loop to 
  ## the last value created.
  mll_old <- mll_new
  ## I bind a 4x4 matrix that has each of the new betas calculated in the loop in each column.
  betas_old_mtx <- as.matrix(cbind(betas_new, betas_new, betas_new, betas_new))
  ## I calculate each derivative, as described above, by taking the likelihood at the various points.
  d1 <- (probit_mll(betas_old_mtx[,1]+e_mtx[,1])-probit_mll(betas_old_mtx[,1]))/e
  d2 <- (probit_mll(betas_old_mtx[,2]+e_mtx[,2])-probit_mll(betas_old_mtx[,2]))/e
  d3 <- (probit_mll(betas_old_mtx[,3]+e_mtx[,3])-probit_mll(betas_old_mtx[,3]))/e
  d4 <- (probit_mll(betas_old_mtx[,4]+e_mtx[,4])-probit_mll(betas_old_mtx[,4]))/e
  ## I bind the derivatives calculated above into a vector so that I can calculate the steep ascent
  ## formula beta_old - alpha*derivatives = beta_new.
  derivs <- rbind(d1, d2, d3, d4)
  ## I run the calculation mentioned above.
  betas_new <- betas_new - alpha*derivs
  ## I find the likelihood value based on the formula I created in the previous step of this problem
  ## of these new betas. This value will then run back to the beginning of the loop, be compared to 
  ## the last likelihood value, and then run through the loop again until the difference in absolute
  ## value is less than the threshold.
  mll_new <- probit_mll(betas_new)
}
print(mll_new)
## After all those steps, I get a maximum likelihood value of 2175.371.

# How different are the parameters from the true parameters?
print(betas_new)
## The parameters created through the steepest ascent process are fairly similar to the true ones.
## To compare, I get 2.9206, 1.2126, -0.8920, and 0.0342 for beta_0 through beta_3 respectively, while
## the true parameters calculated in Exercise 2 are 2.4619, 1.2202, -0.8986, and 0.0741. The most 
## glaring difference is for beta_3 which is about 50% of its true value, likely caused by x3 being
## a binary dummy variable. 

# Exercise 4
# Write and optimize the probit, logit, and linear probability models (can use pre-programmed packages)
## As I created the likelihood function in exercise 3, I now only need the gradient for the "BFGS"
## method in order to use the optim function and get the probit coefficients. As before, I begin by
## creating a function to calculate this gradient.
probit_gr <- function(probit_betas) {
  ## Similar to before, I also need a linear predictor X*beta.
  epsilon <- X %*% probit_betas
  ## As before, the probit is where F(x) is the CDF of a standard normal distribution, so the
  ## probability is also distributed as a standard normal.
  p <- pnorm(epsilon)
  ## In order to take the first order derivative of the likelihood function, I must use the chain
  ## rule and the assumed standard normal distribution of F(X*beta).
  foc <- dnorm(epsilon)*(YDUM-p)/(p*(1-p))
  ## Finally, I calculate the gradient by calculating the cross-product of the matrix X and the first
  ## order conditions as calculated above, matching equation (15) in the slides.
  -crossprod(X, foc)
}
## Finally, using the optim function I can calculate the probit model coefficients and save them for
## future reference.
probit_coeffs <- optim(betas, probit_mll, probit_gr)$par
## I then double-check that my answers are correct by using the built-in generalized linear regression 
## formula for the probit model.
glm_probit <- glm(YDUM~x1+x2+x3, family = binomial(probit))
summary(glm_probit)

## The steps are very similar for the logit. First, I create a function for the logit log likelihood.
logit_mll <- function(logit_betas) {
  ## As before, I need a linear predictor X*beta.
  epsilon <- X %*% logit_betas
  ## The probability distribution differs for the logit. In particular, the CDF function F(x) is now
  ## the logistic function defined as exp(x)/(1+exp(x)). So, p is no longer distributed as a standard
  ## normal and I must account for this difference in distribution.
  p <- exp(epsilon)/(1+exp(epsilon))
  ## The final step is the same as the calculation of the likelihood in exercise 3 above.
  -sum(YDUM*log(p)+(1-YDUM)*log(1-p))
}
logit_mll(betas)
## As I did for probit, I must create a gradient function in order to maximize the log likelihood. The
## steps are very similar below.
logit_gr <- function(logit_betas) {
  ## I create a linear predictor for X*beta.
  epsilon <- X %*% logit_betas
  ## I calculate the probability, accounting for the different CDF function in the logit model.
  p <- exp(epsilon)/(1+exp(epsilon))
  ## I use the chain rule to calculate the first order derivatives.
  foc <- dnorm(exp(epsilon)/(1+exp(epsilon)))*(YDUM-p)/(p*(1-p))
  ## I calculate the cross-product, matching equation (15) in the slides.
  -crossprod(X, foc)
}
## Finally, using the optim function I can calculate the logit model coefficients and save them for
## future reference.
logit_coeffs <- optim(betas, logit_mll, logit_gr)$par
## I then double-check that my answers are correct by using the built-in generalized linear regression 
## formula for the logit model
glm_logit <- glm(YDUM~x1+x2+x3, family = binomial(logit))
summary(glm_logit)

## Finally, the linear probability model follows the same calculation as I did in exercise 2,
## replacing Y with YDUM.
lm_coeffs <- solve(t(X)%*%X)%*%t(X)%*%YDUM
## I then double-check that my answers are correct by using the built-in linear regression formula.
summary(lm(YDUM~x1+x2+x3))

# Interpret and compare the estimated coefficients. How significant are they?
## I first create a data table where I can easily view and compare the coefficients from the three
## different models of regressing X on YDUM.
compare <- cbind(probit_coeffs, logit_coeffs, lm_coeffs)
colnames(compare) <- c("probit", "logit", "lm")
print(compare)
## In all three cases, the coefficients have the same directional interpretation but differ in
## magnitude. In all three, x1 and x3 have a positive relationship with YDUM while x2 has a negative
## relationship with YDUM. For the probit and logit models, one cannot directly interpret the 
## coefficient more than directionally given that I have not yet calculated the marginal effects. For
## the linear probability model, the coefficients tell us how much YDUM is expected to increase when
## each of the x's increases by one unit, holding the others constant. The issue with using the linear
## probability model in discrete choice is that the OLS fitted values may not be between 0 and 1, which
## would not make any sense (although this is not the case here).

## Below, I calculate the standard errors for the probit and logit models using the bootstrap method
## that was used in Exercise 2. The methodology is the exact same, except I create a few new variables
## to reflect the different model being analyzed. I run each model with 49 replications.
data_dum <- as.matrix(cbind(YDUM, X))
B = 49
beta_boot_49_p = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
for(i in 1:B) {
  boot_sample = data_dum[sample(nrow(data_dum), 10000, replace = TRUE),]
  X_sample <- boot_sample[, 2:5]
  YDUM_sample <- boot_sample[, 1]
  ## Since I used the maximum likelihood to calculate the coefficients for the probit and logit models
  ## I have to run the maximum likelihood formulas through the bootstrap loop as I have done below.
  ## As before, I have created sample maximum likelihood and gradient functions to optimize.
  probit_mll_sample <- function(probit_betas) {
    epsilon <- X_sample %*% probit_betas
    p <- pnorm(epsilon)
    -sum(YDUM_sample*log(p)+(1-YDUM_sample)*log(1-p))
  }
  probit_gr_sample <- function(probit_betas) {
    epsilon <- X_sample %*% probit_betas
    p <- pnorm(epsilon)
    foc <- dnorm(epsilon)*(YDUM_sample-p)/(p*(1-p))
    -crossprod(X_sample, foc)
  }
  ## Now that I have created these functions, I can run optim and populate my empty probit coefficient
  ## matrix with the betas from each replication of the bootstrap.
  probit_coeffs_sample <- optim(betas, probit_mll_sample, probit_gr_sample)$par
  beta_boot_49_p[i, ] <- c(probit_coeffs_sample)
}
## Finally, I calcualte the standard errors for each column of betas.
std_err_49_p <- apply(beta_boot_49_p, 2, sd)
## I also calculate the p-values of these coefficients, so I can compare their statistical significance
## later on. I store these p-values.
pvalues_probit <- pt(abs(probit_coeffs/std_err_49_p), df = Inf)

## The bootstrapping process is the exact same for the logit, using the maximum likelihood and gradient
## formulas that I developed for the logit above, calculating the standard errors for each column of
## betas, and then calculating and storing the p-values to view statistical significance.
beta_boot_49_l = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
for(i in 1:B) {
  boot_sample = data_dum[sample(nrow(data_dum), 10000, replace = TRUE),]
  X_sample <- boot_sample[, 2:5]
  YDUM_sample <- boot_sample[, 1]
  logit_mll_sample <- function(logit_betas) {
    epsilon <- X_sample %*% logit_betas
    p <- exp(epsilon)/(1+exp(epsilon))
    -sum(YDUM_sample*log(p)+(1-YDUM_sample)*log(1-p))
  }
  logit_gr_sample <- function(logit_betas) {
    epsilon <- X_sample %*% logit_betas
    p <- exp(epsilon)/(1+exp(epsilon))
    foc <- dnorm(exp(epsilon)/(1+exp(epsilon)))*(YDUM_sample-p)/(p*(1-p))
    -crossprod(X_sample, foc)
  }
  logit_coeffs_sample <- optim(betas, logit_mll_sample, logit_gr_sample)$par
  beta_boot_49_l[i, ] <- c(logit_coeffs_sample)
}
std_err_49_l <- apply(beta_boot_49_l, 2, sd)
pvalues_logit <- pt(abs(logit_coeffs/std_err_49_l), df = Inf)

## Finally, for the linear probability model, I can just calculate the standard errors using the OLS
## formula that I developed in Exercise 2. I replace Y with YDUM and run the same formula for standard
## errors.
sigma_squared_ydum <- sum(((YDUM - X%*%lm_coeffs)^2)/(nrow(X)-ncol(X)))
std_errors_OLS_ydum <- sqrt(diag(sigma_squared_ydum*chol2inv(chol(t(X)%*%X))))
## I store the p-values of these standard errors as well.
pvalues_lm <- pt(abs(lm_coeffs/std_errors_OLS_ydum), df = Inf)
## Finally, I bind the three p-value vectors that I calculated above into one table so that I can 
## compare them. In all three models, the intercept and betas for x1 and x2 are extremely statistically
## significant. However, the betas for x3 are not statistically significant for any of the models at 
## any standard statistical significance level.
pvalue_compare <- cbind(pvalues_probit, pvalues_logit, pvalues_lm)
colnames(pvalue_compare) <- c("probit", "logit", "lm")
print(pvalue_compare)

# Exercise 5
# Compute the marginal effect of X on Y according to the probit and logit models
## To compute standard deviations using the Delta Method, I download a package that easily calculates
## the Jacobian matrix that is used in that step. This does not have an impact on the marginal effects
## which are calculated without the package.
install.packages("numDeriv")
library(numDeriv)
## The marginal effects are calculated using equation (17) in the slides and multiplying the 
## coefficients from each model by the first derivative of X*coefficients. I create a function here
## that will aid in the calculation of the Jacobian matrices below.
jacobian_calc_p <- function(probit_results) 
  dnorm(colMeans(X)%*%probit_results)*c(probit_results)
## I then simply run my probit coefficients calculated in Exercise 4 through the function to find
## the marginal effects.
marg_probit <- jacobian_calc_p(probit_coeffs)
print(marg_probit)

## Similarly, I create a function for the Jacobian of the logit coefficients and then run the logit
## coefficients calculated above through the function to get my marginal effects for the logit model.
jacobian_calc_l <- function(logit_results)
  dlogis(colMeans(X)%*%logit_results)*c(logit_results)
marg_logit <- jacobian_calc_l(logit_coeffs)
print(marg_logit)

# Compute the standard deviations using the delta method
## The delta method entails calculating a Jacobian matrix (matrix of first order partial derivatives)
## using the probit and logit coefficients that have already been calculated and the functions that I
## created above.
jacobian_mtx_p <- jacobian(jacobian_calc_p, probit_coeffs)
## Once the Jacobian is calculated, we can get the standard deviations by simply pre-multiplying by
## the transpose of the 4x4 Jacobian and post-multiplying by the 4x4 Jacobian, with the middle matrix 
## being the 4x4 variance-covariance matrix that I pull from the glm formula in R. 
std_err_probit_me_mtx <- sqrt(t(jacobian_mtx_p)%*%vcov(glm_probit)%*%jacobian_mtx_p)
## The standard deviations of the marginal effects then show up in the diagonals of the new 4x4 matrix, 
## so I pull these out into a vector and display the standard deviations calculated through the delta 
## method.
std_err_probit_me <- rbind(std_err_probit_me_mtx[1,1], std_err_probit_me_mtx[2,2],
                           std_err_probit_me_mtx[3,3], std_err_probit_me_mtx[4,4])
print(std_err_probit_me)

## The process for the delta method for the logit model is the exact same. First, I create the Jacobian
## matrix. Then, I pre- and post-multiply the variance-covariance matrix of the logit model by this
## Jacobian, and finally I pull out the diagonals of the resultant matrix to display the standard
## deviations of the marginal effects.
jacobian_mtx_l <- jacobian(jacobian_calc_l, logit_coeffs)
std_err_logit_me_mtx <- sqrt(t(jacobian_mtx_l)%*%vcov(glm_logit)%*%jacobian_mtx_l)
std_err_logit_me <- rbind(std_err_logit_me_mtx[1,1], std_err_logit_me_mtx[2,2],
                          std_err_logit_me_mtx[3,3], std_err_logit_me_mtx[4,4])
print(std_err_logit_me)

# Compute the standard deviations using Bootstrap
## The bootstrap methodology is almost the exact same as that used in Exercise 4 to calculate the 
## standard errors in that case, with a few extra steps.
B = 49
beta_boot_49_p_me = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
for(i in 1:B) {
  boot_sample = data_dum[sample(nrow(data_dum), 10000, replace = TRUE),]
  X_sample <- boot_sample[, 2:5]
  YDUM_sample <- boot_sample[, 1]
  probit_mll_sample <- function(probit_betas) {
    epsilon <- X_sample %*% probit_betas
    p <- pnorm(epsilon)
    -sum(YDUM_sample*log(p)+(1-YDUM_sample)*log(1-p))
  }
  probit_gr_sample <- function(probit_betas) {
    epsilon <- X_sample %*% probit_betas
    p <- pnorm(epsilon)
    foc <- dnorm(epsilon)*(YDUM_sample-p)/(p*(1-p))
    -crossprod(X_sample, foc)
  }
  ## As before, I create an empty matrix to fill through the for loop, create sample betas using 
  ## sample data from my actual data, create maximum likelihood and gradient functions and then 
  ## calculate the probit coefficients as before.
  probit_coeffs_sample <- optim(betas, probit_mll_sample, probit_gr_sample)$par
  ## Now, there is one additional step, in that I have to calculate the marginal effects for each
  ## sample of betas. Since I have created a function that calculates these effects as seen above, I
  ## simply run the betas through this function and create a new row of marginal effects each time
  ## that the loop runs.
  probit_me_sample <- jacobian_calc_p(probit_coeffs_sample)
  beta_boot_49_p_me[i, ] <- c(probit_me_sample)
}
## Finally, I just take the standard deviation of these marginal effects. They are quite similar to
## the standard deviations calculated through the delta method, and both methods serve as an
## approximation method for the true values.
sd_49_me_p <- apply(beta_boot_49_p_me, 2, sd)
print(sd_49_me_p)

## I run the bootstrap in the same exact way for the logit model, using the function I created above
## to calculate the marginal effects of each run through the loop and then take the standard deviation
## of all of these marginal effects.
beta_boot_49_l_me = as.matrix(cbind(rep(0, B), rep(0, B), rep(0, B), rep(0, B)))
for(i in 1:B) {
  boot_sample = data_dum[sample(nrow(data_dum), 10000, replace = TRUE),]
  X_sample <- boot_sample[, 2:5]
  YDUM_sample <- boot_sample[, 1]
  logit_mll_sample <- function(logit_betas) {
    epsilon <- X_sample %*% logit_betas
    p <- exp(epsilon)/(1+exp(epsilon))
    -sum(YDUM_sample*log(p)+(1-YDUM_sample)*log(1-p))
  }
  logit_gr_sample <- function(logit_betas) {
    epsilon <- X_sample %*% logit_betas
    p <- exp(epsilon)/(1+exp(epsilon))
    foc <- dnorm(exp(epsilon)/(1+exp(epsilon)))*(YDUM_sample-p)/(p*(1-p))
    -crossprod(X_sample, foc)
  }
  logit_coeffs_sample <- optim(betas, logit_mll_sample, logit_gr_sample)$par
  logit_me_sample <- jacobian_calc_l(logit_coeffs_sample)
  beta_boot_49_l_me[i, ] <- c(logit_me_sample)
}
sd_49_me_l <- apply(beta_boot_49_l_me, 2, sd)
print(sd_49_me_l)
