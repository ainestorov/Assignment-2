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
YDUM = as.matrix(ydum)
## To calculate the coefficients of the regression of Y on X, I use the standard econometric formula,
## beta = ((X'X)^(-1))X'Y. I create a function using this formula so I can replicate in the future.
OLS_betas <- function(x, y) {
  solve(t(x)%*%x)%*%t(x)%*%y
}
betas <- OLS_betas(X, Y)
print(betas)
## The coefficients in the output are 2.46, 1.22, -0.89, 0.07 for Beta0, Beta1, Beta2, Beta3
## respectively. Beta1, Beta2, and Beta3 values are very close to the values that we used to create
## y in exercise 1.

## I will now calculate the standard errors of the regression by using standard econometric formulas.
## I create a function so that I can again replicate in the future using other data.
OLS_se <- function(x, y, coeffs) {
  ## First, I calculate the sample variance using the estimated coefficients from above.
  var_OLS <- sum(((y - x%*%coeffs)^2)/(nrow(x)-ncol(x)))
  ## Next, I multiply this variance by (X'X)^-1.
  sqrt(diag(var_OLS*chol2inv(chol(t(x)%*%x))))
}
std_errs <- OLS_se(X, Y, betas)
## Finally, I print the betas and standard errors so that I can view them together.
print(cbind(betas, std_errs))
## I then double-check that my answers are correct by using the built-in linear regression formula.
summary(lm(y~x1+x2+x3))

# Calculate the standard errors using bootstrap with 49 and 499 replications
## First, I bind the Y and X matrices into a data matrix for the bootstrap function.
data_mtx <- as.matrix(cbind(Y, X))
## I set up the bootstrap function. My inputs are "data" = data matrix, "reps" = # of replications.
bootstrap <- function(data, reps) {
  ## I create an empty matrix which the for loop will populate with the different values that I get
  ## from each sample. This matrix has 4 rows and # of columns based on the # of replications.
  beta_boot <- as.matrix(cbind(rep(0, reps), rep(0, reps), rep(0, reps), rep(0, reps)))
  ## Finally, I create the for loop that will populate the matrix above.
  for(i in 1:reps) {
  ## In the loop, I tell R to sample the data from the data matrix 10000 times with replacement, which
  ## is the definition of the bootstrap. This means that some data points will repeat in each column.
    boot_sample <- data[sample(nrow(data), 10000, replace = TRUE),]
  ## I then separate back out the X and Y variables into "sample" matrices so that I can calculate the
  ## betas for each replication.
    X_sample <- boot_sample[, 2:5]
    Y_sample <- boot_sample[, 1]
  ## Now, I use the same ((X'X)^(-1))X'Y formula as above to calculate the betas for each replication.
    beta_sample <- OLS_betas(X_sample, Y_sample)
  ## These betas then populate the empty matrix that I created before the for loop.
    beta_boot[i, ] <- c(beta_sample)
  }
  ## Finally, I calculate the standard errors for each X column.
  std_err_boot <- apply(beta_boot, 2, sd)
  print(std_err_boot)
}
bootstrap(data_mtx, 49)
## Now that the function is built, I print the function using my data and 49 replication. The
## values are fairly similar to the values from the standard error calculation from the OLS formulas 
## that I did above. I see that the standard errors from the bootstrap with 49 replications are:
## 0.04362, 0.01765, 0.00274, 0.02326 vs. 0.04040, 0.01715, 0.00287, 0.02166 for Beta0 through Beta3
## respectively. These are close but not exact, which I expact given that the bootstrap is an
## approximation.

## I then run the bootstrap for 499 replications using the function created above. 
bootstrap(data_mtx, 499)
## As before, the bootstrap gets values that are fairly close but not exactly equal to the values from
## the standard OLS calculation, given it is an approximation.

# Exercise 3
# Write a function that returns the likelihood of the probit estimation of ydum on X.
## First, I set up the function that will be calculating the likelihood of the probit estimation of
## YDUM on X. Inputs are the coefficients and the x and y data matrices.
probit_mll <- function(coeffs, x, y) {
  ## I create the linear predictor x*coeffs that will be used in the log-likelihood function.
  epsilon <- x %*% coeffs
  ## Since the probit model corresponds to the case where F(x) is the CDF of a standard normal
  ## distribution function, the probability p = F(x*beta) is distributed in the form of a standard
  ## normal distribution.
  p <- pnorm(epsilon)
  ## To calculate the likelihood, I must then take the product of the probabilities raised to the 
  ## dependant variables, as seen in equation (13) in the slides. This would be equivalent to taking
  ## the sum of the log likelihood, as in equation (14). Since I need to maximize the log likelihood, 
  ## I multiply the sum by negative 1 as the negative minimization is equivalent to the maximization.
  -sum(y*log(p)+(1-y)*log(1-p))
}
probit_mll(betas, X, YDUM)
## After I create and run this function, I get 2418.051 for the likelihood of the probit.

# Implement the steepest ascent optimization algorithm to maximize the above likelihood.
## This algorithm is used as an approximation of the derivative by first calculating (f(x+e)-f(x))/e.
## Since I am not replicating this process again in the assignment, I have not made this a function.
## Thus, I first create a small value for e, and a matrix that has the e along the diagonals, so I
## can easily add this e to my matrix of betas (this will be used in the while loop).
e = 0.0001
e_mtx = as.matrix(cbind(c(e,0,0,0), c(0,e,0,0), c(0,0,e,0), c(0,0,0,e)))
## I also create a value for alpha. This will be used to calculate each new set of betas using the
## formula beta_old - alpha*derivatives = beta_new, with the derivatives calculated using the formula 
## described above.
alpha = 0.00001
## I also create a value for a threshold, so that the loop gets close enough to the maximum likelihood
## value (i.e. the difference between the maximum likelihood of the old betas and new betas is very 
## tiny). Once this difference is lower than the threshold, the loop will stop running and approximate 
## the maximum likelihood.
threshold = 0.00000001
## Finally, I create a few new values so that the loop does not get stuck when it runs for the very 
## first time. These just represent starting points for the loop and are not important in any way.
betas_new <- rbind(0, 0, 0, 0)
mll_new = 0
mll_old = 1
## As mentioned above, the loop will keep running until the difference between the maximum likelihoods
## is smaller than the threshold. This is how I set up the first line of the loop.
while (abs(mll_new - mll_old) > threshold) {
  ## I then move the new maximum likelihood created through the formulaic approach below into an old
  ## value, so that I can run through the loop and compare the new value created through the loop to 
  ## the last value that was moved into an old value.
  mll_old <- mll_new
  ## I bind a 4x4 matrix that has each of the new betas calculated in the loop in each column.
  betas_old_mtx <- as.matrix(cbind(betas_new, betas_new, betas_new, betas_new))
  ## I calculate each derivative, as described above, by taking the likelihood at the various points.
  d1 <- (probit_mll(betas_old_mtx[,1]+e_mtx[,1], X, YDUM)-probit_mll(betas_old_mtx[,1], X, YDUM))/e
  d2 <- (probit_mll(betas_old_mtx[,2]+e_mtx[,2], X, YDUM)-probit_mll(betas_old_mtx[,2], X, YDUM))/e
  d3 <- (probit_mll(betas_old_mtx[,3]+e_mtx[,3], X, YDUM)-probit_mll(betas_old_mtx[,3], X, YDUM))/e
  d4 <- (probit_mll(betas_old_mtx[,4]+e_mtx[,4], X, YDUM)-probit_mll(betas_old_mtx[,4], X, YDUM))/e
  ## I bind the derivatives calculated above into a vector so that I can calculate the steep ascent
  ## using the formula beta_old - alpha*derivatives = beta_new, and re-label this as beta_new.
  derivs <- rbind(d1, d2, d3, d4)
  betas_new <- betas_new - alpha*derivs
  ## I find the likelihood value based on the function I created in the previous step of this problem
  ## using these new betas. This value will then come back to the beginning of the loop, be compared to 
  ## the last likelihood value, and then run through the loop again until the difference in absolute
  ## value is less than the threshold.
  mll_new <- probit_mll(betas_new, X, YDUM)
}
print(mll_new)
## After all of these steps, I get a maximum likelihood value of 2175.371.

# How different are the parameters from the true parameters?
print(betas_new)
## The parameters created through the steepest ascent process are very similar to the true ones.
## To compare, I get 2.9206, 1.2126, -0.8920, and 0.0342 for Beta0 to Beta3 respectively, while the
## parameters calculated in exercise 2 are 2.4619, 1.2202, -0.8986, and 0.0741 and the true parameters
## based on the Y variable we created in exercise 1 are Beta1=1.2, Beta2=-0.9, Beta3=0.1. The most 
## glaring difference is for Beta3 which is about 34% of its true value, likely affected by x3 being
## a binary dummy variable. 

# Exercise 4
# Write and optimize the probit, logit, and linear probability models (can use pre-programmed packages)
## As I created the likelihood function in exercise 3, I now only need the gradient for the "BFGS"
## method in order to use the optim function and get the probit coefficients. As before, I begin by
## creating a function to calculate this gradient with the same inputs as before.
probit_gr <- function(coeffs, x, y) {
  ## Similar to before, I also need a linear predictor x*coeffs.
  epsilon <- x %*% coeffs
  ## As before, the probit is where F(x) is the CDF of a standard normal distribution, so the
  ## probability is also distributed as a standard normal.
  p <- pnorm(epsilon)
  ## In order to take the first order derivative of the likelihood function, I must use the chain
  ## rule and the assumed standard normal distribution of F(x*beta).
  foc <- dnorm(epsilon)*(y-p)/(p*(1-p))
  ## Finally, I can obtain the gradient by calculating the cross-product of the matrix x and the first
  ## order conditions above, matching equation (15) in the slides.
  -crossprod(x, foc)
}
## Finally, using the optim function I can calculate the probit model coefficients and save them for
## future reference.
probit_coeffs <- optim(betas, probit_mll, x=X, y=YDUM, gr=probit_gr)$par
print(probit_coeffs)
## I then double-check that my answers are correct by using the built-in generalized linear regression 
## formula for the probit model.
glm_probit <- glm(YDUM~x1+x2+x3, family = binomial(probit))
summary(glm_probit)

## The steps are very similar for the logit. First, I create a function for the logit log likelihood.
logit_mll <- function(coeffs, x, y) {
  ## As before, I need a linear predictor x*coeffs.
  epsilon <- x %*% coeffs
  ## The probability distribution differs for the logit. In particular, the CDF function F(x) is now
  ## the logistic function defined as exp(x)/(1+exp(x)). So, p is no longer distributed as a standard
  ## normal and I must account for this difference in distribution.
  p <- exp(epsilon)/(1+exp(epsilon))
  ## The final step is the same as the calculation of the likelihood in exercise 3 above.
  -sum(y*log(p)+(1-y)*log(1-p))
}
logit_mll(betas, X, YDUM)
## As I did for probit, I must create a gradient function in order to maximize the log likelihood. The
## steps are very similar below.
logit_gr <- function(coeffs, x, y) {
  ## I create a linear predictor for x*coeffs.
  epsilon <- x %*% coeffs
  ## I calculate the probability, accounting for the different CDF function in the logit model.
  p <- exp(epsilon)/(1+exp(epsilon))
  ## I use the chain rule to calculate the first order derivatives.
  foc <- dnorm(exp(epsilon)/(1+exp(epsilon)))*(y-p)/(p*(1-p))
  ## I calculate the cross-product, matching equation (15) in the slides.
  -crossprod(x, foc)
}
## Finally, using the optim function I can calculate the logit model coefficients and save them for
## future reference.
logit_coeffs <- optim(betas, logit_mll, x=X, y=YDUM, gr=logit_gr)$par
print(logit_coeffs)
## I then double-check that my answers are correct by using the built-in generalized linear regression 
## formula for the logit model
glm_logit <- glm(YDUM~x1+x2+x3, family = binomial(logit))
summary(glm_logit)

## Finally, the linear probability model uses the function from exercise 2, replacing Y with YDUM.
lm_coeffs <- OLS_betas(X, YDUM)
print(lm_coeffs)
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

## Below, I calculate the standard errors for the probit and logit models by building on the bootstrap
## function created in exercise 2. The function is almost the same, but has the difference in that it 
## uses optim to calculate the coefficients within the loop rather than the OLS formula. 
data_dum <- as.matrix(cbind(YDUM, X))
## I create a new data set binding the YDUM and X matrices into one to be used in the loop. The new
## bootstrap function takes the same inputs as before, but has additional inputs for coefficients
## (probit or logit coefficients) and type of maximum likelihood/gradient function (probit or logit).
bootstrap_dc <- function(data, reps, coeffs, type_mll, type_gr) {
  beta_boot = as.matrix(cbind(rep(0, reps), rep(0, reps), rep(0, reps), rep(0, reps)))
  for(i in 1:reps) {
    boot_sample = data[sample(nrow(data), 10000, replace = TRUE),]
    X_sample <- boot_sample[, 2:5]
    YDUM_sample <- boot_sample[, 1]
    ## Now I run optim and populate my empty coefficient matrix with the betas from each replication
    ## of the bootstrap. As an aside, if I were building this bootstrap at the beginning of the
    ## entire assignment, I would have created one bootstrap function that had an "if else" statement 
    ## with an additional input that would tell the function whether to run the OLS calculation for
    ## coefficients or the discrete choice calculation here.
    coeffs_sample <- optim(coeffs, type_mll, x=X_sample, y=YDUM_sample, gr=type_gr)$par
    beta_boot[i, ] <- c(coeffs_sample)
    }
  std_err <- apply(beta_boot, 2, sd)
}
## I use this function to print the standard errors for each column of betas using probit with 49 reps.
probit_se <- bootstrap_dc(data_dum, 49, betas, probit_mll, probit_gr)
print(probit_se)

## I also calculate the p-values of these coefficients, so I can compare their statistical significance
## later on. I store these p-values.
pvalues_probit <- pt(abs(probit_coeffs/probit_se), df = Inf)

## I also use the function above to print the standard errors for each column of betas using logit 
## with 49 reps, and store the p-values.
logit_se <- bootstrap_dc(data_dum, 49, betas, logit_mll, logit_gr)
print(logit_se)
pvalues_logit <- pt(abs(logit_coeffs/logit_se), df = Inf)

## Finally, for the linear probability model, I can just calculate the standard errors using the OLS
## function that I developed in Exercise 2. I replace Y with YDUM.
lm_se <- OLS_se(X, YDUM, lm_coeffs)
print(lm_se)
## I store the p-values of these standard errors as well.
pvalues_lm <- pt(abs(lm_coeffs/lm_se), df = Inf)
## I bind the three p-value vectors from above into one table so that I can compare them. In all 
## models, the intercept and betas for x1, x2 are extremely statistically significant. However, the 
## betas for x3 aren't statistically significant for any models at standard significance levels.
pvalue_compare <- cbind(pvalues_probit, pvalues_logit, pvalues_lm)
colnames(pvalue_compare) <- c("probit", "logit", "lm")
print(pvalue_compare)

# Exercise 5
# Compute the marginal effect of X on Y according to the probit and logit models
## The marginal effects are calculated using equation (17) in the slides and multiplying the coeffs 
## from each model by the first derivative of X*coeffs. I create a function here that takes as inputs 
## the X data, the coefficients from the specific model (probit or logit), and uses the distribution 
## depending on if I want to see the 'probit' or 'logit' marginal effects.
marg_effect <- function(x, coeffs, type) {
  if (type == 'probit') {
    dnorm(colMeans(x)%*%coeffs)*c(coeffs)
  } else {
    dlogis(colMeans(x)%*%coeffs)*c(coeffs)
  }
}

## I then simply run my probit coefficients calculated in exercise 4 through the function to find
## the marginal effects.
print(marg_effect(X, probit_coeffs, 'probit'))

## Similarly, I run my logit coefficients through the function to find these marginal effects.
print(marg_effect(X, logit_coeffs, 'logit'))

# Compute the standard deviations using the delta method
## I first download a package that easily calculates the Jacobian matrix used in the function below.
install.packages("numDeriv")
library(numDeriv)
## The delta method entails calculating a Jacobian matrix (matrix of first order partial derivatives)
## using the probit and logit coefficients that have already been calculated and the marginal effect.
## I create a function that takes as inputs the coefficients and model I'm interested in (probit or
## logit) and the variance-covariance matrix from the glm of the model I am interested in.
delta_method <- function(coeffs, type, glm) {
  jacobian_calc <- function(coeffs) {
    marg_effect(X, coeffs, type)
  }
  ## The jacobian function from the package only allows for one input, so I must make a function that 
  ## will be used to pull the correct marginal effects below.
  jacobian_mtx <- jacobian(jacobian_calc, coeffs)
  ## Using the package.
  std_err_me_mtx <- sqrt(t(jacobian_mtx)%*%vcov(glm)%*%jacobian_mtx)
  ## Once the Jacobian is calculated, we can get the standard deviations by simply pre-multiplying by
  ## the transpose of the 4x4 Jacobian and post-multiplying by the 4x4 Jacobian, with the middle 
  ## matrix being the 4x4 variance-covariance matrix from the glm formula in R. 
  std_err_me <- rbind(std_err_me_mtx[1,1], std_err_me_mtx[2,2],std_err_me_mtx[3,3], 
                       std_err_me_mtx[4,4])
  ## The standard deviations of the marginal effects show up in the diagonals of the new 4x4 matrix, 
  ## so I pull these out into a vector and display the SDs calculated through the delta method.
  print(std_err_me)
}
delta_method(probit_coeffs, 'probit', glm_probit)

## I now run the delta method function for the logit.
delta_method(logit_coeffs, 'logit', glm_logit)

# Compute the standard deviations using Bootstrap
## The bootstrap function is almost the exact same as in exercise 4. There is an added input "type" 
## which references the "type" input in the marginal effects function. 
bootstrap_me <- function(data, reps, coeffs, type_mll, type_gr, type) {
  beta_boot = as.matrix(cbind(rep(0, reps), rep(0, reps), rep(0, reps), rep(0, reps)))
  for(i in 1:reps) {
    boot_sample = data[sample(nrow(data), 10000, replace = TRUE),]
    X_sample <- boot_sample[, 2:5]
    YDUM_sample <- boot_sample[, 1]
    coeffs_sample <- optim(betas, type_mll, x=X_sample, y=YDUM_sample, gr=type_gr)$par
    me_sample <- marg_effect(X_sample, coeffs_sample, type)
    ## I add one line to the function to calculate the marginal effects of each sample, and then take 
    ## the standard deviation of these samples.
    beta_boot[i, ] <- c(me_sample)
  }
  sd_me <- apply(beta_boot, 2, sd)
  print(sd_me)
}
bootstrap_me(data_dum, 49, probit_coeffs, probit_mll, probit_gr, 'probit')
## I run the function for the probit model. The standard deviations are similar to those calculated 
## through the delta method, and both methods are an approximation method for the true values.

## I run the function for the logit model.
bootstrap_me(data_dum, 49, logit_coeffs, logit_mll, logit_gr, 'logit')

## END.