###################################################################################################################
# ABC's of M-estimation
# 	Code for M-estimator of the mean (Section 1)
# 
# Rachael Ross (2023/03/23)
###################################################################################################################

############################################
# Loading Libraries 

library("numDeriv")
library("rootSolve")
library("geex")

############################################
# Data 

y <- c(7,1,5,3,24)  # Data (observations)
n <- length(y)      # Number of observations

############################################
# Mean by usual / closed-form estimator

mu_closed <- sum(y)/n   # Manually computing mean by sum and dividing by n
# mu_closed = mean(y)   # Built-in function
print("Estimated mean")
print(paste("Closed-form:", round(mu_closed,3)))

############################################
# Defining estimating equation

estimating_function <- function(mu){
  estf = y - mu
  return(estf) 
}

estimating_equation <- function(mu){
  estf = estimating_function(mu)          # Return estimating function
  este = sum(estf)                        # Sum over all estimating functions
  return(este)
}

############################################
# Root-finding

proc <- rootSolve::multiroot(f = estimating_equation,     # Function to find root(s) of
                             start = c(0))                # ... starting values for root-finding procedure

mu_root <- proc$root
print(paste("Root-finder:", round(mu_root,3))) 

############################################
# Baking the bread (approximate derivative)

deriv <- numDeriv::jacobian(func = estimating_equation,   # Function to find derivative of
                            x = mu_root)                  # Array of values to compute derivative at (root of estimating equation)

bread <- -1*deriv / n

############################################
# Cooking the filling (matrix algebra)

outerprod <- sum(estimating_function(mu_root) * estimating_function(mu_root)) 
# outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root) # alternative code using matrix algebra
filling <- outerprod/n 

############################################
# Assembling the sandwich (matrix algebra)

sandwich <- (bread^-1) %*% filling %*% (bread^-1)
se <- sqrt(sandwich / n)

print(paste("95% CI:", round(mu_root - 1.96*se,3), round(mu_root + 1.96*se,3))) 
            
####################################################################################################################
# Using geex instead of by-hand

psi <- function(data){      # Function of estimating functions (to be used in geex::m_estimate)
  y <- data$y
  function(theta){
    return(y-theta)
  }
}

mestr <- geex::m_estimate(estFUN = psi,                                      # Function of estimating functions
                          data = as.data.frame(y),                           # Data to be used (must be data frame)
                          root_control = setup_root_control(start = c(0)))   # Set starting values

mu_geex <- roots(mestr)         # Extract roots
se_geex <- sqrt(vcov(mestr))    # Extract finite sample variance and take sqrt to get se

print(paste("Geex:", round(mu_geex,3)))
print(paste("95% CI:", round(mu_geex - 1.96*se_geex,3), round(mu_geex + 1.96*se_geex,3)))

      