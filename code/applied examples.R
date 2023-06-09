###################################################################################################################
# ABC's of M-estimation
# 	Code for applied examples (Section 2)
# 
# Rachael Ross (2023/03/24)
###################################################################################################################

############################################
# Loading Libraries 

library("tidyverse")
library("numDeriv")
library("rootSolve")
library("geex")

############################################
# EXAMPLE 1: LOGISTIC REGRESSION

############################################
# Data 

dat <- tibble(anemia=c(rep(0,4),rep(1,4)),
              bp=rep(c(rep(0,2),rep(1,2)),2),
              ptb=rep(c(0,1),4),
              n=c(496,74,113,25,85,15,15,3)) %>%
  uncount(n) 

n <- nrow(dat)     # Number of observations

############################################
# Regression by MLE 

reg_mle <- glm(ptb ~ anemia + bp, data=dat, family="binomial")   

ests_mle <-as.data.frame(cbind("beta"=reg_mle$coefficients,
                               "se"=sqrt(diag(vcov(reg_mle))))) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.96*se)
intlist <- row.names(ests_mle)

print("Estimated logistic regression")
print("MLE")
print(round(ests_mle,3))

############################################
# Defining estimating equation

estimating_function <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia + beta[3]*dat$bp)
  ef_1 <- (dat$ptb - p)
  ef_2 <- (dat$ptb - p)*dat$anemia
  ef_3 <- (dat$ptb - p)*dat$bp
  return(cbind(ef_1, ef_2, ef_3)) 
}

estimating_equation <- function(beta){
  estf = estimating_function(beta)        # Return estimating function
  este = colSums(estf)                    # Estimating equations are sum
  return(este)
}

############################################
# Root-finding

proc <- rootSolve::multiroot(f = estimating_equation,     # Function to find root(s) of
                             start = c(-2,0,0))           # Starting values for root-finding procedure

beta_root <- proc$root

############################################
# Baking the bread (approximate derivative)

deriv <- numDeriv::jacobian(func = estimating_equation,   # Function to find derivative of
                            x = beta_root)                # Array of values to compute derivative at (root of estimating equation)

bread <- -1*deriv / n

############################################
# Cooking the filling (matrix algebra)

outerprod <- t(estimating_function(beta_root)) %*% estimating_function(beta_root) # Outer product of the residuals
filling <- outerprod/n 

############################################
# Assembling the sandwich (matrix algebra)

sandwich <- solve(bread) %*% filling %*% t(solve(bread))
se <- sqrt(diag(sandwich / n))

ests_mest <-as.data.frame(cbind("beta"=beta_root,
                                "se"=se)) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.96*se)

row.names(ests_mest) <- intlist
            
print("M-Estimation, by hand")
print(round(ests_mest,3))

####################################################################################################################
# Using geex 

geex_ef <- function(data){               # Function of estimating functions (to be used in geex::m_estimate)
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  function(theta){
    p <- plogis(theta[1] + theta[2]*anemia + theta[3]*bp)
    c((ptb - p),
      (ptb - p)*anemia,
      (ptb - p)*bp)
  }
}

mestr <- geex::m_estimate(estFUN = geex_ef,                                       # Function of estimating functions
                          data = dat,                                             # Data to be used (must be data frame)
                          root_control = setup_root_control(start = c(-2,0,0)))   # Set starting values

beta_geex <- roots(mestr)             # Extract roots
se_geex <- sqrt(diag(vcov(mestr)))    # Extract finite sample variance and take sqrt to get se

ests_geex <-as.data.frame(cbind("beta"=beta_geex,
                                "se"=se_geex)) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.96*se)

row.names(ests_geex) <- intlist

print("M-Estimation, by geex")
print(round(ests_geex,3))

############################################
# Results for logistic regression

print("Estimated logistic regression")
print("MLE")
print(round(ests_mle,3))

print("M-Estimation, by hand")
print(round(ests_mest,3))

print("M-Estimation, by hand")
print(round(ests_geex,3))


################################################
# EXAMPLE 2a: STANDARDIZATION BY G-COMPUTATION

################################################
# Using geex 

geex_ef2A <- function(data){               
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  
  function(theta){
    beta <- theta[1:3]
    mu <- theta[4:5]
    delta <- theta[6:7]
    
    p <- plogis(beta[1] + beta[2]*anemia + beta[3]*bp)
    
    ef_1 <- (ptb - p)
    ef_2 <- (ptb - p)*anemia
    ef_3 <- (ptb - p)*bp
    
    ef_r1 <- plogis(beta[1] + beta[2]*1 + beta[3]*bp) - mu[1]
    ef_r0 <- plogis(beta[1] + beta[2]*0 + beta[3]*bp) - mu[2]
    
    ef_rd <- mu[1] - mu[2] - delta[1]
    ef_lnrr <- log(mu[1]/mu[2]) - delta[2]
    return(c(ef_1,ef_2,ef_3,ef_r1,ef_r0,ef_rd,ef_lnrr))
  }
}

mest_2a <- geex::m_estimate(estFUN = geex_ef2A,                                       
                            data = dat,                                             
                            root_control = setup_root_control(start = c(-2,0,0,.1,.1,0,0)))   

theta2a <- roots(mest_2a)             
se2a <- sqrt(diag(vcov(mest_2a)))    

ests_2a <-as.data.frame(cbind("beta"=theta2a,
                                "se"=se2a)) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.96*se)

row.names(ests_2a) <- c(intlist,"risk1","risk0","rd","lnrr")

print("G-computation by M-Estimation, using geex")
print(round(ests_2a,3))


################################################
# EXAMPLE 2b: STANDARDIZATION BY IPW

################################################
# Using geex 

geex_ef2B <- function(data){               
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  
  function(theta){
    alpha <- theta[1:2]
    mu <- theta[3:4]
    delta <- theta[5:6]
    
    pscore <- plogis(alpha[1] + alpha[2]*bp)
    ef_1 <- (anemia - pscore)
    ef_2 <- (anemia - pscore)*bp
    
    wt <- anemia/pscore + (1-anemia)/(1-pscore)
    ef_r1 <- anemia*wt*ptb - mu[1]
    ef_r0 <- (1 - anemia)*wt*ptb - mu[2]
    
    ef_rd <- mu[1] - mu[2] - delta[1]
    ef_lnrr <- log(mu[1]/mu[2]) - delta[2]
    return(c(ef_1,ef_2,ef_r1,ef_r0,ef_rd,ef_lnrr))
  }
}

mest_2b <- geex::m_estimate(estFUN = geex_ef2B,                                       
                            data = dat,                                             
                            root_control = setup_root_control(start = c(-2,0,.1,.1,0,0)))   

theta2b <- roots(mest_2b)             
se2b <- sqrt(diag(vcov(mest_2b)))    

ests_2b <-as.data.frame(cbind("beta"=theta2b,
                              "se"=se2b)) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.96*se)

row.names(ests_2b) <- c(intlist[1],intlist[3],"risk1","risk0","rd","lnrr")

print("IPW by M-Estimation, using geex")
print(round(ests_2b,3))


############################################
# EXAMPLE 3: DATA FUSION

############################################
# Data 

datfusion <- tibble(r=c(rep(1,950),rep(0,331)),
              y=c(rep(0,950),                           # Missing Y in R=1 set to 0
                  rep(1,242),rep(0,89)),
              w=c(rep(1,680),rep(0,950-680),
                  rep(1,204),rep(0,38),rep(1,18),rep(0,71)))



################################################
# Using geex 

geex_ef3 <- function(data){               
  r <- data$r
  y <- data$y
  w <- data$w
  
  function(theta){
    ef_1 <- r*(w - theta[1])
    ef_2 <- (1 - r)*y*(w - theta[2])
    ef_3 <- (1 - r)*(1 - y)*((1 - w) - theta[3])
    ef_4 <- theta[4]*(theta[2] - (1 - theta[3])) - (theta[1] - (1 - theta[3]))

    return(c(ef_1,ef_2,ef_3,ef_4))
  }
}

mest_3 <- geex::m_estimate(estFUN = geex_ef3,                                       
                           data = datfusion,                                             
                           root_control = setup_root_control(start = c(.7,1,1,0.7)))   

theta3 <- roots(mest_3)             
se3 <- sqrt(diag(vcov(mest_3)))    

ests_3 <-as.data.frame(cbind("beta"=theta3,
                              "se"=se3)) %>%
  mutate(lcl = beta - 1.96*se,
         ucl = beta + 1.95*se)

row.names(ests_3) <- c("pr_w","se","sp","pr_y")

print("ME correction, using geex")
print(round(ests_3,3))
      