####################################################################################################################
# ABC's of M-estimation
#   Code for M-estimator of the mean (Section 1)
#
# Paul Zivich (2025/06/04)
####################################################################################################################

############################################
# Loading dependencies

import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import root
from scipy.optimize import approx_fprime
import delicatessen as deli
from delicatessen import MEstimator

print("versions")
print('NumPy:        ', np.__version__)
print('Pandas:       ', pd.__version__)
print('SciPy:        ', sp.__version__)
print('Delicatessen: ', deli.__version__, '\n')

############################################
# Data

y = np.array([7, 1, 5, 3, 24])   # Data (observations)
n = len(y)                       # Number of observations

############################################
# Mean by usual / closed-form estimator

mu_closed = np.sum(y) / n         # Manually computing mean by sum and dividing by n
# mu_closed = np.mean(y)          # Built-in function in NumPy
print("Estimated mean")
print("Closed-form:", mu_closed)

############################################
# Defining estimating equation


def estimating_function(mu):
    return y - mu


def estimating_equation(mu):
    estf = np.asarray(estimating_function(mu))  # Return estimating function
    return np.sum(estf)                         # Sum over all estimating functions


############################################
# Root-finding

proc = root(estimating_equation,     # Function to find root(s) of
            x0=np.array([0, ]),      # ... starting values for root-finding procedure
            method='lm')             # ... algorithm to use (Levenberg-Marquardt here)
mu_root = proc["x"][0]
print("Root-finder:", mu_root)

############################################
# Baking the bread (approximate derivative)

deriv = approx_fprime(xk=proc["x"],            # Array of values to compute derivative at (root of estimating equation)
                      f=estimating_equation,   # ... function to find derivative of
                      epsilon=1e-9)            # ... distance of points for numerical approximation (should be small)
bread = -1*deriv / n

############################################
# Cooking the filling (matrix algebra)

filling = np.sum(estimating_function(mu_root) * estimating_function(mu_root))
filling = filling / n

############################################
# Assembling the sandwich (matrix algebra)

sandwich = (bread**-1) * filling * (bread**-1) / n
se = np.sqrt(sandwich)

print("95% CI:", np.round([mu_root - 1.96*se[0],
                           mu_root + 1.96*se[0]],
                          3))

####################################################################################################################
# Using delicatessen instead of by-hand


def psi(theta):
    return y - theta


mestr = MEstimator(psi, init=[0, ])
mestr.estimate(solver='lm')

print("Deli:   ", mestr.theta)
print("95% CI: ", np.round(mestr.confidence_intervals(), 3))

####################################################################################################################
# EXPECTED OUTPUT (versions may differ):

# versions
# NumPy:         1.25.2
# Pandas:        1.4.1
# SciPy:         1.11.2
# Delicatessen:  3.2
#
# Estimated mean
# Closed-form: 8.0
# Root-finder: 8.0
# 95% CI: [ 0.772 15.228]
# Deli:    [8.]
# 95% CI:  [[ 0.772 15.228]]

####################################################################################################################
# END
