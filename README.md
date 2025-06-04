# ABC's of Estimating Equations (M-estimation)

Materials for the SER 2025 "ABC's of Estimating Equations" workshop (a continuation of the ABC's of M-estimation
workshop) with support for SAS, R, and Python.

Materials and workshop created by: Paul Zivich, Rachael Ross, Bonnie Shook-Sa, Stephen Cole

[![DOI](https://zenodo.org/badge/618596455.svg)](https://zenodo.org/badge/latestdoi/618596455)

-----------

## Before the workshop

Copy the files from this repository to your computer that will be used during the session. This can be done directly
using your GitHub account. Alternatively, you can click the `<> Code` button on the top of this repository, navigate to
the `Local` tab, and click `Download ZIP` to manually download.

In the `code/` folder, find the file corresponding to your preferred software program and open it in your preferred
editor. Then follow the directions below to ensure your computer is properly set up for the workshop:

### SAS
- Open mean.sas and run the entire script.

Output should:
- Generate output in the Results Viewer that consists of a data set with the columns B, SANDWICH, SE, lcl, and ucl.
The values should be 8, 13.6, 3.68782, 0.77188, and 15.2281, respectively.
- No errors should be produced in the log file.

### R
- Install the following libraries:
    - tidyverse, numDeriv, rootSolve, geex
- Open mean.R and run the entire script

Output should:

```
Estimated mean
Closed-form: 8.000
Root-finder: 8.000
95% CI: [ 0.772 15.228]
Geex:    [8.000]
95% CI:  [ 0.772 15.228]
```

- No errors should be produced.

### Python
- Install the following libraries:
    - NumPy, SciPy, pandas, statsmodels, delicatessen
    - The version of SciPy must be at least 1.9.0 for the provided code to function properly
- Open mean.py (or mean.ipynb) and run the entire script

Output should be:

```
Closed-form: 8.0
Root-finder: 8.0
95% CI: [ 0.772 15.228]
Deli:    [8.]
95% CI:  [[ 0.772 15.228]]
```

- No errors should be produced.

If you have problems or encounter errors, please open an issue on the GitHub page or email me at pzivich@unc.edu

## Schedule

8:30 - 8:45 Introductions and start-up

8:45 - 9:45 Introduction to M-estimation

9:45 - 10:00 Break

10:00 - 11:15 Applications

11:15 - 11:30 Break

11:30 - 12:30 Concluding discussion
