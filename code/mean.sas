/*******************************************************************************************************************
ABC's of M-estimation
	Code for M-estimator of the mean (Section 1)

Paul Zivich & Rachael Ross (2023/04/05)
*******************************************************************************************************************/


/***********************************************
Loading Data */

DATA d;
   INPUT y;
   DATALINES;
7
1
5
3
24
;
RUN;


/***********************************************
Mean by usual / closed-form estimator */
PROC MEANS DATA=d N MEAN;
	VAR y;
RUN;


/***********************************************
M-estimator */
PROC IML;                            /*All steps are completed in PROC IML*/
	*Read data;
	USE d;                           	/*Open input data from above*/
		READ all VAR {y} INTO y;        /*... read the column y in the vector y*/
	CLOSE d;                            /*Close input data from above*/
	n = nrow(y);                        /*Mark number of observations in data (rows of y)*/

	/***********************************************
	Defining estimating equation */
	q = 1;								/*Number of parameters to be estimated*/

	START efunc(theta) global(n, y);	/*Start to define estimating function */
		ef = y - theta;                	/*Estimating function*/
		RETURN(ef);                     /*Return estimating function (n by q matrix)*/
	FINISH efunc;                       /*End definition of estimating function*/

	START eequat(theta);				/*Start to define estimating equation (single argument)*/
		ef = efunc(theta);				/*Call estimating function function*/
		RETURN(ef[+,]);                	/*Return column sums (1 by q vector)*/
	FINISH eequat;                      /*End definition of estimating equation*/

	/***********************************************
	Root-finding */
	mu = {0};                  			* Initial parameter values;
	optn = q || 1;                      * q roots/parameters;
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   mu_hat,                  /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   mu,                      /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                      * Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                mu_hat,                 /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(mu_hat);		      * Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals;   * Outer product of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 			  * Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	b = mu_hat`;                        /*Prepare parameter point estimates for export*/
	se = sqrt(vecdiag(sandwich));       /*Extract corresponding SE for each parameter*/

	* TITLE1 "M-estimator for Mean";    /*Set title for Results Viewer*/
	* PRINT b bread meat sandwich se;   /*Print information to the Results Viewer*/

	CREATE out VAR {b sandwich se};     /*Create an output data set called `out`*/
		APPEND;                         /*... that includes the parameter estimates, variance, and SE*/
	CLOSE out;                          /*Close the output*/
	QUIT;                                   
RUN;

/*Creating confidence intervals*/
DATA out;
	SET out;
	lcl = b - 1.96*se;
	ucl = b + 1.96*se;
RUN;

/*Printing final results to Results Viewer*/
PROC PRINT DATA=out;
RUN;
