/*******************************************************************************************************************
ABC's of M-estimation
	Code for M-estimator of the mean (Section 1)

Paul Zivich (2023/03/22)
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
	q = 1;
	START efunc(theta) global(n, y);	/*Start to define estimating equation (single argument)*/ 
		ee = j(1, 1, 0);                	/*Create storage for estimating equation evaluation*/
		DO i = 1 TO n;                      /*Loop through all observations*/
			ee = ee + (y[i] - theta);       /*... evaluate for observation i, then add to stored estimating equation evaluation*/
		END;                                /*End the loop*/
		RETURN(ee);                         /*Return the estimating equation (\sum_i^n \psi(O_i, \theta)*/
	FINISH efunc;                       /*End definition of estimating equation*/

	/***********************************************
	Root-finding */
	mu = j(1, q, 0.0);                  * Initial parameter values;
	optn = q || 1;                      * q roots/parameters;
	tc = j(1, 12, .);                   * Missing values set to defaults for root-finder;
	tc[6] = 1e-9;                       * Default is 1e-5, but making more stringent here;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   mu_hat,                  /*... name of output parameters that give the root*/
			   "efunc",                 /*... function to find the roots of*/
			   mu,                      /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = j(1, 3, .);                   * par is a length 3 vector of details;
	par[1] = q;                         * tell FD we have q parameters;
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                bread,                  /*... name of output matrix that gives the bread*/
                hess,                   /*... */
                "efunc",                /*... function to approximate the derivative of*/
                mu_hat,                 /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (bread) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	meat = j(q, q, 0);                  * Defining matrix of zeroes to store the calculated filling;
	ef = j(q, 1, 0);                    * Initialize storage for evaluated estimating function;
	DO i = 1 TO n;                      /*Start a loop over the n observations*/
		ef = y[i] - mu_hat;             /*... evaluate the estimating function at parameter estimate*/
		meat = meat + ef * ef`;         /*... take dot product and add back to filling matrix storage*/
	END;                                /*End the loop*/
	meat = meat / n;                    * Average the filling matrix by dividing by n;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * meat * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	b = mu_hat`;                        /*Prepare parameter point estimates for export*/
	se = sqrt(vecdiag(sandwich));       /*Extract corresponding SE for each parameter*/

	* TITLE1 "M-estimator for Mean";      /*Set title for Results Viewer*/
	* PRINT b bread meat sandwich se;     /*Print information to the Results Viewer*/

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
