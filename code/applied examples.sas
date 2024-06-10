/*******************************************************************************************************************
ABC's of M-estimation
	Code for applied examples (Section 2)

Rachael Ross (2023/03/28), Paul Zivich (2024/06/10)
*******************************************************************************************************************/

/***********************************************
Example 1: Logistic Regression */
/***********************************************

/***********************************************
Loading Data */

data dat_cnts;
input anemia bp ptb n;
datalines;
0 0 0 496
0 0 1 74
0 1 0 113
0 1 1 25
1 0 0 85
1 0 1 15
1 1 0 15
1 1 1 3
;
run;

data dat;
set dat_cnts;
do i=1 to n;
	output;
end;
drop n;
run;


/***********************************************
Rregression by MLE */
proc logistic data=dat;
model ptb(ref='0')= anemia bp;
ods output parameterestimates=ests_mle;
run;

proc contents data=ests_mle;
run;

data ests_mle;
set ests_mle(keep=variable estimate stderr);
lcl = estimate - 1.96*stderr;
ucl = estimate + 1.96*stderr;
run;

proc print data=ests_mle noobs;
run;

/***********************************************
M-estimator */
PROC IML;                            /*All steps are completed in PROC IML*/
	*Read data;
	use dat;								/*Open input data from above*/
		read all var {ptb} into ptb;		/*Read in each column as its own vector*/
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        /*Save number of observations in data */

	/***********************************************
	Defining estimating equation */

	q = 3;								/*Save number parameters to be estimated*/

	START efunc(beta) global(ptb, anemia, bp);							/*Start to define estimating function */ 
		p = 1 / (1 + exp(-(beta[1] + beta[2]*anemia + beta[3]*bp)));	
		ef_1 = ptb - p;
		ef_2 = (ptb - p)#anemia;
		ef_3 = (ptb - p)#bp;
		ef_mat = ef_1||ef_2||ef_3;
		RETURN(ef_mat);                         						/*Return n by q matrix for estimating functions*/
	FINISH efunc;                       								/*End definition of estimating equation*/

	START eequat(beta);					 								/*Start to define estimating equation (single argument)*/ 
		ef = efunc(beta);
		RETURN(ef[+,]);                  								/*Return column sums, 1 by q vector)*/
	FINISH eequat;                       								/*End definition of estimating equation*/

	/***********************************************
	Root-finding */
	initial = {-2,0,0};                 * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   beta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	print beta_hat;

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                beta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	print bread;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(beta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	print filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	print sandwich;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","anemia","bp"};  
	est = beta_hat`;                    
	se = sqrt(vecdiag(sandwich));       /*Extract corresponding SE for each parameter*/
	lcl = est - 1.96*se; 				/*Calculated lcl*/
	ucl = est + 1.96*se;				/*Calculate ucl*/

	PRINT variable est se lcl ucl;     	/*Print information to the Results Viewer*/

	CREATE ests_mest VAR {variable est se lcl ucl};   /*Create an output data set called `out`*/
		APPEND;                         		  	  /*... that includes the parameter estimates, variance, and SE*/
	CLOSE ests_mest;                          		  /*Close the output*/
	QUIT;                                   
RUN;


/***********************************************
Results for logistic regression */

title1 "Logistic regression: MLE";
proc print data=ests_mle noobs; run;

title1 "Logistic regression: M-estimation";
proc print data=ests_mest noobs; run;

title1;


/***********************************************
Example 2a: Standardization by g-computation */
/***********************************************

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use dat;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        

	/***********************************************
	Defining estimating equation */

	q = 7;								

	START efunc(theta) global(ptb, anemia, bp, n);							
		p = 1 / (1 + exp(-(theta[1] + theta[2]*anemia + theta[3]*bp)));	
		ef_1 = ptb - p;
		ef_2 = (ptb - p)#anemia;
		ef_3 = (ptb - p)#bp;

		ef_r1 = 1 / (1 + exp(-(theta[1] + theta[2]*1 + theta[3]*bp))) - theta[4];
		ef_r0 = 1 / (1 + exp(-(theta[1] + theta[2]*0 + theta[3]*bp))) - theta[5];

		ef_rd = j(n,1,(theta[4] - theta[5]) - theta[6]);
		ef_lnrr = j(n,1,log(theta[4]/theta[5]) - theta[7]);

		ef_mat = ef_1||ef_2||ef_3||ef_r1||ef_r0||ef_rd||ef_lnrr;
		RETURN(ef_mat);                         						
	FINISH efunc;                       								

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {-2,0,0,.1,.1,0,0};       * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","anemia","bp","risk1","risk0","rd","lnrr"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_2a VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_2a;                          		  
	QUIT;                                   
RUN;


/***********************************************
Example 2b: Standardization by IPW */
/***********************************************

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use dat;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        

	/***********************************************
	Defining estimating equation */

	q = 6;								

	START efunc(theta) global(ptb, anemia, bp, n);							
		pscore = 1 / (1 + exp(-(theta[1] + theta[2]*bp)));	
		ef_1 = anemia - pscore;
		ef_2 = (anemia - pscore)#bp;

		wt = anemia/pscore + (1-anemia)/(1-pscore);
		ef_r1 = anemia#wt#ptb - theta[3];
		ef_r0 = (1-anemia)#wt#ptb - theta[4];

		ef_rd = j(n,1,(theta[3] - theta[4]) - theta[5]);
		ef_lnrr = j(n,1,log(theta[3]/theta[4]) - theta[6]);

		ef_mat = ef_1||ef_2||ef_r1||ef_r0||ef_rd||ef_lnrr;
		RETURN(ef_mat);                         						
	FINISH efunc;                       								

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {-2,0,.1,.1,0,0};       * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","bp","risk1","risk0","rd","lnrr"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_2b VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_2b;                          		  
	QUIT;                                   
RUN;

/***********************************************
Example 3: Transport to External Target */
/***********************************************

/***********************************************
Loading Target Data */

data d_target_cnt;
input zapps anemia bp n;
datalines;
1 0 0 300 
1 0 1 300
1 1 0 300
1 1 1 100
;
run;

data d_target;
set d_target_cnt;
do i=1 to n;
	output;
end;
drop n;
run;

* Adding ZAPPS indicator;
data dat;
set dat;
zapps = 1;
run;

data d_target;
set d_target;
ptb = -999; *setting outcomes to generic missing value;
zapps = 0;
run;

* Stacking data sets;
data df;
set dat d_target;
run;

/***********************************************
Example 3a: Transport by g-computation */
/***********************************************

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use df;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
		read all var {zapps} into s;		
	close df;

	n = nrow(s);                        

	/***********************************************
	Defining estimating equation */

	q = 4;

	START efunc(theta) global(s, ptb, anemia, bp);
		*Outcome nuisance model;
		p = 1 / (1 + exp(-(theta[1] + theta[2]*anemia + theta[3]*bp)));	
		ef_1 = s#(ptb - p);
		ef_2 = s#((ptb - p)#anemia);
		ef_3 = s#((ptb - p)#bp);
		*Risk estimating function;
		ef_r = (1-s)#(p - theta[4]);
		* Stacking estimating functions;
		ef_mat = ef_1||ef_2||ef_3||ef_r;
		RETURN(ef_mat);                         						
	FINISH efunc;       

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {.7,1,1,.7};       * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Beta_0","Beta_1","Beta_2","Risk"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_3 VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_3;                          		  
	QUIT;                                   
RUN;

/***********************************************
Example 3b: Standardization by IOSW */

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use df;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
		read all var {zapps} into s;		
	close df;

	n = nrow(s);                        

	/***********************************************
	Defining estimating equation */

	q = 4;

	START efunc(theta) global(s, ptb, anemia, bp);
		*Outcome nuisance model;
		p = 1 / (1 + exp(-(theta[1] + theta[2]*anemia + theta[3]*bp)));	
		weight = s#(1-p)/p;
		ef_1 = s - p;
		ef_2 = (s - p)#anemia;
		ef_3 = (s - p)#bp;
		*Risk estimating function;
		ef_r = (ptb - theta[4]) # weight;
		* Stacking estimating functions;
		ef_mat = ef_1||ef_2||ef_3||ef_r;
		RETURN(ef_mat);                         						
	FINISH efunc;       

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {.7,1,1,.7};       * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Beta_0","Beta_1","Beta_2","Risk"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_3 VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_3;                          		  
	QUIT;                                   
RUN;


/*END*/
