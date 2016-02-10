/* Some functions needed for the scenario-generation
   All the functions are independent on the SG-algorithm!!!
   */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"       // my matrix class
#include "misc_macros.h"
#include "sg_functions.h"


// check consistency of the target matrices
/* If it finds an error, it exits the program */
void CheckInput(TMatrix const * const p_TarMoms,
                TMatrix const * const p_TgCorrs,
                TVector const * const p_Probs,
                TMatrix const * const p_OutMat, int const TestLevel)
{
	int s;
	int nVar = p_OutMat->nrow;
	int nScen = p_OutMat->ncol;
	int Chol_TestLevel = TestLevel;
	double sum_prob = 0;

	TMatrix Chol_TgCorrs=Mat_0; // Cholesky decomp. of the target correl. matrix

	// Check dimensions
	if ((p_TarMoms->nrow!=4) | (p_TarMoms->ncol!=nVar) | (p_TgCorrs->nrow!=nVar)
	    | (p_TgCorrs->ncol!=nVar) | (p_Probs->size!=nScen)) {
		printf("\n\tHKW: ERROR - Wrong dimension of input matrices!\n\n");
		exit(1);
	}

	// Check if correlation matrix positive semidefinite
	Mat_Init(&Chol_TgCorrs, nVar,nVar);
	if (Mat_Cholesky(p_TgCorrs, &Chol_TgCorrs, Chol_TestLevel) > 0) {
		printf("\nProblem with the target correlation matrix");
		printf(" - Cholesky failed!\n\n");
		exit(1);
	}
	if (TestLevel > 3) {
		Mat_Display(p_TgCorrs, "R_tg");
		Mat_Display(&Chol_TgCorrs, "L_tg");
	}
	Mat_Kill(&Chol_TgCorrs);

	// Check if probabilities sum-up to one
	for (s=0; s<nScen; s++)
		sum_prob += p_Probs->val[s];
	if (fabs(sum_prob-1.0) > EPS) {
		printf("\nERROR - Input probabilities do not sum-up to one!\n\n");
		exit(1);
	}
}


/* Computes moments and correlations of a given outcomes X,
	 given the probability vector P.
	 Results are stored in arrays Moms and Corrs.
	 NO CONTROL OF DIMENSIONS!!!
	 */
void ComputeProperties(int const nVar, int const nScen, double ** const X,
                       double const * const P,
                       double ** const Moms, double ** const Corrs)
{
	int i,j,s;

	// Compute Moments
	for (i=0; i<4; i++)
		for (j=0; j<nVar; j++) {
			Moms[i][j] = 0;
			for (s=0; s<nScen; s++) {
				if (i==0) // mean
					Moms[i][j] += P[s] * X[j][s];
				else
					Moms[i][j] += P[s] * pow(X[j][s]-Moms[0][j],i+1);
			}
			if (i==1) // variance -> stdev
				Moms[i][j] = sqrt(Moms[i][j]);
			if (i>=2) // skewness og kurtosis
				Moms[i][j] /= pow(Moms[1][j],i+1);
		}

	// Compute Correlations
	for (i=0; i<nVar; i++) {
		for (j=0; j<i; j++) {
			Corrs[i][j] = 0;
			for (s=0; s<nScen; s++) {
				Corrs[i][j] += P[s]
					* (X[i][s] - Moms[0][i])*(X[j][s] - Moms[0][j]);
			}
			Corrs[i][j] /= (Moms[1][i]*Moms[1][j]);
			Corrs[j][i] = Corrs[i][j]; // second half of the matrix
		}
		Corrs[i][i] = 1.0; // the diagonal
	}
}


/* Compute the mean square error of moments
	 NO CONTROL OF DIMENSIONS!!!
	 */
double ErrorMoments(int const nVar, double ** const Mom,
                    double ** const TgMom)
{
	int i,j;
	double error = 0;

	for (i=0; i<4; i++)
		for (j=0; j<nVar; j++)
			error += pow(Mom[i][j]-TgMom[i][j],2);
	error = sqrt(error / (4*nVar));  // sqrt(mean_square_error)

	return(error);
}


/* Compute the mean square error of correlations
	 NO CONTROL OF DIMENSIONS!!!
	 */
double ErrorCorrs(int const nVar, double ** const Corr, double ** const TgCorr)
{
	int i,j;
	double error = 0;

	if (nVar>1) { // nVar=1 -> no correlations
		for (i=0; i<nVar; i++)
			for (j=0; j<i; j++)
				error += pow(Corr[i][j]-TgCorr[i][j],2);
		error = sqrt(error / (nVar*(nVar-1)/2));  // sqrt(mean_square_error)
	}

	return(error);
}


/* Transform the original target moments to the right format
	Create the new moments (used inside the algorithm)
	*/
int CreateTargetMoments(int const nVar, int const nSc, unsigned short Format,
                        double ** const TarMom, double ** const TgMom)
{
	// 1. TRANSFORM THE TARGET MOMENTS (TarMom) TO THE ASSUMED FORMAT
	//     : MEAN, STDEV and the higher moments SCALED by StDev
	const unsigned short PopulEstim = 1;
	const unsigned short VarForStDev = 2;
	const unsigned short KurtoseMinus3 = 4;
	const unsigned short NotScaled = 8;
	const unsigned short NotCentralMoms = 16;
	int i,j;

	// MOMENTS are NOT CENTRAL, i.e. TarMom[i-1]=E{X^i} -> CENTRALISE
	if (Format & NotCentralMoms) {
		for (i=0; i<nVar; i++) {
			TarMom[1][i] += - pow(TarMom[0][i],2); // Variance
			TarMom[2][i] += - 3*TarMom[0][i]*TarMom[1][i] - pow(TarMom[0][i],3); // Skewness
			TarMom[3][i] += - 4*TarMom[0][i]*TarMom[2][i] - 6*pow(TarMom[0][i],2)*TarMom[1][i] - pow(TarMom[0][i],4); // Kurtose
			Format = NotScaled + VarForStDev; // Now I have TarMom[i-1]==E{(X-E[X])^i}, i>1
		}
	}
	// CENTRAL MOMENTS
	// Control variance
	for (i=0; i<nVar; i++) {
		if (TarMom[1][i]<=0) {
			printf("\n\nERROR : Variance of variable %d is negative!!!\n\n", i+1);
			exit(1);
		}
	}
	if (Format & VarForStDev) {
		// Variance insted of Standard Deviation
		for (i=0; i<nVar; i++)
			TarMom[1][i] = sqrt(TarMom[1][i]);
	}
	if (Format & NotScaled) {
		// Higher moments are not scaled by Standard Deviation
		for (i=0; i<nVar; i++)
			for (j=2; j<4; j++)
				TarMom[j][i] /= pow(TarMom[1][i],j+1);
	}
	if ((Format & KurtoseMinus3) && !(Format & PopulEstim)) {
		// The 4. moment is Kurtose - 3
		for (i=0; i<nVar; i++)
			TarMom[3][i] += 3;
	}

	// 2. IF REQUIRED, CHANGE THE MOMENTS SO THAT THE USER GETS THE CORRECT
	//    VALUES USING THE STANDARD ESTIMATORS ON SCENARIOS (IN A SPREADSEET)
	//    For formulas, see OpenOffice Calc documentation for skew and kurt:
	//    http://wiki.services.openoffice.org/wiki/Documentation/How_Tos/Calc:_SKEW_function
	//    http://wiki.services.openoffice.org/wiki/Documentation/How_Tos/Calc:_KURT_function
	double N = (double) nSc; // need it as double for division!
	if (Format & PopulEstim) {
		for (i=0; i<nVar; i++) {
			// standard deviation
			TarMom[1][i] /= sqrt(N / (N - 1)); // correct the multiplier
			// skewness:
			TarMom[2][i] *= (N - 1) * (N - 2) / N  / N; // correct sum multiplier
			TarMom[2][i] /= pow((N - 1) / N, 1.5);      // correct std. dev.
			// kurtosis
			TarMom[3][i] += 3 * (N - 1) * (N - 1) / ((N - 2) * (N - 3));
			TarMom[3][i] *= (N - 1) * (N - 2) * (N - 3) / (N * (N + 1))  /N;
			TarMom[3][i] /= pow((N - 1) / N, 2);
		}
	}

	// 3. CHECK CONSISTENCY
	// Control kurtose: kurt > 1 + skew^2
	for (i=0; i<nVar; i++) {
		if (TarMom[3][i] <= 1 + pow(TarMom[2][i],2)) {
			printf("\n\nERROR : Kurtosis of variable %d is too small!!!\n\n", i+1);
			printf("\tKurtosis must be: kurt > var^2 + skew^2/var\n\n");
			printf("\t!!! Control the format of input data !!!\n\n");
			exit(1);
		}
	}

	// 4. CREATE THE MOMENTS FOR INSIDE THE ALGORITHM - TgMom
	for (i=0; i<nVar; i++) {
		TgMom[0][i] = 0;
		TgMom[1][i] = 1;
		for (j=2; j<4; j++)
			TgMom[j][i] = TarMom[j][i];
	}

	return(0);
}


/// Probability density function of N(0,1)
double density_N01(double const x)
{
	return 1/sqrt(2*Pi) * exp(-x*x/2);
}

/// Inverse of the N(0,1) distribution function (polynomial approximation)
double inv_distr_N01(double p)
{
	const double tlow = 1e-307;
	const double thgh = 1 - 1e-16;
	double x,t,n,d;

	if (p<tlow) p = tlow;		//p = tlow*(p<tlow) + p*(p>=tlow);
	if (p>thgh) p = thgh;		//p = thgh*(p>thgh) + p*(p<=thgh);
	t = sqrt(-2*log(fabs((p>0.5)-p)));
  n = 2.515517 + t*(0.802853 + t*0.010328);			//n=2.52+0.80*t+0.01*t^2
  d = 1 + t*(1.432788 + t*(0.189269 + t*0.001308));	//d=1+1.4*t+0.2*t^2+0.001*t^3
  x = t - (n/d);
  if (p<=0.5) x*=-1;		//x=x*(p>0.5) - x*(p<=0.5);
  return(x);
}


/* Fill a given array with numbers from N(0,1)
   Can use either discretizing of random number generation
	*/
void ArrayOfN01(double x[], int const nVal, int const UseDiscretizing)
{
	int i, index;

	if (!UseDiscretizing) {
		for (i=0; i<nVal; i++)
			x[i] = BM_normal();
	}
	else {
		double perc_low, perc_high=1;
		double outc_low, outc_high=1e10;
		double temp;
		// first disretize
		for (i=0; i<nVal; i++) {
			if (i==0) {
				perc_low = 0;
				outc_low = -1e10;
			}
			else {
				// copy values from previous interval
				perc_low = perc_high;
				outc_low = outc_high;
			}
			if (i==nVal-1) {
				perc_high = 1;
				outc_high = 1e10;
			}
			else {
				perc_high = (double) (i+1)/nVal;
				outc_high = inv_distr_N01(perc_high);
			}
			// x is the conditional mean of the interval
			x[i] = (density_N01(outc_high)-density_N01(outc_low)) / (perc_high-perc_low);
		}

		// now randomize the order
		for (i=0; i<nVal; i++) {
			index = i + irand(nVal - i);
			temp = x[i];
			x[i] = x[index];
			x[index] = temp;
		}
	}
}


// Transforms vector to zero mean and a variance of one
void NormalizeVector(int const size, double vector[], double prob_vec[]) {
	int i;
	double mean = 0;
	double stdev = 0;

	// compute moments
	if (prob_vec == NULL) {
		// no probabilities -> equiprobable
		for (i=0; i<size; i++) {
			mean += vector[i];
			stdev += pow(vector[i],2);
		}
		mean /= (double) size;
		stdev /= (double) size;
	} else {
		for (i=0; i<size; i++) {
			mean += prob_vec[i] * vector[i];
			stdev += prob_vec[i] * pow(vector[i],2);
		}
	}
	stdev = sqrt(stdev-pow(mean,2));

	// normalize
	for (i=0; i<size; i++)
		vector[i] = (vector[i] - mean) / stdev;
}
