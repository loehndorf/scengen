/*********************************************************************/
/***  Høyland-Kaut-Wallace algorithm for scenario generation       ***/
/***                                                               ***/
/***  author: Michal Kaut                                          ***/
/*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "HKW_sg.h"
#include "HKW_cubic.h"    // cubic transformation
#include "misc_macros.h"  // misc. C macros used here (like frand)
#include "sg_functions.h" // general functions needed for SG


#define StartWithCubic 1 /// Start with the cubic distribution
#define MaxStartTrial 5  /// Max number of starting points per marginal
#define MaxCubIter 2     /// Max number of cubic transf. per HKW iteration


/* MAIN GENERATION FUNCTION
	 Generates one subtree using the HKW algorithm (heuristics)
	 Parameter FormatOfMoms describes the format of target moments.
		It is a sum of following bits:
			1 -> 2nd moment is Var instead of StDev
			2 -> 4th moment is Kurtose - 3
			4 -> Higher moments are not scaled by StDev
			8 -> TarMom[i-1] = E{X^i} ... lower bits are ignored in this case
		=> The reasonable values are 0,1,2,3,6,7,8
	*/
int HKW_ScenGen(int const FormatOfMoms, TMatrix const * const p_TarMoms,
                TMatrix const * const p_TgCorrs, TVector const * const p_Probs,
                TMatrix * const p_OutMat,
                double const MaxErrMom, double const MaxErrCorr,
                int const TestLevel, int const MaxTrial, int const MaxIter,
                int const UseStartValues,
                double * p_errMom, double * p_errCorr,
                int * p_nmbTrial, int * p_nmbIter)
{
	int i,j,s;
	int err_code;

	int nVar = p_OutMat->nrow;
	int nScen = p_OutMat->ncol;

	int iter;
	int trial = 0, start_trial = 0;
	double errMom, errCorr, bestError;
	double tmp_double;

	TMatrix OutMoms=Mat_0, OutCorrs=Mat_0, Chol_OutCorrs=Mat_0,
	        InvChol_OutCorrs=Mat_0, TrsfMat=Mat_0, TmpOutMat=Mat_0;
	TMatrix Chol_TgCorrs=Mat_0;
	TMatrix TgMoms=Mat_0; // The transformed (normalised) moments

	int Chol_TestLevel = TestLevel;

	// Variables for the cubic transformation
	// input data - defined in HKW_cubic.h
	//extern double InMom[13]; // Input moments + In[0]:=1
	//extern double TgMom[4];  // Target moments
	double CubParam[4];        // This is what we want to find
	int cub_iter;
	double cub_error;
	InMom[0] = 1;


	// Check consistency of the input
	CheckInput(p_TarMoms, p_TgCorrs, p_Probs, p_OutMat, TestLevel);

	// Compute the Cholesky decomposition
	// (done in CheckInput(), but not returned)
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

	// Transform the target moments
	Mat_InitAsBigAs(&TgMoms, p_TarMoms);
	CreateTargetMoments(nVar, nScen, (unsigned short) FormatOfMoms,
	                    p_TarMoms->val, TgMoms.val);
	if (TestLevel > 3) {
		Mat_Display(p_TarMoms, "TARMOM");
		Mat_Display(&TgMoms, "MOM");
	}


	// Matrices for the algorithm
	Mat_Init(&OutMoms, 4,nVar);  // moments of outcomes
	Mat_Init(&OutCorrs, nVar,nVar);  // correlations of outcomes
	Mat_Init(&Chol_OutCorrs, nVar,nVar);  // Cholesky decomp of OutCorrs
	Mat_Init(&InvChol_OutCorrs, nVar,nVar);  // inverse of Chol_OutCorrs
	Mat_Init(&TrsfMat, nVar,nVar);  // matrix of the linear transformation
	Mat_InitAsBigAs(&TmpOutMat, p_OutMat);  // temporary outcome matrix


	// TRIALS
	do {

		trial++;
		if (TestLevel > 1)
			printf("\ntrial %d\n", trial);

		if (UseStartValues) {
			// Using starting values
			// !!! HAVE TO TRANSFORM THEM TO MEAN=0, STDEV=1 !!!
			for (i=0; i<nVar; i++)
				for (s=0; s<nScen; s++)
					p_OutMat->val[i][s] = (p_OutMat->val[i][s] - p_TarMoms->val[0][i])
					                      / p_TarMoms->val[1][i];
		}
		else {
			// GENERATING STARTING VALUES
			for (i=0; i<nVar; i++) {

				if (!StartWithCubic) {
					// NOT starting with cubic -> first is the matrix transformation!
					// Start with N(0,1) values, discretizing (FIXED values)
					ArrayOfN01(p_OutMat->val[i], nScen, 1);
				}
				else {
					// Start with cubic
					//  - do it by marginals, test if we can get the right moments directly!!!

					// target moments
					TgMom[0] = TgMoms.val[0][i];
					TgMom[1] = pow(TgMoms.val[1][i],2);
					TgMom[2] = TgMoms.val[2][i] * pow(TgMoms.val[1][i],3);
					TgMom[3] = TgMoms.val[3][i] * pow(TgMoms.val[1][i],4);

					if (TestLevel > 6) {
						printf("CUBIC - Target:");
						for (j=0; j<4; j++)
							printf("  %.4f", TgMom[j]);
						printf("\n");
					}

					start_trial = 0;
					bestError = 1e10;
					do { // trials - different starting distributions
						start_trial++;

						switch (start_trial) {
							case 1:
								// Start with N(0,1) values, discretizing (FIXED values)
								ArrayOfN01(TmpOutMat.val[i], nScen, 1);
								break;
							case MaxStartTrial:
								// Last trial - try uniform instead of N(0,1)!
								for (s=0; s<nScen; s++)
									TmpOutMat.val[i][s] = frand;
								break;
							default:
								// Start with N(0,1) values, RANDOM
								ArrayOfN01(TmpOutMat.val[i], nScen, 0);
								break;
						}

						// CUBIC TRANSFORMATION
						// Repeat the cubic transformation several times if needed
						cub_iter = 0;
						do {
							cub_iter++;

							// compute the 12 moments
							for (j=1; j<=12; j++) {
								InMom[j] = 0;
								for (s=0; s<nScen; s++)
									InMom[j] += p_Probs->val[s]
									            * pow(TmpOutMat.val[i][s],j);
							}
							if (TestLevel > 6) {
								printf("CUBIC - Start:");
								for (j=1; j<=12; j++)
									printf("  %.4f", InMom[j]);
								printf("\n");
							}

							// find the parameters
							cub_error = cubic_solve(CubParam);

							if (TestLevel > 2) {
								printf("Error of the cubic transformation of r.v.");
								printf(" %2d = %g\n", i+1, cub_error);
							}

							if ((cub_error > EPS) && (cub_iter < MaxCubIter)) {
								// transform for new cub_iter; can not use OutMat!
								// using a+bX+cX^2+dX^3 = ((d*X+c)*X+b)*X+a
								for (s=0; s<nScen; s++) {
									tmp_double = CubParam[3]; // d
									for (j=2; j>=0; j--) {
										tmp_double *= TmpOutMat.val[i][s]; // *X
										tmp_double += CubParam[j]; // + c, b, a
									}
									TmpOutMat.val[i][s] = tmp_double;
								} // next s
							} // end if
						} while ((cub_error>EPS) && (cub_iter<MaxCubIter));

						// out of the cub_iter loop
						// cub_error < EPS  OR  cub_iter = MaxCubIter
						if (cub_error < bestError) {
							bestError = cub_error;
							// new best solution -> copy to OutMat
							// do the transformation - directly to OutMat
							// using a+bX+cX^2+dX^3 = ((d*X+c)*X+b)*X+a
							for (s=0; s<nScen; s++) {
								p_OutMat->val[i][s] = CubParam[3]; // d
								for (j=2; j>=0; j--) {
									p_OutMat->val[i][s] *= TmpOutMat.val[i][s]; // *X
									p_OutMat->val[i][s] += CubParam[j]; // + c, b, a
								}
							}
						} // end if
					} while ((cub_error>EPS) && (start_trial<MaxStartTrial));

					if (TestLevel > 0) {
						if ((cub_error<EPS) && (start_trial==MaxStartTrial)) {
							printf("Start values for r.v. %d", i+1);
							printf(" were generated with uniform distribution!\n");
						}
						if ((cub_error>EPS) && (start_trial==MaxStartTrial)) {
							printf("Warning: Did NOT find a good starting point for");
							printf(" r.v. %2d (error = %g)!\n", i+1, bestError);
						}
					}
				} // end: if (StartWithCubic)

			} // next i (starting values for next marginals)
		} // end: if (UseStartValues)


		// Compute starting properties and errors
		ComputeProperties(nVar, nScen, p_OutMat->val, p_Probs->val,
		                  OutMoms.val, OutCorrs.val);
		errMom = ErrorMoments(nVar, OutMoms.val, TgMoms.val);
		errCorr = ErrorCorrs(nVar, OutCorrs.val, p_TgCorrs->val);

		if (TestLevel > 6) {
			if (TestLevel > 10) Mat_DisplayTransp(p_OutMat, "outcomes");
			Mat_DisplayTransp(&OutMoms, "moments");
			Mat_DisplayTransp(&TgMoms, "tg_moms");
			Mat_DisplayTransp(&OutCorrs, "corrs");
			Mat_DisplayTransp(p_TgCorrs, "tg_corrs");
		}
		if (TestLevel > 2) {
			printf("Starting properties:\n");
			printf("  Error in moments is %.4f:\n", errMom);
			if (TestLevel > 3) Mat_Display(&OutMoms, "Moments");
			printf("  Error in correlations is %.4f:\n", errCorr);
			if (TestLevel > 3) Mat_Display(&OutCorrs, "Corrs");
		}


		// MAIN ITERATIONS OF HKW-alg
		iter=0;
		while (((errMom>MaxErrMom) || (errCorr>MaxErrCorr)) && (iter<MaxIter)) {
			iter++;
			if (TestLevel > 2)
				printf("\niter %d\n", iter);

			// Cholesky decomp. of correlation matrix of outcomes
			err_code = Mat_Cholesky(&OutCorrs, &Chol_OutCorrs, Chol_TestLevel);
			if (err_code > 0) {
				printf("\nError in Cholesky decomposition of correlation matrix");
				printf(" of outcomes!\n");
				printf(" - exit code was %d\n\n", err_code);
				goto next_trial;
			}
			Mat_LowTriangInverse(&Chol_OutCorrs, &InvChol_OutCorrs);  // inverse
			// this creates the matrix we multiply with
			Mat_Mult_LeftTriang(&Chol_TgCorrs, 'L', &InvChol_OutCorrs, &TrsfMat);

			if (TestLevel > 3) {
				Mat_Display(&OutCorrs, "R");
				Mat_Display(&Chol_OutCorrs, "L");
				Mat_Display(&InvChol_OutCorrs, "invL");
				Mat_Display(&Chol_TgCorrs, "L_tg");
				Mat_Display(&TrsfMat, "L_iter");
			}

			// THE LINEAR (MATRIX) TRANSFORMATION
			if (Mat_Mult_LeftTriang(&TrsfMat, 'L', p_OutMat, &TmpOutMat) > 0) {
				printf("ERROR during the linear transformation!\n\n");
				goto next_trial;
			}

			// Compute properties and errors
			ComputeProperties(nVar, nScen, TmpOutMat.val, p_Probs->val,
			                  OutMoms.val, OutCorrs.val);
			errMom = ErrorMoments(nVar, OutMoms.val, TgMoms.val);
			errCorr = ErrorCorrs(nVar, OutCorrs.val, p_TgCorrs->val);

			if (TestLevel > 2) {
				printf("After the linear (matrix) transformation:\n");
				printf("  Error in moments is %.4f:\n", errMom);
				if (TestLevel > 3) Mat_Display(&OutMoms, "Moments");
				printf("  Error in correlations is %.4f:\n", errCorr);
				if (TestLevel > 3) Mat_Display(&OutCorrs, "Corrs");
				printf("\n");
			}


			// THE CUBIC TRANSFORMATION - separately for every r.v.
			for (i=0; i<nVar; i++) {

				// set the target moments
				TgMom[0] = TgMoms.val[0][i];
				TgMom[1] = pow(TgMoms.val[1][i],2);
				TgMom[2] = TgMoms.val[2][i] * pow(TgMoms.val[1][i],3);
				TgMom[3] = TgMoms.val[3][i] * pow(TgMoms.val[1][i],4);

				if (TestLevel > 6) {
					printf("CUBIC - Target:");
					for (j=0; j<4; j++)
						printf("  %.4f", TgMom[j]);
					printf("\n");
				}

				cub_iter = 0;
				do { // Repeat the cubic transformation several times if needed
					cub_iter++;

					// compute the 12 moments
					for (j=1; j<=12; j++) {
						InMom[j] = 0;
						for (s=0; s<nScen; s++)
							InMom[j] += p_Probs->val[s] * pow(TmpOutMat.val[i][s],j);
					}
					if (TestLevel > 6) {
						printf("CUBIC - Start:");
						for (j=1; j<=12; j++)
							printf("  %.4f", InMom[j]);
						printf("\n");
					}

					// find the parameters
					cub_error = cubic_solve(CubParam);

					if (TestLevel > 2) {
						printf("Error of the cubic transformation of r.v.");
						printf(" %2d = %g\n", i+1, cub_error);
					}

					// do the transformation
					// using a+bX+cX^2+dX^3 = ((d*X+c)*X+b)*X+a
					for (s=0; s<nScen; s++) {
						p_OutMat->val[i][s] = CubParam[3]; // d
						for (j=2; j>=0; j--) {
							p_OutMat->val[i][s] *= TmpOutMat.val[i][s]; // *X
							p_OutMat->val[i][s] += CubParam[j]; // + c, b, a
						}
					}

					// if we are going to repeat cubic, have to copy data!!!
					if ((cub_error>EPS) && (cub_iter<MaxCubIter))
						for (s=0; s<nScen; s++)
							TmpOutMat.val[i][s] = p_OutMat->val[i][s];

				} while ((cub_error>EPS) && (cub_iter<MaxCubIter));
			} // next variable (to cubic transform)
			if (TestLevel > 2) printf("\n");

			// Compute properties and errors
			ComputeProperties(nVar, nScen, p_OutMat->val, p_Probs->val,
			                  OutMoms.val, OutCorrs.val);
			errMom = ErrorMoments(nVar, OutMoms.val, TgMoms.val);
			errCorr = ErrorCorrs(nVar, OutCorrs.val, p_TgCorrs->val);

			if (TestLevel > 2) {
				printf("After the cubic transformation:\n");
				printf("  Error in moments is %.4f:\n", errMom);
				if (TestLevel > 3) Mat_Display(&OutMoms, "Moments");
				printf("  Error in correlations is %.4f:\n", errCorr);
				if (TestLevel > 3) Mat_Display(&OutCorrs, "Corrs");
			}

			if (TestLevel > 1)
				printf("Distance after iteration %2d = %9.6f\n", iter,
				       errMom + errCorr);
		} // next iter

		if (TestLevel > 0) {
			printf("Distance in trial %2d = %9.6f\n", trial, errMom+errCorr);
			if (TestLevel > 2) printf("\n");
		}
		next_trial:;
	} while (((errMom>MaxErrMom) || (errCorr>MaxErrCorr)) && (trial<MaxTrial));


	// De-allocate matrices
	Mat_Kill(&TgMoms);
	Mat_Kill(&Chol_TgCorrs);
	Mat_Kill(&OutMoms);
	Mat_Kill(&OutCorrs);
	Mat_Kill(&Chol_OutCorrs);
	Mat_Kill(&InvChol_OutCorrs);
	Mat_Kill(&TrsfMat);
	Mat_Kill(&TmpOutMat);


	// Re-scale the outcomes to the original moments
	// Assumes that the second target moment is STANDARD DEVIATION
	//     and that the higher target moments are scaled!!!
	for (i=0; i<nVar; i++)
		for (s=0; s<nScen; s++)
			p_OutMat->val[i][s] = p_TarMoms->val[0][i]
			 + p_TarMoms->val[1][i] * p_OutMat->val[i][s];

	// copy results to the output pointers
	if (p_errMom != NULL) *p_errMom = errMom;
	if (p_errCorr != NULL) *p_errCorr = errCorr;
	if (p_nmbTrial != NULL) *p_nmbTrial = trial;
	if (p_nmbIter != NULL) *p_nmbIter = iter;

	return ( (errMom>MaxErrMom) | (errCorr>MaxErrCorr) ? 1 : 0 );
}


int scengen_HKW(double ** const tgMoms, int const FormatOfMoms,
                double ** const tgCorrs, double * const probs,
                int const nVar, int const nScen, double ** outSc,
                double const MaxErrMom, double const MaxErrCorr,
                int const TestLevel, int const MaxTrial, int const MaxIter,
                int const UseStartValues,
                double * p_errMom, double * p_errCorr,
                int * p_nmbTrial, int * p_nmbIter)
{
	// create temp. matrices and vectors, directly using the provided arrays
	// -> no allocation, unless we get NULL pointers
	TMatrix const TarMoms={4, nVar, tgMoms};
	TMatrix const TgCorrs={nVar, nVar, tgCorrs};
	TVector Probs={nScen, probs};
	if (probs == NULL) {
		Vec_Init(&Probs, nScen);
		double p = 1.0 / (double) nScen;
		int s;
		for (s = 0; s < nScen; ++s) {
			Probs.val[s] = p;
		}
	}
	TMatrix OutMat={nVar, nScen, outSc};
	if (outSc == NULL) {
		Mat_Init(&OutMat, nVar, nScen);
		outSc = OutMat.val;
	}

	int retVal = HKW_ScenGen(FormatOfMoms, &TarMoms, &TgCorrs, &Probs, &OutMat,
	                         MaxErrMom, MaxErrCorr, TestLevel, MaxTrial,
	                         MaxIter, UseStartValues,
	                         p_errMom, p_errCorr, p_nmbTrial, p_nmbIter);

	/* De-allocate the temp. vector of probabilities
	   Note that we cannot kill OutMat, since this would de-allocate the data
	   in outSc! It follows that it is up to the user to free this memory! */
	if (!probs)
		Vec_Kill(&Probs);

	return retVal;
}
