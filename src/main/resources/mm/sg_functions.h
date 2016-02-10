// header file for sg_functions.c

#ifndef SG_FUNCTIONS_H
#define SG_FUNCTIONS_H

#include "matrix.h"


/* Some functions needed for the scenario-generation
   All the functions are independent on the SG-algorithm!!!
   */

/// check consistency of the target matrices
/** \warning If it finds an error, it exits the program **/
void CheckInput(TMatrix const * const p_TarMoms,
                TMatrix const * const p_TgCorrs,
                TVector const * const p_Probs,
                TMatrix const * const p_OutMat, int const TestLevel);

/// computes moments and correlations of a given outcomes X,
/**
	\param[in] nVar number of random variables
	\param[in] nScen number of scenarios
	\param[in] X array (2D) of scenario values
	\param[in] P array (1D) of scenario probabilities
	\param[out] Moms output array (2D) of moments
	\param[out] Corrs output array (2D) of correlations
	\warning NO CONTROL OF DIMENSIONS!!!
**/
void ComputeProperties(int const nVar, int const nScen, double ** const X,
                       double const * const P,
                       double ** const Moms, double ** const Corrs);


/// compute the mean square error of moments
/**
	\param[in] nVar number of random variables
	\param[in] Mom array (2D) of actual moments
	\param[in] TgMom array (2D) of target moments
	\return the error, i.e. the mean-square difference
	\warning NO CONTROL OF DIMENSIONS!!!
**/
double ErrorMoments(int const nVar, double ** const Mom,
                    double ** const TgMom);

/// compute the mean square error of moments
/**
	\param[in] nVar number of random variables
	\param[in] Corr array (2D) of actual correlations
	\param[in] TgCorr array (2D) of target correlations
	\return the error, i.e. the mean-square difference
	\warning NO CONTROL OF DIMENSIONS!!!
**/
double ErrorCorrs(int const nVar, double ** const Corr, double ** const TgCorr);


/// transform the original target moments to new ones with the right format
/**
	\param[in] nVar number of random variables
	\param[in] nSc number of scenarios - used only for population estimators
	\param[in] Format format of the target moments - see \c HKW_ScenGen
	\param[in] TarMom input array (2D) of the target moments in format \a Format
	\param[out] TgMom output array (2D) in the "right" format
	\return zero on success
	\warning On problems exits the code!
	\warning Changes values of \a Format and \a TarMom!
**/
int CreateTargetMoments(int const nVar, int const nSc, unsigned short Format,
                        double ** const TarMom, double ** const TgMom);


/// fill a given array with numbers from N(0,1)
/**
	\param x array to fill
	\param[in] nVal number of values
	\param[in] UseDiscretizing if 1, use fixed discretization (in random order),
	           otherwise generate random numbers
	\warning assumes \a x has been allocated - no control!
**/
void ArrayOfN01(double x[], int const nVal, int const UseDiscretizing);


/// transform vector to zero mean and a variance of one (in-place)
/**
	\param[in] size size/length of the vector
	\param vector the array of samples
	\param prob_vec array of probabilities; NULL means equiprobable
	\warning assumes both vectors have been allocated - no control!
**/
void NormalizeVector(int const size, double vector[], double prob_vec[]);


#endif
