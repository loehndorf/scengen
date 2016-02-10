/*********************************************************************/
/***  Høyland-Kaut-Wallace algorithm for scenario generation       ***/
/***                                                               ***/
/***  \author: Michal Kaut                                         ***/
/*********************************************************************/

#ifndef HKW_SG_H
#define HKW_SG_H

#include "matrix.h"
#include "dll_export_def.h"


/// main scenario-generation procedure
/**
	Generates one subtree using the HKW algorithm (heuristics)
	\param[in] FormatOfMoms - format of the target moments
		It is a sum of following bits:
		-  1 -> use population estimators (as in spreadsheets)
		-  2 -> 2nd moment is Var instead of StDev
		-  4 -> 4th moment is Kurtose - 3
		-  8 -> Higher moments are not scaled by StDev
		- 16 -> TarMom[i-1] = E{X^i} ... lower bits are ignored in this case \n
		The first bit should be only combined with the second one.
	\param[in] p_TarMoms pointer to [4 x N] matrix of target moments
	\param[in] p_TgCorrs pointer to [N x N] matrix of target correlations
	\param[in] p_Probs pointer to [S] vector of target probabilities
	\param     p_OutMat pointer to [N x S] matrix where the scenarios go
	\param[in] MaxErrMom maximum allowed error for moments
	\param[in] MaxErrCorr maximum allowed error for correlations
	\param[in] TestLevel level of output during the algorithm
	\param[in] MaxTrial maximum number of trials
	\param[in] MaxIter maximum number of iterations in each trial
	\param[in] UseStartValues use values in p_OutMat as a starting point \n
		This will only be used in the first trial. If this does not converge
		and \a MaxTrial > 1, the other trials start from a random starting point.
	\param[out] p_errMom where to put final error in moments (or NULL)
	\param[out] p_errCorr where to put final error in correlations (or NULL)
	\param[out] p_nmbTrial where to put the number of trials used (or NULL)
	\param[out] p_nmbIter where to put the number of iterations used (or NULL)
	\return zero on convergence, one otherwise

	The main structure of the algorithm is as follows:
	- 1	trial = 1
	- 2	if UseStartValues = true and trial = 1, use the current values in
	   	\a p_OutMat \n otherwise generate a random starting point
	- 3	iter = 1
	- 4	if error-in-moments > \a MaxErrMom, correct moments
	- 5	if error-in-correlations > \a MaxErrCorr, correct correlations
	- 6	if one of the errors > max and iter < \a MaxIter,
	   	then iter++ and go to 4
	- 7	if one of the errors > max and trial < \a MaxTrial,
	   	then trial++ and go to 2
	- 8	if both errors < max, report the results (we have convergence),
		otherwise report failure and return the best scenarios found
**/
DLL_PUBLIC int HKW_ScenGen(int const FormatOfMoms,
                           TMatrix const * const p_TarMoms,
                           TMatrix const * const p_TgCorrs,
                           TVector const * const p_Probs,
                           TMatrix * const p_OutMat,
                           double const MaxErrMom, double const MaxErrCorr,
                           int const TestLevel, int const MaxTrial,
                           int const MaxIter, int const UseStartValues,
                           double * p_errMom, double * p_errCorr,
                           int * p_nmbTrial, int * p_nmbIter);


/// wrapper for HKW_ScenGen(), using only standard C arrays
/**
	For description of parameters see HKW_ScenGen(). \n
	The only difference, apart from the types, is that \a probs is allowed
	to be NULL, meaning equiprobable scenarios, and \a outSc is allowed to
	be NULL, in which case it gets allocated.
**/
DLL_PUBLIC int scengen_HKW(double ** const tgMoms, int const FormatOfMoms,
                           double ** const tgCorrs, double * const probs,
                           int const nVar, int const nScen,
                           double ** const outSc,
                           double const MaxErrMom, double const MaxErrCorr,
                           int const TestLevel, int const MaxTrial,
                           int const MaxIter, int const UseStartValues,
                           double * p_errMom, double * p_errCorr,
                           int * p_nmbTrial, int * p_nmbIter);

#endif

/** \mainpage Moment-matching scenario generation heuristic

	This code implements the moment-matching heuristic from paper
	'<em>A Heuristic for Moment-matching Scenario Generation</em>'
	by Kjetil Høyland, Michal Kaut and Stein W. Wallace, published in
	<em>Computational Optimization and Applications</em>, 24 (2-3), pp. 169–185,
	2003; <a href="http://dx.doi.org/doi:10.1023/A:1021853807313">
	doi:10.1023/A:1021853807313</a>. Most of the code was written by
	Michal Kaut, except for the cubic_solve() function written by Diego Mathieu.

	\section License
	The code is freely distributed under the
	<a href="http://www.eclipse.org/legal/epl-v10.html">
	Eclipse Public License</a>. For information about the license, including
	compatibility with other licenses, see the official
	<a href="http://www.eclipse.org/legal/eplfaq.php">FAQ</a> or its
	<a href="http://en.wikipedia.org/wiki/Eclipse_Public_License">Wikipedia</a>
	entry.

**/
