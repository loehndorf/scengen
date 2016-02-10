#ifndef MOMENTS_H
#define MOMENTS_H

#include "misc_macros.h"


#define NMOM 12
const double NormalMom_Central[12] = {0,1,0,3,0,15,0,105,0,945,0,10395};


/// Computes a combinatorial number "n over k"
int CombNumber(int const n, int const k);

// all the converting procedures were declared as DLL_PUBLIC - why???

/// Transforms non-central moments (E[X]^k) to central moments (E[(X-E[X])^k])
/**
	\param[in] nmb_mom number of moments
	\param[in] nc_mom array of non-central moments (input)
	\param[out] c_mom array of central moments (output)
	\warning returns mean as the first central moment (instead of zero)
**/
void NonCentral_TO_Central(int const nmb_mom, double const nc_mom[],
                           double c_mom[]);

/// Transforms central moments (E[(X-E[X])^k]) to non-central moments (E[X]^k)
/**
	\param[in] nmb_mom number of moments
	\param[in] c_mom array of central moments (input)
	\param[out] nc_mom array of non-central moments (output)
	\warning assumes mean as the first central moment (instead of zero)!
**/
void Central_TO_NonCentral(int const nmb_mom, double c_mom[], double nc_mom[]);

/// Normalizes central moments (i.e. makes mean=0 and variance=1)
/**
	\param[in] nmb_nom number of moments
	\param c_mom the array of moments to normalize (in-place)
**/
void Normalize_Central(int const nmb_mom, double c_mom[]);

/// Normalizes non-central moments (i.e. makes mean=0 and variance=1)
/**
	\param[in] nmb_nom number of moments
	\param nc_mom the array of moments to normalize (in-place)
**/
void Normalize_NonCentral(int const nmb_mom, double nc_mom[]);


/*
// Computes non-central moments of a mix of normal N(0,1) distributions
// Every distribution is specified by its position (mean) and a probability
DLL_PUBLIC void NCMoments_Of_N01_Mix(int const nmb_distr, int nmb_mom,
                                     double const position[], double prob[],
                                     double nc_mom[]);


// Computes central moments of a mix of normal N(0,1) distributions
// Every distribution is specified by its position (mean) and a probability
DLL_PUBLIC void CMoments_Of_N01_Mix(int const nmb_distr, int const nmb_mom,
                                    double const position[], double prob[],
                                    double c_mom[]);


// Computes both central and non-central moments of a mix of N(0,1) distrib.
// Every distribution is specified by its position (mean) and a probability
DLL_PUBLIC void Moments_Of_N01_Mix(int const nmb_distr, int const nmb_mom,
                                   double const position[], double prob[],
                                   double nc_mom[], double c_mom[]);
*/


#endif
