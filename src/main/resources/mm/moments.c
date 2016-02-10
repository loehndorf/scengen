#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "moments.h"

// Computes a combinatorial number "n over k"
int CombNumber(int const  n, int const k) {
	int i;
	double result = 1;

	for (i=0; i<k; i++)
		result *= (double) (n-i)/(i+1);

	return (int) (result+0.5);
}

// Transforms non-central moments (E[X]^k) to central moments (E[(X-E[X])^k])
// returns mean as the first central moment (instead of zero)
void NonCentral_TO_Central(int const nmb_mom, double const nc_mom[],
                           double c_mom[]) {
	int j,k;

	// !!! nc_mom[0] = c_mom[0] = mean !!!
	c_mom[0] = nc_mom[0];

	for (k=2; k<=nmb_mom; k++) { // moments numbered from 1 (as in the book)!
		c_mom[k-1] = 0;
		for (j=0; j<k; j++) // j=k needs 0-th moment, which we do not have
			c_mom[k-1] += CombNumber(k,j) * pow(-nc_mom[0],j) * nc_mom[k-1-j];
		c_mom[k-1] += CombNumber(k,j) * pow(-nc_mom[0],j); // nc_mom[-1] = 1
	}
	return;
}

// Transforms central moments (E[(X-E[X])^k]) to non-central moments (E[X]^k)
// assumes mean as the first central moment (instead of zero)!!!
void Central_TO_NonCentral(int const nmb_mom, double c_mom[],
                           double nc_mom[]) {
	int j,k;

	// We assume that c_mom[0] includes nc_mom[0] = mean
	nc_mom[0] = c_mom[0];

	// But we need c_mom[0]=0 (that is the definition) in the formula:
	c_mom[0] = 0;

	for (k=2; k<=nmb_mom; k++) { // moments numbered from 1 (as in the book)!
		nc_mom[k-1] = 0;
		for (j=0; j<k; j++) // j=k needs 0-th moment, which we do not have
			nc_mom[k-1] += CombNumber(k,j) * pow(nc_mom[0],j) * c_mom[k-1-j];
		nc_mom[k-1] += CombNumber(k,j) * pow(nc_mom[0],j); // c_mom[-1] = 1
	}

	// return the mean to c_mom[0]
	c_mom[0] = nc_mom[0];

	return;
}

// Normalizes central moments
// i.e. transforms them to zero mean and variance of one
void Normalize_Central(int const nmb_mom, double c_mom[]) {
	int i;

	for (i=2; i<nmb_mom; i++)
		c_mom[i] /= pow(sqrt(c_mom[1]),i+1);

	c_mom[0] = 0;
	c_mom[1] = 1;

	return;
}

// Normalizes non-central moments
// i.e. transforms them to zero mean and variance of one
void Normalize_NonCentral(int const nmb_mom, double nc_mom[]) {
	double *c_mom;
	int i;

	allocate(c_mom, double, nmb_mom);

	NonCentral_TO_Central(nmb_mom, nc_mom, c_mom);
	Normalize_Central(nmb_mom, c_mom);

	// once normalized, central and non-central are equal (E[X]=0)
	for (i=0; i<nmb_mom; i++)
		nc_mom[i] = c_mom[i];

	free(c_mom);
	return;
}


/*
// Computes non-central moments of a mix of normal N(0,1) distributions
// Every distribution is specified by its position (mean) and a probability
void NCMoments_Of_N01_Mix(int const nmb_distr, int nmb_mom,
                          double const position[], double prob[],
                          double nc_mom[]) {
	int d,i;
	double temp;
	double *c_mom;
	double *nc_mom_part;

	allocate(c_mom, double, nmb_mom);
	allocate(nc_mom_part, double, nmb_mom);

	if (nmb_mom > NMOM) {
		printf("\nError - Moments_Of_N01_Mix() can compute only the first");
		printf(" %d moments!", NMOM);
		nmb_mom = NMOM;
	}

	// control probabilities
	temp = 0;
	for (d=0; d<nmb_distr; d++)
		temp += prob[d];
	if (fabs(temp-1) > 1e-10) {
		printf("\nWarning - probabilities of the mix do not sum up to one");
		printf(" - rescaling...\n\n");
		for (d=0; d<nmb_distr; d++)
			prob[d] /= temp;
	}

	// init
	for (i=0; i<nmb_mom; i++) {
		nc_mom[i] = 0;
		c_mom[i] = NormalMom_Central[i];
	}

	// compute non-central moments - just a weighted sum!
	for (d=0; d<nmb_distr; d++) {
		c_mom[0] = position[d]; // c_mom = central moments of N(position[d],1)
		Central_TO_NonCentral(nmb_mom, c_mom, nc_mom_part); // convert to non-c.
		for (i=0; i<nmb_mom; i++)
			nc_mom[i] += prob[d]*nc_mom_part[i];
	}

	free(nc_mom_part);
	free(c_mom);
}

// Computes central moments of a mix of normal N(0,1) distributions
// Every distribution is specified by its position (mean) and a probability
void CMoments_Of_N01_Mix(int const nmb_distr, int const nmb_mom,
                         double const position[], double prob[],
                         double c_mom[]) {
	double * nc_mom;
	allocate(nc_mom, double, nmb_mom);

	NCMoments_Of_N01_Mix(nmb_distr, nmb_mom, position, prob, nc_mom);
	NonCentral_TO_Central(nmb_mom, nc_mom, c_mom);

	free(nc_mom);
}

// Computes both central and non-central moments of a mix of N(0,1) distrib.
// Every distribution is specified by its position (mean) and a probability
void Moments_Of_N01_Mix(int const nmb_distr, int const nmb_mom,
                        double const position[], double prob[],
                        double nc_mom[], double c_mom[]) {
	NCMoments_Of_N01_Mix(nmb_distr, nmb_mom, position, prob, nc_mom);
	NonCentral_TO_Central(nmb_mom, nc_mom, c_mom);
}
*/
