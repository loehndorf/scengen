/*********************************************************************/
/***  Front-end for the HKW scenario generation algorithm          ***/
/***                                                               ***/
/***  author: Michal Kaut                                          ***/
/*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "matrix.h"
#include "HKW_sg.h"


// Parameters that can be changed from command line
static int NmbScen = 4;
static int MaxTrial = 10;
static int HKW_MaxIter = 20;
static int FormatOfTgMoms = 0; // see PrintUsageAndExit for details
static double MaxErrMom = 1e-3;
static double MaxErrCorr = 1e-3;
static int TestLevel = 2;
static char OutFileName[] = "out_scen.txt";
static char ProbFileName[] = "";
static char *StartFileName = NULL;
static unsigned char UseStartDistrib = 0;
static int RandomSeed = 0;


void PrintUsageAndExit(char ExecName[])
{
	printf("\nScenario generation code based on paper by K. Høyland, M. Kaut & S.W. Wallace.\n");
	printf("Code by Michal Kaut (michal.kaut@iot.ntnu.no) & Diego Mathieu.\n");

	printf("\nUsage: %s nmb_scen [options]\n", ExecName);
	printf("\nList of options                                                 %15s\n", "default value");
	printf(" -t trials ... number of Trials (attempts to generate)          %15d\n", MaxTrial);
	printf(" -i iters  ... number of Iterations in HKW alg.                 %15d\n", HKW_MaxIter);
	printf(" -f format ... Format of input data - see below                 %15d\n", FormatOfTgMoms);
	printf(" -m mError ... maximal error for Moments (scaled to var=1)      %15g\n", MaxErrMom);
	printf(" -c cError ... maximal error for Correlations (scaled to var=1) %15g\n", MaxErrCorr);
	printf(" -l levOut ... output Level (amount of output)                  %15d\n", TestLevel);
	printf(" -o oFile  ... Output filename (for scenarios)                  %15s\n", OutFileName);
	printf(" -p pFile  ... read Probabilities from a given file             %15s\n", "p[s]=1/nmb_scen");
	printf(" -s sFile  ... read Starting distribution from a given file     %15s\n", "sampled values");
	printf(" -r rSeed  ... set a random seed                                %15s\n", "use comp. time");
	printf(" -h        ... display an additional Help message\n");

	printf("\n%s:\n  %s\n  %s\n  %s\n  %s\n  %s\n",
		"option -f: format is a sum of following bits (any number from 0 to 16)",
		" 1 -> use population estimators (as in spreadsheets)",
		" 2 -> 2nd moment is Var instead of StDev",
		" 4 -> 4th moment is Kurtosis - 3",
		" 8 -> Higher moments are not scaled by StDev",
		"16 -> Use non-central moments, E{X^i} ... lower bits are ignored"
	);
	exit(1);
}


void PrintHelpAndExit(char ExecName[])
{
	printf("\n");
	printf("\nScenario generation code from paper 'A Heuristic for Moment-matching Scenario");
	printf("\nGeneration' by K. Høyland, M. Kaut & S.W. Wallace; Computational Optimization");
	printf("\nand Applications, 24 (2-3), pp. 169–185, 2003; doi:10.1023/A:1021853807313.\n");
	printf("Code by Michal Kaut (michal.kaut@iot.ntnu.no) & Diego Mathieu.\n");
	printf("\n");
	printf("The code generates scenarios for multivariate random variables.\n");
	printf("Distribution is described by the first four moments and correlations.\n");
	printf("Iterative algorithm -> reports distance from the target properties.\n");
	printf("No convergence guarantee -> runs in several trials, if needed.\n");

	printf("\nINPUT FILES\n");
	printf("target moments must be in file named: 'tg_moms.txt'\n");
	printf("target correlations must be in file: 'tg_corrs.txt'\n");
	printf("These files must include a matrix of numbers in the following format:\n");
	printf(" : number of rows\n");
	printf(" : number of columns\n");
	printf(" : data (by rows)\n");
	printf(" : rest of the file is ignored\n");
	printf("\n");
	printf("If the probabilities are given, the file has to be in a vector format:\n");
	printf(" : number of elements\n");
	printf(" : data\n");
	printf(" : rest of the file is ignored\n");
	printf("\n");
	exit(1);
}


int main(int argc, char* argv[])
{
	int i;
	int errorlevel;

	TMatrix TgMoms=Mat_0, TgCorrs=Mat_0, OutMat=Mat_0;
	TVector Probs=Vec_0;

	// Proces command-line arguments
	if ((argc==2) && (argv[1][0]=='-') && (argv[1][1]=='h')) PrintHelpAndExit(argv[0]);
	if ((argc==1) || (argv[1][0]=='-')) PrintUsageAndExit(argv[0]);
	NmbScen = atoi(argv[1]); // First argument is number of scenarios
	if (NmbScen<4) {
		printf("\nERROR: Wrong number of scenarios\n");
		PrintUsageAndExit(argv[0]);
	}
	for (i=2; i<argc; i=i+2) {
		if ((argv[i][0]!='-') || (argc<=i+1)) PrintUsageAndExit(argv[0]);
		switch (argv[i][1]) {
			case 't': MaxTrial = atoi(argv[i+1]); break;
			case 'i': HKW_MaxIter = atoi(argv[i+1]); break;
			case 'f': FormatOfTgMoms = atoi(argv[i+1]); break;
			case 'm': MaxErrMom = atof(argv[i+1]); break;
			case 'c': MaxErrCorr = atof(argv[i+1]); break;
			case 'l': TestLevel = atoi(argv[i+1]); break;
			case 'o': strcpy(OutFileName, argv[i+1]); break;
			case 'p': strcpy(ProbFileName, argv[i+1]); break;
			case 'r': RandomSeed = atoi(argv[i+1]); break;
			case 's': StartFileName = argv[i+1]; break;
			case 'h': PrintHelpAndExit(argv[0]); break;
			default: PrintUsageAndExit(argv[0]);
		}
	}

	// Read the target files
	if (Mat_GetFromFile(&TgMoms, "tg_moms.txt") > 0)
		exit(1);
	if (Mat_GetFromFile(&TgCorrs, "tg_corrs.txt") > 0)
		exit(1);
	// Check the dimensions
	if ((TgMoms.nrow!=4) | (TgMoms.ncol!=TgCorrs.nrow) | (TgCorrs.nrow!=TgCorrs.ncol)) {
		printf("\n\tERROR - Wrong dimension of input matrices!\n\n");
		exit(1);
	}

	// Allocate the output matrix
	Mat_Init(&OutMat, TgMoms.ncol, NmbScen);
	if (StartFileName != NULL) {
		// read start distribution from file
		// it is TRANSPOSED compared to the natural form!!!
		if (Mat_GetFromFile(&OutMat, StartFileName) > 0)
			exit(1);
  	UseStartDistrib = 1;
	}

	// Set the random seed
	if (RandomSeed > 0)
		srand(RandomSeed);
	else
		srand((int) time(NULL));

	// Create the probabilities
	Vec_Init(&Probs, NmbScen); // Initiate -> if we read it from the file, we can control the size
	if (strlen(ProbFileName) > 0) {
		// read probabilities from file
		if (Vec_GetFromFile(&Probs, ProbFileName) > 0)
			exit(1);
	} else {
		// use default values
		for (i=0; i<NmbScen; i++)
			Probs.val[i] = 1/(double) NmbScen;
	}


	// GENERATION ROUTINE
	errorlevel = HKW_ScenGen(FormatOfTgMoms, &TgMoms, &TgCorrs, &Probs, &OutMat,
	                         MaxErrMom, MaxErrCorr, TestLevel,
	                         MaxTrial, HKW_MaxIter, UseStartDistrib,
	                         NULL, NULL, NULL, NULL);
	if (errorlevel>0) {
		printf("\nWarning: the algorithm did not converge!\n\n");
		return(1);
	}

	printf("\n\tDONE\n\n");
	if (TestLevel > 4)
		Mat_DisplayTransp(&OutMat, "Scenarios");


	{
		FILE *outFile;
		int s;
		// Open the file
		if((outFile = fopen(OutFileName, "w")) == NULL) {
			printf("\n\tERROR: Can not open the output file >%s<!\n\n", OutFileName);
			exit(1);
		}
		fprintf(outFile, "# Scenario generation method from a paper by Hoyland, Kaut & Wallace\n");
		fprintf(outFile, "# Code by Michal Kaut (michal.kaut@iot.ntnu.n) & Diego Mathieu\n\n");
		fprintf(outFile, "# %6s", "prob");
		for (i=0; i<OutMat.nrow; i++)
			fprintf(outFile, "   var %2d", i+1);
		for (s=0; s<NmbScen; s++) {
			fprintf(outFile, "\n%8.6f", Probs.val[s]);
			for (i=0; i<OutMat.nrow; i++)
				fprintf(outFile, " %8.5f", OutMat.val[i][s]);
		}
		fprintf(outFile, "\n");
		fclose(outFile);
	}


	// De-allocate matrices
	Mat_Kill(&TgMoms);
	Mat_Kill(&TgCorrs);
	Mat_Kill(&OutMat);
	Vec_Kill(&Probs);

	return(0);
}
