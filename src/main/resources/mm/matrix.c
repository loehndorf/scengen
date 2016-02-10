/*********************************************************************/
/***  Definition of the TVector and TMatrix structures             ***/
/***                                                               ***/
/***  author: Michal Kaut                                          ***/
/*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "misc_macros.h"

#define Ident(i,j) (i==j ? 1 : 0)


int Vec_Init(TVector * const p_vec, int const size)
{
	// Control if already allocated
	if (p_vec->val != NULL) {
		printf("\n\tERROR: Attempt to Initialise vector with val != NULL\n\n");
		exit(1);
	}

	// Allocate memory for the vector
	if ((p_vec->val = (double *) malloc(size * sizeof(double))) == NULL) {
		printf("\n\tERROR: Not Enough Memory To Allocate Vector!\n\n");
		exit(1);
	}
	p_vec->size = size;

	return(0);
}

void Vec_InitAsBigAs(TVector * const p_New, TVector const * const p_Old)
{	Vec_Init(p_New, p_Old->size); }

int Vec_GetFromFile(TVector * const p_A, char const DatFileName[])
{
	int i;
	int size;
	FILE *VecFile;

	// Open the file
	if((VecFile = fopen(DatFileName, "r")) == NULL) {
		printf("\n\tERROR: Can not open datafile >%s<!\n\n", DatFileName);
		exit(1);
	}
	// Read dimensions
	if (fscanf(VecFile, "%d", &size) != 1) {
		printf("\n\tERROR: Not enough data in datafile >%s<!\n\n", DatFileName);
		exit(1);
	}

	// Initialisation / Control of dimensions
	if (p_A->val == NULL) // p_A not initialised -> do it
		Vec_Init(p_A, size);
	else // p_A is initialised -> control dimensions
		if (p_A->size!=size) {
			printf ("\n\tERROR: Tries to read to a vector with wrong size!\n\n");
			exit(1);
		}

	// Read data
	for (i=0; i<size; i++)
		if (fscanf(VecFile, "%lf", &(p_A->val[i])) != 1) {;
			printf("\n\tERROR while reading datafile <%s>!\n\n", DatFileName);
			exit(1);
		}
	fclose(VecFile);

	return(0);
}

void Vec_Kill(TVector * const p_A)
{
	// DeAllocate memory
	if (p_A->val == NULL) return;
	free((void *) p_A->val);
	p_A->val = NULL;
	p_A->size = 0;
}

void Vec_Display(TVector const * const p_A, char const name[])
{
	int i;

	if (p_A->val == NULL) {
		printf("\nWarning: Attemp to print a vector that has not been");
		printf(" initialised!\n\n");
		return;
	}
	printf("Vector %s:\n", name);
	for (i=0; i<p_A->size; i++)
		printf(" %7.4f", p_A->val[i]);
	printf("\n\n");
}

void Vec_DisplayTransp(TVector const * const p_A, char const name[])
{
	int i;

	if (p_A->val == NULL) {
		printf("\nWarning: Attemp to print a vector that is not initialised!\n\n");
		return;
	}
	printf("Vector %s (tr):\n", name);
	for (i=0; i<p_A->size; i++)
		printf(" %7.4f\n", p_A->val[i]);
	printf("\n");
}

void Vec_Copy(TVector const * const p_From, TVector * const p_To)
{
	int i;

	// Test dimensions
	if (p_To->val == NULL) // p_To not initialised -> do it
		Vec_InitAsBigAs(p_To, p_From);
	else // p_To is initialised -> control dimensions
		if (p_To->size!=p_From->size) {
			printf ("\n\tERROR: Wrong size of vectors for copying!\n\n");
			exit(1);
		}

	for (i=0; i<p_From->size; i++)
		p_To->val[i] = p_From->val[i];
}


// ***************************************************************** //

int Mat_Init(TMatrix * const p_mat, int const nmb_rows, int const nmb_cols)
{
	int i;

	// Control if already allocated
	if (p_mat->val != NULL) {
		printf("\n\tERROR: Attempt to Initialise matrix with val != NULL\n\n");
		exit(1);
	}

	// Allocate memory for the matrix
	if ((p_mat->val = (double **) malloc(nmb_rows * sizeof(double *))) == NULL) {
		printf("\n\tERROR: Not Enough Memory To Allocate Matrix!\n\n");
		exit(1);
	}
	for (i=0; i<nmb_rows; i++)
		if ((p_mat->val[i] = (double *) malloc(nmb_cols * sizeof(double)))
		    == NULL) {
			printf("\n\tERROR: Not Enough Memory To Allocate Matrix!\n\n");
			exit(1);
		}
	p_mat->ncol = nmb_cols;
	p_mat->nrow = nmb_rows;

	return(0);
}

void Mat_InitAsBigAs(TMatrix * const p_New, TMatrix const * const p_Old)
{	Mat_Init(p_New, p_Old->nrow, p_Old->ncol); }

int Mat_GetFromFile(TMatrix * const p_A, char const DatFileName[])
{
	int i,j;
	int nmb_rows, nmb_cols;
	FILE *MatFile;

	// Open the file
	if((MatFile = fopen(DatFileName, "r")) == NULL) {
		printf("\n\tERROR: Can not open datafile >%s<!\n\n", DatFileName);
		exit(1);
	}
	// Read dimensions
	if (fscanf(MatFile, "%d %d", &nmb_rows, &nmb_cols) != 2) {
		printf("\n\tERROR: Not enough data in datafile >%s<!\n\n", DatFileName);
		exit(1);
	}

	// Initialisation / Control of dimensions
	if (p_A->val == NULL) // p_A not initialised -> do it
		Mat_Init(p_A, nmb_rows, nmb_cols);
	else // p_A is initialised -> control dimensions
		if ((p_A->nrow!=nmb_rows) | (p_A->ncol!=nmb_cols)) {
			printf ("\n\tERROR: Tries to read to a matrix with wrong size!\n\n");
			printf ("Matrix '%s': expected %dx%d, found %dx%d\n\n",
   			DatFileName, p_A->nrow, p_A->ncol, nmb_rows, nmb_cols);
			exit(1);
		}

	// Read data
	for (i=0; i<p_A->nrow; i++)
		for (j=0; j<p_A->ncol; j++)
			if (fscanf(MatFile, "%lf", &(p_A->val[i][j])) != 1) {;
				printf("\n\tERROR while reading datafile <%s>!\n\n", DatFileName);
				exit(1);
			}
	fclose(MatFile);

	return(0);
}

void Mat_Kill(TMatrix * const p_A)
{
	int i;

	// DeAllocate memory
	if (p_A->val == NULL) return;
	for (i=0; i<p_A->nrow; i++)
	{
		free((void *) p_A->val[i]);
		p_A->val[i] = NULL;
	}
	free((void *) p_A->val);
	p_A->val = NULL;
	p_A->nrow = 0;
	p_A->ncol = 0;
}

void Mat_Display(TMatrix const * const p_A, char const name[])
{
	int i,j;

	if (p_A->val == NULL) {
		printf("\nWarning: Attemp to print a matrix that has not been");
		printf("initialised!\n\n");
		return;
	}
	printf("Matrix %s:\n", name);
	for (i=0; i<p_A->nrow; i++) {
		for(j=0; j<p_A->ncol; j++)
			printf(" %13.10f", p_A->val[i][j]); // was: %7.4f
		printf("\n");
	}
	printf("\n");
}

void Mat_DisplayTransp(TMatrix const * const p_A, char const name[])
{
	int i,j;

	if (p_A->val == NULL) {
		printf("\nWarning: Attemp to print a matrix that has not been");
		printf("initialised!\n\n");
		return;
	}
	printf("Matrix %s (tr):\n", name);
	for (i=0; i<p_A->ncol; i++) {
		for(j=0; j<p_A->nrow; j++)
			printf(" %7.4f", p_A->val[j][i]);
		printf("\n");
	}
	printf("\n");
}

void Mat_Copy(TMatrix const * const p_From, TMatrix * const p_To)
{
	int i,j;

	// Test dimensions
	if (p_To->val == NULL) // p_To not initialised -> do it
		Mat_InitAsBigAs(p_To, p_From);
	else // p_To is initialised -> control dimensions
		if ((p_To->nrow!=p_From->nrow) | (p_To->ncol!=p_From->ncol)) {
			printf ("\n\tERROR: Wrong size of matrices for copying!\n\n");
			exit(1);
		}

	for (i=0; i<p_From->nrow; i++)
		for (j=0; j<p_From->ncol; j++)
			p_To->val[i][j] = p_From->val[i][j];
}

void Mat_CopyTransp(TMatrix const * const p_From, TMatrix * const p_To)
{
	int i,j;

	// Test dimensions
	if (p_To->val == NULL) // p_To not initialised -> do it
		Mat_Init(p_To, p_From->ncol, p_From->nrow);
	else // p_To is initialised -> control dimensions
		if ((p_To->ncol!=p_From->nrow) | (p_To->nrow!=p_From->ncol)) {
			printf ("\n\tERROR: Wrong size of matrices for copy_transp!\n\n");
			exit(1);
		}

	for (i=0; i<p_From->nrow; i++)
	for (j=0; j<p_From->ncol; j++)
		p_To->val[j][i] = p_From->val[i][j];
}

int Mat_AreEqual(TMatrix const  * const p_A, TMatrix const * const p_B)
{
	int i,j;

	if ((p_A->nrow!=p_B->nrow) | (p_A->ncol!=p_B->ncol))
		return(0);
	for (i=0; i<p_A->nrow; i++)
		for (j=0; j<p_A->ncol; j++)
			if (fabs(p_A->val[i][j]-p_B->val[i][j]) > EPS)
				return(0);

	return(1);
}


// Matrix multiplication
int Mat_Mult(TMatrix const * const p_Left, TMatrix const * const p_Right,
             TMatrix * const p_Result)
{
	int i,j,k;

	// Control of dimensions
	if (p_Left->ncol!=p_Right->nrow)
		return(1);
	if (p_Result->val == NULL) // if not initialised, do it
		Mat_Init(p_Result, p_Left->nrow, p_Right->ncol);
	else // if initialised, control dimensions
		if ((p_Left->nrow!=p_Result->nrow) | (p_Right->ncol!=p_Result->ncol))
			return(1);

	for (i=0; i<p_Result->nrow; i++)
		for (j=0; j<p_Result->ncol; j++) {
			p_Result->val[i][j] = 0;
			for (k=0; k<p_Left->ncol; k++)
				p_Result->val[i][j] += p_Left->val[i][k] * p_Right->val[k][j];
		}

	return(0);
}


// Matrix multiplication where the LEFT matrix is TRIANGULAR
int Mat_Mult_LeftTriang(TMatrix const * const p_Left, char WhichTriang,
                        TMatrix const * const p_Right, TMatrix * const p_Result)
{
	int i,j,k;

	if (WhichTriang > 90) WhichTriang -= 32; // convert to upper case
	if ((WhichTriang != 'L') && (WhichTriang != 'U'))
		return(1);

	// Control of dimensions
	if (p_Left->ncol!=p_Right->nrow)
		return(1);
	if (p_Result->val == NULL) // if not initialised, do it
		Mat_Init(p_Result, p_Left->nrow, p_Right->ncol);
	else // if initialised, control dimensions
		if ((p_Left->nrow!=p_Result->nrow) | (p_Right->ncol!=p_Result->ncol))
			return(1);

	if (WhichTriang=='L') { // control it only once!
		for (i=0; i<p_Result->nrow; i++)
			for (j=0; j<p_Result->ncol; j++) {
				p_Result->val[i][j] = 0;
				for (k=0; k<=i; k++)
					p_Result->val[i][j] += p_Left->val[i][k] * p_Right->val[k][j];
			}
	} else {
		for (i=0; i<p_Result->nrow; i++)
			for (j=0; j<p_Result->ncol; j++) {
				p_Result->val[i][j] = 0;
				for (k=i; k<p_Left->ncol; k++)
					p_Result->val[i][j] += p_Left->val[i][k] * p_Right->val[k][j];
			}
	}

	return(0);
}


// Inverse of a LOWER-TRINAGULAR matrix
int Mat_LowTriangInverse(TMatrix const * const p_L, TMatrix * const p_invL)
{
	int i,j,k;
	int n = p_L->ncol;

	int doTest = 1;

	if (p_L->nrow!=n) {
		printf("\n\tERROR : Attemp to invert non-square matrix!\n\n");
		return(1);
	}
	if (p_invL->val == NULL) // target not initialised -> do it
		Mat_Init(p_invL, n,n);
	else // target initialised -> test dimensions
		if ((p_invL->ncol!=n) | (p_invL->nrow!=n)) {
			printf("\n\tERROR: Wrong size of target for inverse matrix!\n\n");
			return(1);
		}

	for (k=0; k<n; k++) {
		for (i=n-1; i>=0; i--) {
			p_invL->val[k][i] = Ident(k,i); // 1's on diagonal
			for (j=i+1; j<n; j++)
				p_invL->val[k][i] -= p_L->val[j][i] * p_invL->val[k][j];
			if (fabs(p_L->val[i][i]) < EPS)
				return(2);
			p_invL->val[k][i] /= p_L->val[i][i];
		}
	}

	if (doTest) {
		TMatrix Test=Mat_0;
		Mat_Init(&Test, n,n);
		Mat_Mult_LeftTriang(p_L, 'L', p_invL, &Test);
		for (i=0; i<n; i++)
			for (j=0; j<n; j++) {
				if (fabs(Test.val[i][j] - (i==j)) > EPS) {
					printf("\nError in the inversion of a matrix!\n\n");
					exit(1);
				}
			}
		Mat_Kill(&Test);
	}

	return(0);
}


/* Cholesky decomposition of a symmetric matrix
	This is the "straightforward" implementation, which should be
	faster and more stable(?) than LU (LDL') decomposition
	*/
int Mat_Cholesky(TMatrix const * const p_R, TMatrix * const p_CholR,
                 int const TestLevel)
{
	int n; // size of the matrices
	int i,j,k;
	double sum;
	int CriticalError = 0; // error indicator
	// test level:
	// 	0 = no testing
	// 	1 = test the final result
	// 	2 = test also intermediate results
	int DoTests = TestLevel / 2; //
	int TestPrint = (TestLevel < 6 ? 0 : 1); // 1 = print tested matrices

	// Test dimensions
	if (p_R->ncol!=p_R->nrow)
		return 1; // not a square matrix
	n = p_R->ncol;
	// Test if symmetric
	for (i=0; i<n; i++)
		for (j=0; j<i; j++)
			if (p_R->val[i][j] != p_R->val[j][i])
				return 2;
	// Test/Set dimensions of the output matrix
	if (p_CholR->val == NULL) // p_CholR not initialised -> do it
		Mat_Init(p_CholR, n,n);
	else // p_CholR is initialised -> control dimensions
		if ((p_CholR->ncol!=p_CholR->nrow) | (p_CholR->ncol!=p_R->ncol))
			return 1;

	for (i=0; i<n; i++) {
		for (j=i; j<n; j++) {
			sum = p_R->val[i][j];
			for (k=0; k<i; k++)
				sum -= p_CholR->val[i][k] * p_CholR->val[j][k];
			if (i == j) {
				if (sum < EPS) {
					return(4);
				}
				p_CholR->val[i][i] = sqrt(sum);
			} else {
				p_CholR->val[j][i] = sum / p_CholR->val[i][i];
			}
		}
		// Fill the rest with zeros
		for (j=i+1; j<n; j++)
			p_CholR->val[i][j] = 0;
	}

	// TEST if  p_R = L L'
	TMatrix Test=Mat_0, U=Mat_0;
	if (DoTests) {
		Mat_CopyTransp(p_CholR, &U);
		Mat_Mult_LeftTriang(p_CholR, 'L', &U, &Test);
		if (!Mat_AreEqual(&Test, p_R)) {
			printf("\n\t!!!  ERROR in the Cholesky Code : R != L L'  !!!\n\n");
			CriticalError = 1;
			TestPrint = 1; // print the matrices
		}
		if (TestPrint) {
			Mat_Display(p_CholR, "L=Chol(R)");
			Mat_Display(&Test, "L L'");
			Mat_Display(p_R, "R");
		}
		Mat_Kill(&Test);
		Mat_Kill(&U);
		if (CriticalError) exit(1);
	}
	return 0;
}
