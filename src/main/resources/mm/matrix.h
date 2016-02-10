/*********************************************************************/
/***  Definition of the TVector and TMatrix structures             ***/
/***                                                               ***/
/***  author: Michal Kaut                                          ***/
/*********************************************************************/

#ifndef MATRIX_H
#define MATRIX_H


/// a simple structure representing a vector of doubles
typedef struct {
	int size;
	double *val;
} TVector;

/// default empty vector
#define Vec_0 {0,NULL}

/// initialize a vector - works only on empty vectors
/**
	\param p_vec pointer to vector to initialize
	\param size size (length) of the vector
	\return zero on success
**/
int Vec_Init(TVector * const p_vec, int const size);

/// initialize one vector to be as big as another one
/**
	\param     p_New pointer to vector to initialize
	\param[in] p_Old pointer to vector we take the size from
**/
void Vec_InitAsBigAs(TVector * const p_New, TVector const * const p_Old);

/// read vector from a given data file
/**
	\param p_A pointer to vector we write the data to
	\param DatFileName file name to read from
	The format is size (int) values (double), the rest is ignored.
	\return zero on success
**/
int Vec_GetFromFile(TVector * const p_A, char const DatFileName[]);

/// delete the vector's data, free the allocated memory
void Vec_Kill(TVector * const p_A);

/// print the vector to stdout, as one line
/**
	\param[in] p_A pointer to vector to display
	\param[in] name string used to identify the vector
**/
void Vec_Display(TVector const * const p_A, char const name[]);

/// print the vector to stdout, as a column (one value per line)
/**
	\param[in] p_A pointer to vector to display
	\param[in] name string used to identify the vector
**/
void Vec_DisplayTransp(TVector const * const p_A, char const name[]);

/// copy one vector to another
/**
	\param[in]  p_From pointer to the source vector
	\param[out] p_To pointer to the target vector
**/
void Vec_Copy(TVector const * const p_From, TVector * const p_To);


// ***************************************************************** //

/// a simple structure representing a matrix of doubles
typedef struct {
	int nrow, ncol;
	double ** val;
} TMatrix;

/// default empty vector
#define Mat_0 {0,0,NULL}

/// initialize a matrix - works only on empty matrices
/**
	\param p_mat pointer to matrix to initialize
	\param nmb_rows number od rows (dim 1)
	\param nmb_cols number od columns (dim 2)
	\return zero on success
**/
int Mat_Init(TMatrix * const p_mat, int const nmb_rows, int const nmb_cols);

/// initialize one matrix to be as big as another one
/**
	\param     p_New pointer to matrix to initialize
	\param[in] p_Old pointer to matrix we take the size from
**/
void Mat_InitAsBigAs(TMatrix * const p_New, TMatrix const * const p_Old);

/// read matrix from a given data file
/**
	\param p_A pointer to matrix we write the data to
	\param DatFileName file name to read from
	The format is dimension (2 x int) values (double), the rest is ignored.
	\return zero on success
**/
int Mat_GetFromFile(TMatrix * const p_A, char const DatFileName[]);

/// delete the matrix's data, free the allocated memory
void Mat_Kill(TMatrix * const p_A);

/// print the matrix to stdout
/**
	\param[in] p_A pointer to the matrix to display
	\param[in] name string used to identify the matrix
**/
void Mat_Display(TMatrix const * const p_A, char const name[]);

/// print the matrix to stdout, transposed
/**
	\param[in] p_A pointer to the matrix to display
	\param[in] name string used to identify the matrix
**/
void Mat_DisplayTransp(TMatrix const * const p_A, char const name[]);

/// copy one matrix to another
/**
	\param[in]  p_From pointer to the source matrix
	\param[out] p_To pointer to the target matrix
**/
void Mat_Copy(TMatrix const * const p_From, TMatrix * const p_To);

/// copy one matrix to another, transposed
/**
	\param[in]  p_From pointer to the source matrix
	\param[out] p_To pointer to the target matrix
**/
void Mat_CopyTransp(TMatrix const * const p_From, TMatrix * const p_To);

/// check if two matrices are equal (all values, up to EPS)
/** \return 1 if equal, 0 otherwise **/
int Mat_AreEqual(TMatrix const * const p_A, TMatrix const * const p_B);


/// matrix multiplication
/**
	\param[in] p_Left pointer to the left input matrix
	\param[in] p_Right pointer to the right input matrix
	\param[out] p_Result pointer to the output matrix
		must be either empty or have correct dimensions
	\return zero on success, 1 otherwise (wrong dimensions)
**/
int Mat_Mult(TMatrix const * const p_Left, TMatrix const * const p_Right,
             TMatrix * const p_Result);

/// matrix multiplication where the LEFT matrix is TRIANGULAR
/**
	\param[in] p_Left pointer to the left (triangular) input matrix
	\param[in] WhichTriang form of triangularity: either 'l'/'L', 'u'/'U'
	\param[in] p_Right pointer to the right input matrix
	\param[out] p_Result pointer to the output matrix
	\return zero on success, 1 otherwise (wrong dimensions etc)
	\warning Does not control triangularity, just ignores the rest
	 */
int Mat_Mult_LeftTriang(TMatrix const * const p_Left, char WhichTriang,
                        TMatrix const * const p_Right,
                        TMatrix * const p_Result);


/// inverse of a LOWER-TRINAGULAR matrix
/**
	\param[in] p_L pointer to the matrix to invert
	\param[out] p_invL pointer to the inverted matrix
	\return zero on success, 1 otherwise (wrong dimensions etc)
**/
int Mat_LowTriangInverse(TMatrix const * const p_L, TMatrix * const p_invL);

/// Cholesky decomposition of a symmetric matrix
/**
	Some of the routines are taken from "Num Math and Computing"
	alg. is: R = L U = L D L' = CholR CholR'
	\param[in] p_R pointer to the matrix to decompose
	\param[out] p_CholR pointer to the matrix with the result (lower-triang.)
	\param[in] TestLevel controls the amount of testing and output
	\return zero on success, an error-code otherwise
**/
int Mat_Cholesky(TMatrix const * const p_R, TMatrix * const p_CholR,
                 int const TestLevel);

#endif
