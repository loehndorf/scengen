/*********************************************************************/
/***  Finding parameters of the cubic transformation               ***/
/***  using the Newton algorithm applied to                        ***/
/***  the least-square function                                    ***/
/***                                                               ***/
/*********************************************************************/

/* The data are as such:
	 InMom[0..12] ... input moments + InMom[0]:=1
	 TgMom[0..3]  ... target moments
	 moment[0..4][0..6] ... 'momF(i,j)' used inside the code
	 */

#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "HKW_cubic.h"
#include "misc_macros.h"

void momF1(double value[7],double *xk);
void momF2(double value[7],double *xk);
void momF3(double value[7],double *xk);
static double momF4(double *xk);
double objval( double *xk);
void gradhessian( double c[N],double Q[N][N],double *xk);
int migs (double a[N][N],double x[N][N],int indx[N]);
int elgs (double a[N][N],int indx[N]);
void nextstep(double invhk[N][N],double gk[N],double xk[N]);
double norminf(double vecteur[N]);
int spofa(double mat[N][N]);

double moment[5][7]={{0}};

// function that implements the newton method
double cubic_solve(double *xk){

	const int IncrDevStep = 50; // how often we increase the START_DEV

	int i=0,j,nb_iter=0;
	int indx[N];
	double gk[N];//gradient
	double hessk[N][N], invhk[N][N];//hessian and inverse hessian
	double ni;//infinite norm for the gradient
	double error1,error2;//L2-error on the target moments, at initial and final point
	double x_ini[N]={IDENT1,IDENT2,IDENT3,IDENT4};

	//Initialization of the starting point
	for (j=0;j<N;j++){
		xk[j]=x_ini[j];
	}

	gradhessian(gk,hessk,xk);
	ni=norminf(gk);
	etiq1:
	if (ni>=EPSILON){
		error1=sqrt(objval(xk));
		if (migs(hessk,invhk,indx)>0) return(error1); // ERROR
		for (i=nb_iter;i<NBITMAX;i++){
			nextstep(invhk,gk,xk);//compute the current step dk giving by Hk*dk=-gk
			gradhessian(gk,hessk,xk);
			ni=norminf(gk);
			if (ni<EPSILON){
				break;
			}
			if (migs(hessk,invhk,indx)>0) return(error1); // ERROR

		}
		error2=sqrt(objval(xk));
		if ((i<NBITMAX) && ((error2>error1)|(spofa(hessk)==0)) ){
			//convergence to a local maximizer or to a local minimizer greater than the initial value
			//take a new initial point randomly generated
			for (j=0;j<N;j++){
				xk[j]=x_ini[j]+(frand-0.5)*START_DEV*pow(2,i/IncrDevStep);
			}
			gradhessian(gk,hessk,xk);
			ni=norminf(gk);
			nb_iter=i;
			//re-execute the algorithm
			goto etiq1;
		}
	}
	else{
		//check whether the point is a local minimizer or a local maximizer
		if ((i<NBITMAX) && (spofa(hessk)==0)){//it was a local maximizer
			//take an new initial point randomly generated
			for (j=0;j<N;j++){
				xk[j]=x_ini[j]+(frand-0.5)*START_DEV*pow(2,nb_iter/IncrDevStep);
			}
			gradhessian(gk,hessk,xk);
			ni=norminf(gk);
			//re-execute the algorithm
			goto etiq1;
		}
		else{//local minimizer
			error2=sqrt(objval(xk));
		}
	}
return error2;
}


void momF1(double value[7],double *xk){
	int i,k;
	double c;

	for(i=0;i<7;i++){
			value[i]=0;
	}
	for (k=0;k<N;k++){
		c=xk[k];
		for(i=0;i<7;i++){
			value[i]+=c*InMom[i+k];
		}
	}
}

void momF2(double value[7],double *xk){
	int i,k,l,ml1,ml2;
	double c;

	for(i=0;i<7;i++){
			value[i]=0;
	}
	for (k=0;k<7;k++){
			c=0;
			ml1=min(3,k);
			ml2=max(0,k-3);
			for (l=ml2;l<=ml1;l++){
				c+=xk[l]*xk[k-l];
			}
			for(i=0;i<7;i++){
			value[i]+=c*InMom[k+i];
			}
	}
}

void momF3(double value[7],double *xk){
	int i,k,l,m,ml1,ml2,mm1,mm2;
	double c;

	for(i=0;i<N;i++){
				value[i]=0;
	 }
	for (k=0;k<10;k++){
			c=0;
			ml1=min(3,k);
			ml2=max(0,k-6);
			for(l=ml2;l<=ml1;l++){
				mm1=min(k-l,3);
				mm2=max(k-l-3,0);
				for(m=mm2;m<=mm1;m++){
					c+=xk[k-l-m]*xk[l]*xk[m];
				}
			}
			for(i=0;i<4;i++){
				value[i]+=c*InMom[k+i];
			}
	}
}

static double momF4(double *xk)
{
	double value = 0;
	int k,l,m,n,ml1,ml2,mm1,mm2,mn1,mn2;
	double c;


	for (k=0;k<13;k++){
		c=0;
		ml1=min(3,k);
		ml2=max(0,k-9);
		for(l=ml2;l<=ml1;l++){
			mm1=min(k-l,3);
			mm2=max(k-l-6,0);
			for(m=mm2;m<=mm1;m++){
				mn1=min(k-l-m,3);
				mn2=max(k-l-m-3,0);
				for(n=mn2;n<=mn1;n++){
					c+=xk[k-l-m-n]*xk[l]*xk[m]*xk[n];
				}
			}
		}
		value+=c*InMom[k];
	}
return value;
}

/* Function objval(x) computes and returns obj(x). */
double objval( double *xk)
{
	double value = 0;
	int i;

	for (i=0;i<N;i++)
		value += pow(moment[i+1][0]-TgMom[i],2);
	return value;
}

/* Function gradhessian computes the gradient and the Hessian of obj(x) at the current point xk */
void gradhessian( double c[N],double Q[N][N],double *xk)
{	int i,j,k;

	//first row initialisation
	for(i=0;i<7;i++){
		moment[0][i]=InMom[i];
	}

	//initilisation of the other rows

	momF1(moment[1],xk);
	momF2(moment[2],xk);
	momF3(moment[3],xk);
	moment[4][0]=momF4(xk);

	//gradient
	//initialisation
	for(i=0;i<N;i++){
		c[i]=0;
	}
	//value
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		c[i]+=2*(j+1)*moment[j][i]*(moment[j+1][0]-TgMom[j]);
		}
	}

	// hessian - no initialization needed, done for k=0
	for(i=0;i<N;i++){
		for (j=0; j<=i; j++){
			// updated on 2010-06-17 - gcc complained about indexing error for k=0
			Q[i][j] = 2*(moment[0][i]*moment[0][j]); // version for k=0
			for (k=1; k<N; k++){
				Q[i][j] += 2*(k+1)*((k+1)*moment[k][i]*moment[k][j]
				                    +k*moment[k-1][i+j]*(moment[k+1][0]-TgMom[k]));
			}
			//symetrisation of the Hessian
			if (j<i) {
				Q[j][i]=Q[i][j];
			}

		}
	}
}

/* Updated 10/24/2001. */

/*************************    Program 4.4    ****************************/
/*                                                                      */
/************************************************************************/
/* Please Note:                                                         */
/*                                                                      */
/* (1) This computer program is written by Tao Pang in conjunction with */
/*     his book, "An Introduction to computational Physics," published  */
/*     by cambridge University Press in 1997.                           */
/*                                                                      */
/* (2) No warranties, express or implied, are made for this program.    */
/*                                                                      */
/************************************************************************/

int migs (double a[N][N],double x[N][N],int indx[N]){

/* Function to invert matrix a[][] with the inverse stored
   in x[][] in the output.  copyright (c) Tao Pang 2001. */
  int i,j,k;
  double b[N][N];
//  int elgs();

  for (i = 0; i < N; ++i){
    for (j = 0; j < N; ++j){
      b[i][j] = 0;
    }
  }
  for (i = 0; i < N; ++i){
    b[i][i] = 1;
  }

  if (elgs(a,indx)>0) return(1);

  for (i = 0; i < N-1; ++i){
    for (j = i+1; j < N; ++j){
      for (k = 0; k < N; ++k){
              b[indx[j]][k] = b[indx[j]][k]-a[indx[j]][i]*b[indx[i]][k];
      }
    }
  }

  for (i = 0; i < N; ++i){
    x[N-1][i] = b[indx[N-1]][i]/a[indx[N-1]][N-1];
    for (j = N-2; j >= 0; j = j-1){
      x[j][i] = b[indx[j]][i];
      for (k = j+1; k < N; ++k){
         x[j][i] = x[j][i]-a[indx[j]][k]*x[k][i];
      }
      x[j][i] = x[j][i]/a[indx[j]][j];
    }
  }

  return(0);
}

int elgs (double a[N][N],int indx[N]){

/* Function to perform the partial-pivoting Gaussian elimination.
   a[][] is the original matrix in the input and transformed
   matrix plus the pivoting element ratios below the diagonal
   in the output.  indx[] records the pivoting order.
   copyright (c) Tao Pang 2001. */

  int i, j, k=4, itmp;
  double c1, pi, pi1, pj;
  double c[N];


/* Initialize the index */

  for (i = 0; i < N; ++i){
    indx[i] = i;
  }

/* Find the rescaling factors, one from each row */

  for (i = 0; i < N; ++i){
    c1 = 0;
    for (j = 0; j < N; ++j){
      if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
    }
    c[i] = c1;
  }

/* Search the pivoting (largest) element from each column */

  for (j = 0; j < N-1; ++j){
    pi1 = 0;
    for (i = j; i < N; ++i){
      pi = fabs(a[indx[i]][j])/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

    // ADDED BY MK
    if (k>=N) return (1); // ERROR


/* Interchange the rows via indx[] to record pivoting order */

    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < N; ++i){
      pj = a[indx[i]][j]/a[indx[j]][j];

/* Record pivoting ratios below the diagonal */

      a[indx[i]][j] = pj;

/* Modify other elements accordingly */

      for (k = j+1; k < N; ++k){
        a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
      }
    }
  }

  return (0);
}

/* Function that do a multiplication matrix-vector, give the opposite vector of the result
   and add to the initial xk this direction vector. Used at the beginning of each iteration of
   the newton algorithm */

void nextstep(double invhk[N][N],double gk[N],double xk[N]){
	int i,j;
	double di=0;

	for (i=0;i<N;i++){
		for (j=0;j<N;j++){
			di+=-invhk[i][j]*gk[j];
		}
		xk[i]+=di;
		di=0;
	}
}

//norminf computes the infinite norm of an array
double norminf(double vecteur[N]){
	double val=0;
	int i;
	for(i=0;i<N;i++){
		if (vecteur[i]>val){
			val=vecteur[i];
		}
		else if (-vecteur[i]>val){
				val=-vecteur[i];
		}
	}
	return val;
}

int spofa(double mat[N][N]){
	/* spofa tests if the matrix is positive definite*/
	/*returns 1 if the matrix is positive definite and 0 otherwise*/
	int i,j,k;
	double sum;
	double p[N];

	for(i=0;i<N;i++){
		for(j=i;j<N;j++){
			for (sum=mat[i][j],k=i-1;k>=0;k--){
				sum-=mat[i][k]*mat[j][k];
			}
			if(i==j){
				if(sum<=0.0){
					return 0;
				}
				else{
					p[i]=sqrt(sum);
				}

			}
			else {
				mat[j][i]=sum/p[i];
			}
		}
	}
	return 1; // Cholesky run to the end -> matrix is PD
}
