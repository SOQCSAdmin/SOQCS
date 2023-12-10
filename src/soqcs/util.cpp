//======================================================================================================
// File util.cpp
//
// CONSTANTS and UTILITIES LIBRARY
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "util.h"
#include <chrono>


// Initialization of extern variables
int def_nph=4;        // Default value of the maximum photon occupation by level.


//-----------------------------------------------
//
// Calculates the power p of an integer x.
// C++ pow only works with floats.
// Recipe obtained from:
// https://stackoverflow.com/questions/1505675/power-of-an-integer-in-c
// by Matthieu M. and edited by Johan
//
//-----------------------------------------------
int intpow(int x, unsigned int p)
{
//  int x;              // Base
//  unsigned int p;     // Power


  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = intpow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}


//--------------------------------------------
//
//  Calculate the factorial of n
//  Taken from literature/web. Various sources (common knowledge)
//
//--------------------------------------------
long int factorial(long int n){
//  long int n;      // Integer to calculate the factorial


    // Single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}


//--------------------------------------------
//
//  Generate a real random number between 0 and 1.
//
//--------------------------------------------
double urand(){
//  Variables
    uniform_real_distribution<double> rng(0.0, 1.0);    // Uniform distribution


    return rng(gen);
}


//--------------------------------------------
//
//  Generate an integer random number following
//  a Poisson distribution of most probable average
//  value lambda.
//
//--------------------------------------------
int prand(double lambda){
//  double lambda;   // Most probable value
//  Variables
    poisson_distribution<int> rng(lambda);          // Poisson distribution


    return rng(gen);
}


//--------------------------------------------
//
//  Generate an random number following from
//  a Normal/Gaussian distribution.
//
//--------------------------------------------
double grand(double mu, double stdev){
//  double mu;      // Mean value
//  double stdev;   // Standard deviation
//  Variables
    normal_distribution<double> rng(mu,stdev);


    return rng(gen);
}


//--------------------------------------------
//
//  Sets the maximum number of photons by level in the
//  simulation and the debug flag. It establishes a base
//  for the indexing mechanism
//
//--------------------------------------------
void cfg_soqcs(int nph){
//  int nph     // Number of photons by level


    def_nph=nph;
    // Random seed for the random number generator that feeds the distributions
    gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
}


//--------------------------------------------
//
//  Calculate a hash value for a given vector of numbers in a
//  for a maximum number of photons. (We get an unique number for a given occupation vector)
//
//--------------------------------------------
long long int hashval(int *chainv,int n,int nph){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv
//  int nph           // Maximum number of photons


    return decval(chainv,n,nph+1);
}


//--------------------------------------------
//
//  Calculate a decimal value for a given vector of numbers in a
//  specified base.
//
//--------------------------------------------
long long int decval(int *chainv,int n,int base){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv
//  int base          // Base of the numbers in chainv
//  Variables
    long long int value;   // Value to be returned
//  Auxiliary index
    int i;            // Aux index


    value=0;
    for(i=0;i<n;i++){
        if(chainv[i]>=0){
            value=value*base+chainv[i];
        }
    }
    return value;
}

//--------------------------------------------------
//
// Calculates the coupling of two Gaussian wave packets.
//
//--------------------------------------------------
cmplx gauss_coup(double ti,double wi,double dwi, double tj,double wj,double dwj){
//  double ti;      // Packet i central time
//  double tj;      // Packet j central time
//  double wi;      // Packet i central frequency
//  double wj;      // Packet j central frequency
//  double dwi;     // Packet i width
//  double dwj;     // Packet j width
//  Variables
    cmplx dt;       // Difference of time between packets
    cmplx dt2;      // Difference of time squared dt^2
    cmplx dwi2;     // Width of the packet i squared wi^2
    cmplx dwj2;     // Width of the packet i squared wj^2
    cmplx coef;     // Coefficient/factor
    cmplx ce;       // Exponent
    cmplx result;   // Final result: ce*exp(coef)


    // Calculate derived variables
    dt  = ti-tj;
    dt2 = conj(dt)*dt;
    dwi2= conj(dwi)*dwi;
    dwj2= conj(dwj)*dwj;

    // Expression calculation
    coef=(sqrt(2.0*dwi*dwj))/sqrt(dwi2 + dwj2);
    ce=-(dt2*dwi2*dwj2  + pow((wi - wj),2) + 2.0*jm*dt*(dwj2*wi + dwi2*wj))/(2.0*(dwi2 + dwj2));
    result=coef*exp(ce);

    // Return result
    return result;
}


//--------------------------------------------------
//
// Calculates the coupling of two exponential wave packets.
//
//--------------------------------------------------
cmplx exp_coup(double ti,double wi, double txi, double tj, double wj, double txj){
//  double ti;      // Packet i characteristic time
//  double tj;      // Packet j characteristic time
//  double wi;      // Packet i characteristic frequency
//  double wj;      // Packet j characteristic frequency
//  double txi;     // Packet i characteristic decay time
//  double txj;     // Packet j characteristic decay time
//  Variables
    int   sgn;      // Do we should return the complex conjugate? 1=Yes/0=No
    cmplx dt;       // Difference of time between packets
    cmplx dw;       // Difference of frequency between packets
    cmplx txm;      // Minimum characteristic time
    cmplx wm;       // Minimum frequency of the two packets
    cmplx denom;    // Denominator value of the result
    cmplx result;   // Final result


    // Calculate derived variables
    if((tj-ti)>0){
        dt  = tj-ti;
        dw  = wj-wi;
        txm = txi;
        wm  = wi;
        sgn = 0;
    }else{
        dt  = ti-tj;
        dw  = wi-wj;
        txm = txj;
        wm  = wj;
        sgn = 1;
    }

    // Expression calculation
    denom=(txi+txj+2.0*jm*txi*txj*dw)/(2*sqrt(txi*txj));
    result=(exp(-0.5*dt/txm + jm*wm*dt))/denom;
    if(sgn==1) result=conj(result);

    // Return result
    return result;
}

//-----------------------------------------------
//
// Gram-Schmidt orthonormalization procedure by
// means of a modified Cholesky decomposition
//
//-----------------------------------------------
matc GSP(matc S){
// matc C;      // Coupling matrix S.
// Variables
    matc L;     // Cholesky decomposed matrix S=LL^*
    matc D;     // Diagonal matrix S=U D U^*
    matc U;     // Unitary  matrix S=U D U^*
    matc SA;    // Approximate positive definite coupling matrix
    ComplexEigenSolver<matc> es; // Eigensolver
// Auxiliary index
    int i;      // Aux index.


    // Obtain eigenvalues of S
    es.compute(S);
    D = es.eigenvalues().asDiagonal();
    U= es.eigenvectors() ;

    //Fix the negative ones
    for(i=0;i<D.rows();i++){
        if(real(D(i,i))<xcut){
          D(i,i)=min(abs(D(i,i)),xcut);
        }
    }

    // Reconstruct the coupling matrix
    SA=U*D*U.adjoint();

    // Perform Cholesky decomposition
    L = SA.llt().matrixL();

    // Return matrix with the transformation coefficients.
    return L;
}


//-----------------------------------------------
//
//  Calculation of the permanent of a matrix using
//  the Glynn method in Gray code. Translation to C++
//  of the algorithm in python coded by "xnor" in:
//  https://codegolf.stackexchange.com/questions/97060/calculate-the-permanent-as-quickly-as-possible
//
//  Glynn Formula:
//  https://en.wikipedia.org/wiki/Computing_the_permanent
//
//-----------------------------------------------
cmplx glynn(matc M){
//  matc M Square matrix to calculate the permanent.
//  Variables
//  Iteration variables
    long int n;         // Number of row/columns of the square matrix M.
    cmplx total;        // Total value of the permanent.
    long int old_gray;  // Old gray codification
    long int new_gray;  // New gray codification
    long int num_loops; // Number of loops
    long int diff;      // old_gray - new_gray
    long int gray_diff; // old_gray^new_grays
    long int bin_index; // Binary index
    cmplx sign;         // Sing value
    cmplx direction;    // Direction of the difference.
    cmplx reduce;       // Reduced row combination
    cmplx *row_comb;    // Row combination
    cmplx *new_vector;  // New vector
    long int gray_diff_index;                   // Gray code index of differences
    thash binary_power_dict;                    // Hash table of binary power differences.
    thash::const_iterator hash_gray_diff_index; // Iterator for the has table of binary power differences.
//  Auxiliary index.
    long int i;         // Aux index
    long int j;         // Aux index


    // Configuration
    n=M.cols();
    if(n==0) return 1.0;

    // Initializations
    row_comb=new cmplx[n];
    for(j=0;j<n;j++){
        row_comb[j]=0.0;
        for(i=0;i<n;i++){
            row_comb[j]=row_comb[j]+M(i,j);
        }
    }

    total = 0;
    old_gray = 0;
    sign = +1;

    new_vector=new cmplx[n];
    for(i=0;i<n;i++) binary_power_dict[pow(2,i)]=i;
    num_loops=pow(2,n-1);


    //  Main loop
    for(bin_index=1;bin_index<=num_loops;bin_index++){

        reduce=1.0;
        for(i=0;i<n;i++) reduce=reduce*row_comb[i];
        total=total+(sign*reduce);

        new_gray = bin_index^(bin_index/2);
        gray_diff = old_gray^new_gray;
        gray_diff_index = binary_power_dict[gray_diff];
        hash_gray_diff_index=binary_power_dict.find(gray_diff);
        gray_diff_index=hash_gray_diff_index->second;

        for(i=0;i<n;i++) new_vector[i]=M(gray_diff_index,i);
        diff=old_gray-new_gray;
        direction=0;
        if(diff>0) direction=2;
        if(diff<0) direction=-2;

        for(i=0;i<n;i++){
            row_comb[i] = row_comb[i]+ (new_vector[i] * direction);
        }

        sign = -sign;
        old_gray = new_gray;
    }

    // Free memory
    delete[] row_comb;
    delete[] new_vector;

    // Return value
    return total/(double)num_loops;

}


//-----------------------------------------------
//
//  Calculation of the permanent of a matrix using
//  the parallelized Ryser formula as published in:
//  P.H. Lundow, K. Markström, Journal of Computational Physics, Volume 455, 2022,110990.
//  https://doi.org/10.1016/j.jcp.2022.110990.
//
//  The code is a translation into C++ of the Fortran code for the calculation of a
//  permanent of a complex matrix provided with the previous article. Some small modifications
//  for convenience of use have been added.
//
//  Ryser Formula:
//  https://en.wikipedia.org/wiki/Computing_the_permanent
//
//-----------------------------------------------
cmplx ryser_omp( matc M, int nthreads){
//  matc   M;            // Square matrix to calculate the permanent.
//  int    nthreads;     // Number of threads
//  Variables
    int    n;            // Number of rows and columns of the square matrix M
    cmplx  f;            // Permanent value
    cmplx *f_priv;       // Partial values of the permanent calculated by each thread.
    // Auxiliary index
    int    i;            // Aux index

    // Obtain the number of rows
    n=M.rows();

    // Discard trivial solutions
    switch (n)
    {
        case 0: // 0x0 matrix
            f=1.0;
            break;
        case 1: // 1x1 matrix
            f= M(0,0);
            break;
        case 2: // 2x2 matrix
            f = M( 0, 0 )*M( 1, 1 ) +  M( 0, 1 )*M( 1, 0 );
            break;
        case 3: // 3x3 matrix
            f = M( 0, 0 )*M( 1, 1 )*M( 2, 2 )  // 123
              + M( 0, 0 )*M( 1, 2 )*M( 2, 1 )  // 132
              + M( 0, 1 )*M( 1, 0 )*M( 2, 2 )  // 213
              + M( 0, 1 )*M( 1, 2 )*M( 2, 0 )  // 231
              + M( 0, 2 )*M( 1, 0 )*M( 2, 1 )  // 312
              + M( 0, 2 )*M( 1, 1 )*M( 2, 0 ); // 321
            break;
        default: // nxn matrix
            // Reserve memory for the thread results.
            // prama reduce directive does not work with complex numbers.
            f_priv = new cmplx[nthreads]();

            // Create threads
            #pragma omp parallel
            {
                #pragma omp for nowait schedule (static)
                for(i=0;i<nthreads;i++){
                    f_priv[i] = sub_permanent(M, nthreads, i);
                }
            }

            // Join results and calculate the permanent
            f = 0.0;
            for(i=0;i<nthreads;i++) f=f+f_priv[i];
            f = f * 2.0* pow((-1.0),n);

            // Free memory
            delete[] f_priv;
            break;
    }

    // Return permanent
    return f;
}


//-----------------------------------------------
//
//  Auxiliary routine for the calculations of sub-permanents using the Ryser formula.
//  The code is a translation into C++ of the Fortran code for the calculation of a
//  permanent of a complex matrix provided by:
//  P.H. Lundow, K. Markström, Journal of Computational Physics, Volume 455, 2022,110990.
//  https://doi.org/10.1016/j.jcp.2022.110990.
//
//  Ryser Formula:
//  https://en.wikipedia.org/wiki/Computing_the_permanent
//
//-----------------------------------------------
cmplx sub_permanent( matc a, int num, int inx){
//  matc   a;            // Square matrix to calculate the permanent.
//  int    num;          // Number of threads
//  int    inx;          // Threads index.
//  Variables
    int    n;            // Number of rows and columns of a
    int    m;            // Algorithm variable
    int    q;            // Algorithm variable
    int    r;            // Algorithm variable
    int    s;            // Algorithm variable
    int    istart;       // Algorithm variable
    int    icount;       // Algorithm variable
    int    ii;           // Algorithm index
    cmplx  p;            // Algorithm variable
    bool  *x;            // Algorithm vector variable
    cmplx *w;            // Algorithm vector variable
    cmplx  f;            // Sub-permanent result
//  Auxiliary index
    int    i;            // Aux index
    int    j;            // Aux index


    // Init variables
    f=0.0;
    n=a.rows();

    // Work share
    m = pow(2,n-1);
    q = m / num;
    r = m % num;
    istart = inx*q + min (inx, r);
    icount =q;
    if(inx<r) icount=icount+1;
    if(icount==0) return 0.0;

    // Init loop
    s=pow(-1,istart);
    x=unrank_gray(istart,n);
    x[n-1]=false;
    w=new cmplx[n]();
    for(i=0;i<n;i++){
        w[i]=0.0;
        for(j=0;j<n;j++){
            w[i]=w[i]+a(i,j);
        }
        w[i]=a(i,n-1) - 0.5 *w[i];
    }

    for(j=0;j<n;j++){  // Original alfgorith j<n-1 but x[n-1]=false by definition
        if ( x[j] ==true ) for(i=0;i<n;i++) w[i] = w[i] + a( i, j );
    }

    // Main loop
    for(ii=0;ii<icount;ii++){
        p=1.0;
        for(i=0;i<n;i++) p=p*w[i];
        j=0;
        s=-s;
        if( s==1){
           while (not x[j]) j=j+1;
           j=j+1;
        }else{
            p=-p;
        }
        f = f + p;
        x[j]=not x[j];
        if(x[j]==true) for(i=0;i<n;i++) w[i]=w[i]+a(i,j);
        else for(i=0;i<n;i++) w[i]=w[i]-a(i,j);
    }

    // Free memory
    delete[] x;
    delete[] w;

    // Return sub-permanent
    return f;

}


//-----------------------------------------------
//
//  Auxiliary routine that returns  the r:th gray code as a logical array x.
//  The code is a translation into C++ of the corresponding Fortran code provided by the library.
//  P.H. Lundow, K. Markström, Journal of Computational Physics, Volume 455, 2022,110990.
//  https://doi.org/10.1016/j.jcp.2022.110990.
//
//-----------------------------------------------
bool *unrank_gray(int r, int n){
//  int r;      // Input integer
//  int n;      // Vector size
//  Variables
    int s;      // Bit vector
    int auxr;   // Bit mask
    bool *x;    // Boolean vector
//  Auxiliary index
    int i;      // Aux index

    auxr = r << 1;    // Left shift
    s= r^auxr;        // IEOR
    x= new bool[n];
    for(i=0;i<n;i++){ // BTEST
        // It should be i according to original code but it does not work. Some mismatch with respect indexing in C++?
        x[i]= s & (1 << (i+1));
    }

    // Return boolean vector
    return x;
}


//-----------------------------------------------
//
//  Estimation of the confidence we have in a triangular
//  matrix that represents a Grand Schmidt orthonormalization.
//
//-----------------------------------------------
double mat_confidence(matc L){
//  matc L;             // Matrix with the transformation between non-orthogonal to orthonormal states
//  Variables.
    double norm;        // Norm of a row
    double maxdev;      // Deviation from the ideal value one
    double conf;        // Confidence 1-norm between 0 an 1
//  Auxiliary  index
    int    i;           // Aux index
    int    j;           // Aux index


    maxdev=0.0;
    for(i=0;i<L.rows();i++){
        norm=0;
        for(j=0;j<L.cols();j++){
            norm=norm+abs(conj(L(i,j))*L(i,j));
        }
        maxdev=max(maxdev,abs(1.0-norm));
    }

    conf=1-maxdev;
    if(conf<0.0) conf=0.0;
    return conf;
}


//-----------------------------------------------
//
//  Inverse of the exponential function
//
//-----------------------------------------------
double expi(double u){
//  double u;   // Input value of the function.


    return -log(1-u);
}


//-----------------------------------------------
//
//  Approximation of the inverse of the Error
//  function (erfi) as read in:
//  https://en.wikipedia.org/wiki/Error_function
//
//-----------------------------------------------
double erfi( double u){
//  double u;           // Input value of the function.
//  Variables
    double res;         // Approx to erfi(u) up to DEFLIMERF elements
    vecd   c;           // Constants Cm
//  Auxiliary index
    int    k;           // Aux index
    int    m;           // Aux index


    res=0;
    c.resize(DEFLIMERF);
    for(k=0; k<DEFLIMERF;k++){
        c(k)=0.0;
        for(m=0; m<k-1;m++){
            c(k)=c(k)+c(m)*c(k-1-m)/((m+1)*(2*m+1));
        }

        if(k<2) c(k)=1.0;
        res=res+c(k)*pow(cerf*u,2*k+1)/(2*k+1);

    }
    return res;
}
