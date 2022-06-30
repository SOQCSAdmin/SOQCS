//======================================================================================================
// File util.cpp
//
// CONSTANTS and UTILITIES LIBRARY
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "util.h"
#include <limits>
#include <chrono>


// Initialization of extern variables
int base=4;        // Maximum photon occupation by level. (default value 4).


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


    base=nph;
    // Random seed for the random number generator that feeds the distributions
    gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
}


//--------------------------------------------
//
//  Calculate a hash value for a given vector of numbers.
//  (We get an unique number for a given occupation vector)
//
//--------------------------------------------
long long int hashval(int *chainv,int n){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv


    return hashval(chainv,n,base);
}


//--------------------------------------------
//
//  Calculate a hash value for a given vector of numbers in a
//  specified base. (We get an unique number for a given occupation vector)
//
//--------------------------------------------
long long int hashval(int *chainv,int n,int chbase){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv
//  int cnbase        // Base of the numbers in chainv
//  Variables
    long long int value;   // Value to be returned
//  Auxiliary index
    int i;            // Aux index


    value=0;
    for(i=0;i<n;i++){
        if(chainv[i]>=0){
            value=value*(chbase+1)+chainv[i];
        }
    }
    return value;
}


//--------------------------------------------------
//
// Calculates the coupling of two Gaussian wave packets.
//
//--------------------------------------------------
cmplx gauss_coup(double ti,double wi,double dwi, double tj,double wj,double dwj, double tri,double trj){
//  double ti;      // Packet i time
//  double tj;      // Packet j time
//  double wi;      // Packet i frequency
//  double wj;      // Packet j frequency
//  double dwi;     // Packet i width
//  double dwj;     // Packet j width
//  double tri;     // Packet i phase time
//  double trj;     // Packet j phase time
//  Variables
    cmplx dt;       // Difference of time between packets
    cmplx dt2;      // Difference of time squared dt^2
    cmplx dwi2;     // Width of the packet i squared wi^2
    cmplx dwj2;     // Width of the packet i squared wj^2
    cmplx dtr;      // Difference of phase time between packets
    cmplx coef;     // Coefficient/factor
    cmplx ce;       // Exponent
    cmplx result;   // Final result: ce*exp(coef)


    // Calculate derived variables
    dt  = ti-tj;
    dt2 = conj(dt)*dt;
    dwi2= conj(dwi)*dwi;
    dwj2= conj(dwj)*dwj;
    dtr=tri-trj;

    // Expression calculation
    coef=(sqrt(2.0*dwi*dwj))/sqrt(dwi2 + dwj2);
    ce=-(dt2*dwi2*dwj2  + pow((wi - wj),2) + 2.0*jm*dtr*(dwj2*wi + dwi2*wj))/(2.0*(dwi2 + dwj2));
    result=coef*exp(ce);
    if(abs(result)<xcut) result=0.0;

    // Return result
    return result;
}


//--------------------------------------------------
//
// Calculates the coupling of two exponential wave packets.
//
//--------------------------------------------------
cmplx exp_coup(double ti,double wi, double txi, double tj, double wj, double txj,double tri,double trj){
//  double ti;      // Packet i time
//  double tj;      // Packet j time
//  double wi;      // Packet i frequency
//  double wj;      // Packet j frequency
//  double txi;     // Packet i characteristic decay time
//  double txj;     // Packet j characteristic decay time
//  double tri;     // Packet i phase time
//  double trj;     // Packet j phase time
//  Variables
    cmplx dt;       // Difference of time between packets
    cmplx dtr;      // Difference of phase time between packets
    cmplx dw;       // Difference of frequency between packets
    cmplx wm;       // Minimum frequency of the two packets
    int   sgn;      // Do we should return the complex conjugate? 1=Yes/0=No
    cmplx result;   // Final result


    // Calculate derived variables
    if(tj>=ti){
        dt  = (ti/txi)-(tj/txj);
        dw  = txi*wi-txj*wj;
        dtr = (tri/txi)-(trj/txj);
        wm=txi*wi;
        sgn=0;
    }else{
        dt  = (tj/txj)-(ti/txi);
        dw  = txj*wj-txi*wi;
        dtr = (trj/txj)-(tri/txi);
        wm=txj*wj;
        sgn=1;
    }

    // Expression calculation
    result=(exp(0.5*dt - jm*dtr*wm))/(1.0 - jm*dw);
    if(sgn==1) result=conj(result);
    if(abs(result)<xcut) result=0.0;

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
    matc SA;    // Aproximate positive definite coupling matrix
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
    long int num_loops; // Numver of loops
    long int diff;      // old_gray - new_gray
    long int gray_diff; // old_gray^new_grays
    long int bin_index; // Binary index
    cmplx sign;         // Sing value
    cmplx direction;    // Direction of the difference.
    cmplx reduce;       // Rei¡duced row combination
    cmplx *row_comb;    // Row combination
    cmplx *new_vector;  // New vector
    long int gray_diff_index;                   // Gray code index of differences
    thash binary_power_dict;                    // Hash table of binary power differences.
    thash::const_iterator hash_gray_diff_index; // Iterator for the has table of binary power differences.
//  Auxiliary index.
    long int i;         // Aux index
    long int j;         // Aux index.


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
//  Estimation of the confidence we have in a triangular
//  matrix that represents a Grand Schmidt orthonormalization.
//
//-----------------------------------------------
double mat_confidence(matc L){
    int i,j;
    double norm;
    double maxdev;
    double conf;

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

