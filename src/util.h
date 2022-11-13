/**************************************************************************
* @file util.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Constants and utilities library.
*
***************************************************************************/

#include <iostream>      // In/out
#include <iomanip>       // Formatting in/out
#include <complex>       // Complex numbers
#include <random>        // Random number generators
#include <string>        // Strings management
#include <algorithm>     // Permutations
#include <unordered_map> // Hash tables
#include "Eigen/Dense"   // Eigen3 library

// Name spaces
using namespace std;
using namespace Eigen;

//Types definitions.
typedef complex<double>  cmplx;            ///< Complex double precision definition.
typedef MatrixXi  mati;                    ///< Simplified type name of integer matrix.
typedef MatrixXd  matd;                    ///< Simplified type name of double matrix.
typedef MatrixXcd matc;                    ///< Simplified type name of complex double matrix.
typedef VectorXi  veci;                    ///< Simplified type name of integer vector.
typedef VectorXd  vecd;                    ///< Simplified type name of double vector.
typedef VectorXcd vecc;                    ///< Simplified type name of complex double vector.

// Hash table definition
typedef unordered_map<long long int,long long int> thash; ///< Simplified type name for a hash table


//Extern variables.
extern int maxnph;                         ///< Maximum photon occupation by level. (default value 4).
                                           ///< It is also the base of our index system.
// Static variables (a single instance for all the library) Random Number Generator (RNG) variables
// Recipe taken from literature (cpp reference guide, common knowledge).
static std::random_device rd;              ///< Variable to be used to obtain a seed for the random number engine.
static std::mt19937_64 gen(rd());          ///< 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000 generator set up.

// Constant to be used across the library
const double xcut = 1.0e-10;               ///< Value below which a real number is truncated to zero.
const double pi   = std::acos(-1);         ///< Value of Pi.
const std::complex<double> jm(0, 1);       ///< The pure imaginary number i.

// Erfi constants
const int DEFLIMERF = 100;                 /// Number of elements in erfi approximation.
const double cerf   = 0.5*sqrt(pi);        /// Constant used in erfi calculation

// Color codes for input in bash terminal.
// Taken from literature web (common knowledge).
// the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"                  ///< RESET.
#define BLACK   "\033[30m"                 ///< Black.
#define RED     "\033[31m"                 ///< Red.
#define GREEN   "\033[32m"                 ///< Green.
#define YELLOW  "\033[33m"                 ///< Yellow.
#define BLUE    "\033[34m"                 ///< Blue.
#define MAGENTA "\033[35m"                 ///< Magenta.
#define CYAN    "\033[36m"                 ///< Cyan.
#define WHITE   "\033[37m"                 ///< White.
#define BOLDBLACK   "\033[1m\033[30m"      ///< Bold Black.
#define BOLDRED     "\033[1m\033[31m"      ///< Bold Red.
#define BOLDGREEN   "\033[1m\033[32m"      ///< Bold Green.
#define BOLDYELLOW  "\033[1m\033[33m"      ///< Bold Yellow.
#define BOLDBLUE    "\033[1m\033[34m"      ///< Bold Blue.
#define BOLDMAGENTA "\033[1m\033[35m"      ///< Bold Magenta.
#define BOLDCYAN    "\033[1m\033[36m"      ///< Bold Cyan.
#define BOLDWHITE   "\033[1m\033[37m"      ///< Bold White.


//**************************************************************************
//
// Extra functions
//
//*************************************************************************

/**
* Calculation of the power p of an integer x.
*
* @param int x an integer base.
* @param int p an integer exponent.
* @return The power p of an integer x: x^p.
*/
int intpow(int x, unsigned int p);

/**
* Calculation of the factorial of integer n.
*
* @param int n an integer argument.
* @return The factorial of n.
*/
long int factorial(long int n);

/**
* Generates an uniformly distributed random number in the interval between 0 and 1.
*
* @return A random number uniformly distributed between 0 and 1.
*/
double urand();

/**
* Generates a Poisson distributed random number  with most probable average value lambda.
*
* @param double lambda  Most probable outcome.
* @return A random number from a Poisson distribution.
*/
int prand(double lambda);

/**
* Generates a Normally distributed random number of mean value mu and stdev deviation.
*
* @param double mu  Mean value.
* @param double stdev  Standard deviation.
* @return A random number from a Normal distribution.
*/
double grand(double mu, double stdev);

/**
* Sets the maximum occupation by level in the simulation to the value defined by nph.
* The default number of photons if this instruction is not used is four.
*
* @param int nph New maximum number of photons by level.
*/
void cfg_soqcs(int nph);

/**
* Calculation of a hash value for a given occupation vector using maxnph as base.
*
* @param int *chainv  Occupation vector
* @param int  n       Length of the occupation vector chainv.
* @return             A hash index value of chainv.
* \xrefitem know  "KnowIss" "Known Issues" If the number of degrees of freedom is very large a long int may not be enough to represent the return value. Nevertheless
* if the density of indexed elements is small it will also work provided there is no collisions.
*/
long long int hashval(int *chainv,int n);

/**
* Calculation of a hash value for a given occupation vector. The maximum number of photons is specified explicitly to be used as a base.
*
* @param int *chainv  Occupation vector.
* @param int  n       Length of the occupation vector chainv.
* @param int  nph     Maximum number of photons.
* @return             A hash indexing value of chainv.
* \xrefitem know  "KnowIss" "Known Issues" If the number of degrees of freedom is very large a long int may not be enough to represent the return value. Nevertheless
* if the density of indexed elements is small it will also work provided there is no collisions.
*/
long long int hashval(int *chainv,int n,int nph);

/**
* Calculation of a decimal value for a given vector treating each entry as a digit. The most significant digits are those of lower index.
*
* @param int *chainv  Integer vector
* @param int  n       Length of chainv.
* @param int  base    Base in which the numbers of the vector are given.
* @return             A number in decimal base.
*/
long long int decval(int *chainv,int n,int base);


/**
* Calculates the coupling between two Gaussian wave packets of defined parameters.
*
*  @param double ti   Packet i central time.
*  @param double tj   Packet j central time.
*  @param double wi   Packet i central frequency.
*  @param double wj   Packet j central frequency.
*  @param double dwi  Packet i width.
*  @param double dwj  Packet j width.
*  @param double tri  Packet i phase time.
*  @param double trj  Packet j phase time.
*  @return The coupling between two Gaussian wave packets.
*/
cmplx gauss_coup(double ti,double wi,double dwi, double tj,double wj,double dwj, double tri,double trj);

/**
* Calculates the coupling between two exponential wave packets of defined parameters.
*
*  @param double ti   Packet i characteristic time.
*  @param double tj   Packet j characteristic time.
*  @param double wi   Packet i characteristic frequency.
*  @param double wj   Packet j characteristic frequency.
*  @param double txi  Packet i characteristic decay time.
*  @param double txj  Packet j characteristic decay time.
*  @param double tri  Packet i phase time.
*  @param double trj  Packet j phase time.
*  @return The coupling between two exponential wave packets.
*/
cmplx exp_coup(double ti,double wi, double txi, double tj, double wj, double txj,double tri,double trj);

/**
* Gram-Schmidt orthonormalization procedure.
* This method orthonormalizes a set of non orthogonal states using an approximate Cholesky method and returns the transformation matrix between
* the non-orthogonal and the orthogonal set of states:
*
*         |t1>: U11 |nt1>
*         |t2>: U21 |nt1> + U22 |nt2>
*         .
*         .
*         |tn>: Un1 |nt1> + Un2 |nt2> + ... Unn |ntn>
*
*  where |ti> are the elements of the original non orthogonal base and |nti> the new
*  orthonormal ones.
*  as input.
*
* @param matc S   Overlapping coefficients between the non-orthogonal initial states.
* @return         Transformation matrix between the non-orthogonal and orthogonal states ordered by rows as a lower diagonal matrix.
*/
matc GSP(matc S);

/**
* Calculates an estimation of the confidence of the Gram-Schmidt orthonormalized basis.
*
* @param matc L   Triangular matrix result of a Gram-Schmidt orthonormalization procedure.
* @return         Estimation between 1.0 and 0.0 of the correctness of the basis represented by L being 1.0 the best score.
*/
double mat_confidence(matc L);

/**
* Calculates the permanent of a square matrix using the Balasubramanian/Bax/Franklin/Glynn formula
* in gray code (this is a compact and fast bit representation).
*
* @param matc M   Square matrix.
* @return         The permanent of a square complex matrix.
*/
cmplx glynn(matc M);

/**
* String to integer converter that can be initialized with constant strings.
*
* @param const char* str  String.
* @param int h  Base of the conversion.
* @return An integer that is unique to the given string.
*/
constexpr unsigned int str2int(const char* str, int h = 0){
//-----------------------------------------------
//
// String to constant converter
// Solution to use it witch c++ switch reported by Serhiy in
// https://stackoverflow.com/questions/16388510/evaluate-a-string-with-a-switch-in-c
//
//-----------------------------------------------
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

/**
* Method that implements the inverse of the exponential function.
*
* @param double u Function input value.
* @return Inverse of the exponential function for the value of u.
*/
double expi(double u);

/**
* Calculates an approximation to the inverse of the error function as described in: <br>
* https://en.wikipedia.org/wiki/Error_function <br>
* The error function is the cumulated integral of the normal function. The approximation
* to its inverse is created by calculating the coefficients of a series expansion up to the nth power.
* The cut off of the expansion is determined by the constant DEFLIMERF.
*
* @param double u Function input value.
* @return Inverse of the error function for the value u.
*/
double erfi(double u);

