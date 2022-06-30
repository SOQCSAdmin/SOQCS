/**************************************************************************
* @file util.h
* @version 3.7
* @date 19/06/2022
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Constants and utilities library.
* @brief Library of commonly used mathematical constant and functions.
*
***************************************************************************/
#include <iostream>      // In/out
#include <iomanip>       // In/out files
#include <complex>       // Complex numbers
#include <random>        // Random number generators
#include <string>        // Strings management
#include <algorithm>     // Permutations
#include <unordered_map> // Hash tables
#include "Eigen/Dense"   // Eigen library

// Name spaces
using namespace std;
using namespace Eigen;

//Types definitions.
typedef complex<double>  cmplx;           ///< Complex double precision definition.
typedef MatrixXi  mati;                   ///< Simplified type name of integer matrix.
typedef MatrixXd  matd;                   ///< Simplified type name of double matrix.
typedef MatrixXcd matc;                   ///< Simplified type name of complex double matrix.
typedef VectorXi  veci;                   ///< Simplified type name of integer vector.
typedef VectorXd  vecd;                   ///< Simplified type name of double vector.
typedef VectorXcd vecc;                   ///< Simplified type name of complex double vector.

typedef unordered_map<long long int,long long int> thash; ///< Simplified type name for a hash table

// Constant to be used across the library
const double xcut=1.0e-10;                ///< Value below which a real number is considered/truncated to zero.
const double pi = std::acos(-1);          ///< Value of Pi.
const std::complex<double> jm(0, 1);      ///< The pure imaginary number i.


//Static variables.
extern int base;                                        ///< Maximum photon occupation by level. (default value 4).
                                                        ///< It is also the base of our index system (therefore the name).


// Static (a single instance for all the library) Random Number Generator (RNG) variables
// Recipe taken from literature (cpp reference guide, common knowledge).
static std::random_device rd;                           ///< Variable to be used to obtain a seed for the random number engine.
static std::mt19937_64 gen(rd());                       ///< 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000 generator set up.


// Color codes for input in bash terminal.
// Taken form literature web (common knowledge).
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
* @param int x an integer base.
* @param int p an integer exponent.
* @return The power p of an integer x: x^p.
*/
int intpow(int x, unsigned int p);

/**
* Calculation of the factorial of integer n.
* @param int n an integer argument.
* @return The factorial of n.
*/
long int factorial(long int n);

/**
* Calculation of an uniformly distributed random number generated in the interval between 0 and 1.
* The kind of generator is controlled by gen(rd()).
* @see static std::mt19937_64 gen(rd())
* @return A random number uniformly distributed between 0 and 1.
*/
double urand();

/**
* Calculation random number generated from a Poisson distribution  with most probable average value lambda.
* The kind of generator is controlled by gen(rd()).
* @param double lambda  Most probable outcome.
* @see static std::mt19937_64 gen(rd())
* @return A random number from a Poisson distribution.
*/
int    prand(double lambda);

/**
* Calculation random number generated from a Normal distribution.
* The kind of generator is controlled by gen(rd()).
* @param double mu  Mean value.
* @param double stdev  Standard deviation.
* @see static std::mt19937_64 gen(rd())
* @return A random number from a Normal distribution.
*/
double grand(double mu, double stdev);

/**
* Sets the maximum number of photons by level in the simulation
* to the value defined by nph. The default number of photons is this instruction is not used is four.
* This quantity is the base of our index mechanism.
* @param int nph New maximum number of photons by level.
*/
void cfg_soqcs(int nph);

/**
* Calculation of a hash value in b base for a given vector of numbers.
* The base is the number of photons given to the soqcs library by cfg_soqcs.
* If no number has been given a default value of four is taken.
* (We get an unique number for a given occupation vector with b as the maximum
*  occupation by level)
* @param int *chainv  Chain of values, a list of numbers in base b. For us occupations.
* @param int  n       Length of chainv.
* @return             A hash indexing value of chainv
* \xrefitem know "KnowIss" "Known Issues" If the number of degrees of freedom is very large a long int maybe not enough to represent the return value. Nevertheless
* if the density of indexed elements is small (as is our case) it will work provided there is no collisions (same number for two
* representations)
*/
long long int hashval(int *chainv,int n);

/**
* Calculation of a hash value in b base for a given vector of numbers.
* The base is specified explicitly in the call.
* (We get an unique number for a given occupation vector with b as the maximum
*  occupation by level)
* @param int *chainv  Chain of values, a list of numbers in base b. For us occupations.
* @param int  n       Length of chainv.
* @param int chbase   Base of the values in the chain of values.
* @return             A hash indexing value of chainv
* \xrefitem know "KnowIss" "Known Issues" If the number of degrees of freedom is very large a long int maybe not enough to represent the return value. Nevertheless
* if the density of indexed elements is small (as is our case) it will work provided there is no collisions (same number for two
* representations)
*/
long long int hashval(int *chainv,int n,int base);


/**
* Coupling between two Gaussian wave packets of defined parameters.
*  @param double ti   Packet i time.
*  @param double tj   Packet j time.
*  @param double wi   Packet i frequency.
*  @param double wj   Packet j frequency.
*  @param double dwi  Packet i width.
*  @param double dwj  Packet j width.
*  @param double tri  Packet i phase time.
*  @param double trj  Packet j phase time.
*  @return The coupling between two Gaussian wave packets.
*/
cmplx gauss_coup(double ti,double wi,double dwi, double tj,double wj,double dwj, double tri,double trj);

/**
* Coupling between two exponential wave packets of defined parameters.
*  @param double ti   Packet i time.
*  @param double tj   Packet j time.
*  @param double wi   Packet i frequency.
*  @param double wj   Packet j frequency.
*  @param double txi  Packet i characteristic decay time.
*  @param double txj  Packet j characteristic decay time.
*  @param double tri  Packet i phase time.
*  @param double trj  Packet j phase time.
*  @return The coupling between two exponential wave packets.
*/
cmplx exp_coup(double ti,double wi, double txi, double tj, double wj, double txj,double tri,double trj);

/**
* Gram-Schmidt orthonormalization procedure.
* Orthonormalizes a set of non orthogonal states and returns the transformation matrix between
* the non-orthogonal and the orthogonal set of states;
*
*         |t1>: U11 |nt1>
*         |t2>: U21 |nt1> + U22 |nt2>
*         .
*         .
*         |tn>: Un1 |nt1> + Un2 |nt2> + ... Unn |ntn>
*
*  where |ti> are the elements of the original non orthornomal base and |nti> the new
*  orthornormal ones. The overlapping coefficients between the non orthogonal states is provided
*  as input.
* @param matc S   Overlapping coefficients between the non-orthogonal initial states.
* @return         Transformation matrix between the non-orthogonal and orthogonal states ordered by rows as a lower diagonal matrix.
*/
matc GSP(matc S);

/**
* Estimation of the confidence of the Gram-Schmidt orthonormalized basis  obtained by an approximate Cholesky method.
*
* @param matc L   Triangular matrix result of a Gram-Schmidt orthonormalization procedure.
* @return         Estimation between 1.0 and 0.0 of the correctness of the basis represented by L being 1.0 the best score.
*/
double mat_confidence(matc L);

/**
* Calculation of the permanent of a square matrix using the Balasubramanian/Bax/Franklin/Glynn formula
* in gray code (this is a compact and fast bit representation).
* @param matc M   Square matrix
* @return         The permanent of a square complex matrix.
*/
cmplx glynn(matc M);

/**
* String to integer converter that can be initialized with constant strings.
* @param const char* str  String
  @param int h  Base of the conversion
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
