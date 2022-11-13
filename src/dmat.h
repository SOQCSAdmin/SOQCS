/**************************************************************************
* @file dmat.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Density matrix library
*
***************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************

class dmatrix{
public:
    // Management functions
    dmatrix();                                                    // Create density matrix. The quantity of reserved memory is set by default.
    dmatrix(int i_mem);                                           // Create density matrix. The matrix row number is specified.
    ~dmatrix();                                                   // Destroy density matrix
    dmatrix *clone();                                             // Copies a density matrix
    void clear();                                                 // Clears a density matrix


    // Basic matrix operations
    double trace();                                               // Calculates trace of the density matrix
    void   normalize();                                           // Normalize to trace=1 density matrix
    void   add(dmatrix *addm);                                    // Sum two compatible density matrix
    double fidelity(state* input);                                // Calculate the Fidelity of a density matrix with respect to a state.
    double get_result(mati def,qocircuit *qoc);                   //Gets the probability of an event defined by def.
    p_bin *get_pbin();                                            // Returns the diagonal elements of a density matrix as a set of probability bins

    // Update matrix operations
    void add_state(state *newrun, qocircuit *qoc);                // Adds new state to the density matrix according with the detector definitions in qoc
    void add_state(state *newrun, qodev *dev);                    // Adds new state to the density matrix according with the detector definitions in dev
    void add_reduced_state(int ndec,mati def,veci chlist,state* input, qocircuit *qoc); // Adds a conditional detection and traces out channels
    void add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc);             // Adds a conditional detection
    void sum_state(state *newstate);                              // Adds new state to the density matrix
    dmatrix *calc_measure(qocircuit *qoc);                        // Adds a conditional detection from sampling results.
    dmatrix *calc_measure(qodev *dev);                            // Adds a conditional detection from sampling results (uses qodev instead of qocircuit)
    dmatrix *get_counts(int npack,qocircuit *qoc);                // Returns a density matrix independent of the wave packet degrees of freedom
    dmatrix *partial_trace(mati pack_def,qocircuit *qoc);         // Calculates the partial trace of the density matrix (the result is another matrix)

    // Print functions
    void prnt_mtx();                                              // Print the matrix ( numerically )
    void prnt_mtx(double thresh);                                 // Print the matrix ( numerically ) xcut value overriden by thresh.
    void prnt_mtx(int format,double thresh,qocircuit *qoc);       // Print the matrix ( human readable form)
    void prnt_mtx(int format, double thresh, qodev *dev);         // Print the matrix ( human readable form) using qodev instead of qocircuit
    void aux_prnt_mtx(int format, double thresh, qocircuit *qoc); // Auxiliary method to print matrix
    void prnt_results();                                          // Prints the diagonal elements of a density matrix
    void prnt_results(int format, qocircuit *qoc);                // Prints the diagonal elements of a density matrix

protected:
    // Auxiliary function. (This ones can not be private).
    void create_dmtx(int i_mem);                                                        // Create density matrix auxiliary function
    int ketcompatible(state* A, state*B,mati pack_idx,qocircuit *qoc);                  // Check ket "compatibility"
};
***********************************************************************************/


#include "pbin.h"

const int DEFMATDIM= 100; ///< Default density matrix dimension
const int DEFWIDTH = 6;   ///< Default number of spaces to print matrix entries


/** @defgroup Dens Density matrix
*   Density matrix classes and methods
*/

/** \class dmatrix
*   \brief Contains all the information
*    and methods to create, update, manage and print
*    a density matrix.
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*    @ingroup Dens
*/
class dmatrix{
public:
    // Public variables.
    int N;              ///< Number of states used to create the matrix.
    int mem;            ///< Reserved memory for the matrix dictionaries.

    ket_list *dicc;     ///< Base elements that correspond with each of the rows in the matrix.
    matd     dens;      ///< Density matrix coefficients.


    // Management functions
    /** @defgroup Dens_management Density matrix management
    *   @ingroup Dens
    *   Creation and management of the density matrix object.
    */
    /**
    *  Creates a density matrix object. The quantity of reserved memory for the matrix is set by default.
    *
    *  @ingroup Dens_management
    */
    dmatrix();
    /**
    *  Creates a density matrix object where the maximum dimension of the density matrix is set by its number of rows.
    *
    *  @param int i_mem Row number of the density matrix. (Internal memory).
    *  @ingroup Dens_management
    */
    dmatrix(int i_mem);
    /**
    *  Destroys a density matrix object.
    *
    *  @ingroup Dens_management
    */
    ~dmatrix();
    /**
    *  Copies a density matrix object.
    *
    *  @ingroup Dens_management
    */
    dmatrix *clone();
    /**
    *   Re-initializes a density matrix to zero entries.
    *
    *  @ingroup Dens_management
    */
    void clear();

    // Basic matrix operations
    /** @defgroup Dens_basic Density matrix basic operations
    *   @ingroup Dens
    *   Basic operations of the density matrix.
    */
    /**
    *  Calculates the trace of the density matrix.
    *
    *  @return Returns the trace of the density matrix.
    *  @ingroup Dens_basic
    */
    double trace();
    /**
    *  Normalizes the density matrix to trace = 1.
    *
    *  @ingroup Dens_basic
    */
    void   normalize();
    /**
    *  Sums two compatible density matrices. In this context compatible means
    *  that they are built with states related to the same circuit ( and with
    *  the same post-selection procedures if any).
    *
    *  @param dmatrix *addm Density matrix to be added.
    *  @ingroup Dens_basic
    */
    void   add(dmatrix *addm);
    /**
    *  Calculates the fidelity of a density matrix with respect a reference state.
    *
    *  @param state *input Reference state to calculate the fidelity of the density matrix.
    *  @return Fidelity value.
    *  @ingroup Dens_basic
    */
    double fidelity(state* input);
    /**
    *  Returns the probability of an outcome stored in the density matrix described by the provided definition.
    *
    *  @param mati def  Matrix that defines the outcome. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create a definition depending on the number of rows
    *                            of the matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    * @param qocircuit *qoc Circuit to which the outcomes are referred.
    * @return It is returned the probability of the event defined by mati def.
    * @see get_counts(int npack,qocircuit *qoc);
    * @ingroup Dens_basic
    */
    double get_result(mati def,qocircuit *qoc);
    /**
    *  Returns the diagonal elements of a density matrix as a set of probability bins.
    *
    *  @return A set of probability bins with the values of the diagonal elements of the density matrix.
    *  @ingroup Dens_basic
    */
    p_bin *get_pbin();

    // Update matrix operations
    /** @defgroup Dens_update Density matrix update operations
    *   @ingroup Dens
    *   Update operations of the density matrix.
    */
    /**
    *  Adds a new state to the density matrix using the post-selection condition defined by the circuit detectors.
    *
    *  @param state *newrun State to be added to the density matrix.
    *  @param qocircuit *qoc Circuit where the detectors are defined.
    *  @ingroup Dens_update
    */
    void add_state(state *newrun, qocircuit *qoc);
    /**
    *  Adds a new state to the density matrix using the post-selection condition defined by the device detectors.
    *
    *  @param state *newrun State to be added to the density matrix.
    *  @param qodev *dev Device where the detectors are defined.
    *  @ingroup Dens_update
    */
    void add_state(state *newrun, qodev *dev);
    /**
    *  Adds a new state only if the detector conditions in def are fulfilled. The resulting states are then reduced by tracing out
    *  the channels in chlist. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param int ndec Number of detectors. It has to be the number of columns in def.
    *  @param mati def Definition of the conditional detection. "def" is a 2xn matrix where in the first
    *  row are defined the channels where the conditional detection takes place while in the second row is the number
    *  of photons to be detected.
    *  @param veci chlist. List of channels to be traced out from the states added to the density matrix.
    *  @param state *input State to be added to the density matrix.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @ingroup Dens_update
    */
    void add_reduced_state(int ndec,mati def,veci chlist,state* input, qocircuit *qoc);
    /**
    *  Adds a new state only if the detection conditions in def are fulfilled. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param int ndec Number of detectors. It has to be the number of columns in def.
    *  @param mati def Definition of the conditional detection. "def" is a 2xn matrix where in the first
    *  row are defined the channels where the conditional detection takes place while in the second row is the number
    *  of photons to be detected.
    *  @param state *in_state State to be added to the density matrix.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @ingroup Dens_update
    */
    void add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc);
    /**
    *  Adds a new state to the density matrix just summing the new state and increasing the amount of accounted states in one.
    *  This routine is the core of the process of adding states.
    *
    *  @param state *newstate State to be added to the density matrix.
    *  @ingroup Dens_update
    */
    void sum_state(state *newstate);
    /**
    *  Calculate measure. It means to trace out the frequency degrees of freedom or both the frequency and time degrees of freedom depending on the circuit configuration of the detectors.
    *
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @return A new smaller matrix with the integrated degrees of freedom.
    *  @ingroup Dens_update
    */
    dmatrix *calc_measure(qocircuit *qoc);
    /**
    *  Calculate measure. It means to trace out the frequency degrees of freedom or both the frequency and time degrees of freedom depending on the device configuration of the detectors.
    *
    *  @param qodev *dev Device to which the density matrix is related.
    *  @return A new smaller matrix with the integrated degrees of freedom.
    *  @ingroup Dens_update
    */
    dmatrix *calc_measure(qodev *dev);
    /**
    *  Integrates the frequency and time degrees of freedom to have only the total probability count by channel (and polarization).
    *
    *  @param int npack Number of defined packets. Note that the actual number of defined packets may be inferior to the circuit maximum of them.
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @return A new smaller matrix with the integrated degrees of freedom.
    *  @ingroup Dens_update
    */
    dmatrix *get_counts(int npack,qocircuit *qoc);
    /**
    *  Calculates the partial trace of the elements of the density matrix. Usually is used to integrate the frequency degree of
    *  freedom to have only the total probability count by channel, polarization and time.
    *
    *  @param mati pack_def Integer matrix of three rows and as many columns as photon packets where:<br>
    *               1st-Row: is the packet number.<br>
    *               2nd-Row: is the number assigned to designate the collective set of degrees of freedom not traced out. <br>
    *               3rd-Row: is the number assigned to the collective set of degrees of freedom that are being traced out. <br>
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @return Returns a smaller partially traced density matrix.
    *  @ingroup Dens_update
    */
    dmatrix *partial_trace(mati pack_def,qocircuit *qoc);

    // Print functions
    /** @defgroup Dens_print Density matrix output
    *   @ingroup Dens
    *   Methods to print the content of the density matrix.
    */

    /**
    *  Prints a density matrix. States are printed as a list of occupations by level.
    *
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_mtx();
    /**
    *  Prints a relevant density submatrix. Only the rows and columns with sum values
    *  greater than "thresh" are printed.
    *
    *  @param double thresh Threshold value to print.
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_mtx(double thresh);
    /**
    *  Prints a relevant density submatrix. Only the rows and columns with sum values
    *  greater than "thresh" are printed. ( Circuit version ).
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Threshold value to print.
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_mtx(int format, double thresh, qocircuit *qoc);
    /**
    *  Prints a relevant density submatrix. Only the rows and columns with sum values
    *  greater than "thresh" are printed. ( Device version ).
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Threshold value to print.
    *  @param qodev *dev Device to which the density matrix is related.
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_mtx(int format, double thresh, qodev *dev);
    /**
    *  Prints the diagonal elements of the density matrix.
    *
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_results();
    /**
    *  Prints the diagonal elements of the density matrix.
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param qocircuit *qoc    Circuit to which the density matrix is related.
    *  @see state::prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Dens_print
    */
    void prnt_results(int format, qocircuit *qoc);

protected:
    // Auxiliary functions.
    /** @defgroup Dens_aux Density matrix auxiliary methods
    *   @ingroup Dens
    *   Auxiliary methods not intended for external use.
    */

    /**
    *  Auxiliary method to create a density matrix. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param int i_mem Row number of the density matrix. (Internal memory).
    *  @ingroup Dens_aux
    */
    void create_dmtx(int i_mem);
    /**
    *  Auxiliary method to check if a density matrix element has an non-zero partial trace.<br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state *A Row dictionary state.
    *  @param state *B Column dictionary state.
    *  @param qocircuit *qoc    Circuit to which the density matrix is related.
    *  @see partial_trace(mati pack_def,qocircuit *qoc);
    *  @ingroup Dens_aux
    */
    int ketcompatible(int A, int B,mati pack_idx,qocircuit *qoc);
    /**
    *  Auxiliary method to print density matrix. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket
    *  @param double thresh Threshold value to print.
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @see state::prnt_ket(int format, int column, bool loss, qocircuit *qoc)
    *  @ingroup Dens_aux
    */
    void aux_prnt_mtx(int format, double thresh, qocircuit *qoc);
};
