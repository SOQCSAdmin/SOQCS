/**************************************************************************
* @file dmat.h
* @version 3.7
* @date 9/06/2022
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
    dmatrix();                                              // Create density matrix. The quantity of reserved memory is set by default.
    dmatrix(int i_mem);                                     // Create density matrix. The matrix row number is specified.
    ~dmatrix();                                             // Destroy density matrix
    dmatrix *clone();                                       // Copies a density matrix


    // Basic matrix operations
    double trace();                                         // Calculates trace of the density matrix
    void   normalize();                                     // Normalize to trace=1 density matrix
    void   add(dmatrix *addm);                              // Sum two compatible density matrix
    double fidelity(state* input);                          // Calculate the Fidelity of a density matrix with respect to a state.
    double get_result(mati def,qocircuit *qoc);             //Gets the probability of an event defined by def.
    p_bin *get_pbin();                                      // Returns the diagonal elements of a density matrix as a set of probability bins

    // Update matrix operations
    void add_state(state *newrun, qocircuit *qoc);          // Adds new state to the density matrix according with the detector definitions in qoc
    void sum_state(state *newstate);                        // Adds new state to the density matrix
    void add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc);       // Adds a conditional detection
    void add_state_loss(int ndec, mati def,state *in_state, qocircuit *qoc);      // Adds a conditional detection from sampling results.
    dmatrix *calc_measure(qocircuit *qoc);                  // Adds a conditional detection from sampling results.
    dmatrix *get_counts(int npack,qocircuit *qoc);          // Returns a density matrix independent of the wave packet degrees of freedom
    dmatrix *partial_trace(mati pack_def,qocircuit *qoc);   // Calculates the partial trace of the density matrix (the result is another matrix)

    // Print functions
    void prnt_mtx();                                        // Print the matrix ( numerically )
    void prnt_mtx(qocircuit *qoc);                          // Print the matrix ( human readable form)
    void prnt_mtx(double thresh, qocircuit *qoc);                   // Print the matrix ( human readable form). xcut value overriden by thresh.
    void prnt_results(qocircuit *qoc);                      // Prints the diagonal elements of a density matrix

protected:
    // Auxiliary function. (This ones can not be private).
    void create_dmtx(int i_mem);                                                  // Create density matrix auxiliary function
    int ketcompatible(state* A, state*B,mati pack_idx,qocircuit *qoc);            // Check ket "compatibility"
};
***********************************************************************************/


#include "state.h"

/** @defgroup Dens Density matrix
*   Density matrix classes and methods
*/

const int DEFMATDIM=100;    ///< Default density matrix dimension
const int DEFWIDTH=6;    ///< Default density matrix dimension

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
    *  Creates a density matrix object.
    *  The quantity of reserved memory for the matrix is set by default.
    *  @ingroup Dens_management
    */
    dmatrix();
    /**
    *  Creates a density matrix object where the maximum dimension
    *  of the density matrix is set by its number of rows.
    *  @param int i_mem Row number of the density matrix.
    *  @ingroup Dens_management
    */
    dmatrix(int i_mem);
    /**
    *  Destroys a density matrix object.
    *  @ingroup Dens_management
    */
    ~dmatrix();
    /**
    *  Copies a density matrix object.
    *  @ingroup Dens_management
    */
    dmatrix *clone();

    // Basic matrix operations
    /** @defgroup Dens_basic Density matrix basic operations
    *   @ingroup Dens
    *   Basic operations of the density matrix
    */
    /**
    *  Calculates the trace of the density matrix.
    *  @return Returns the trace of the density matrix.
    *  @ingroup Dens_basic
    */
    double trace();
    /**
    *  Normalizes the density matrix to trace = 1.
    *  @ingroup Dens_basic
    */

    void   normalize();
    /**
    *  Sums two compatible density matrices. In this context compatible means
    *  that they are built with states related to the same circuit ( and with
    *  the same post-selection procedures if any).
    *  @param dmatrix *addm Density matrix to be added.
    *  @ingroup Dens_basic
    */
    void   add(dmatrix *addm);
    /**
    *  Calculates the fidelity of a density matrix with respect a reference state.
    *  @param state *input Reference state to calculate the fidelity of the density matrix.
    *  @return Fidelity value.
    *  @ingroup Dens_basic
    */
    double fidelity(state* input);
    /**
    *  Gets the probability of a counting event defined by def. In other words, we obtain the probability
    *  to find a specific number of photons, in a specific channel, polarization and with a determined frequency
    *  and time characteristics.
    *  @param mati def  Is a 4xn matrix with the definition of the counting event of which the probability wants to be obtained where in each column:<br>
    *  <br>
    *           1st-Row: Is the channels definition. <br>
    *           2nd-Row: Is the polarization definition. <br>
    *           3rd-Row: The wavepacket numbers. To obtain the probability of a counting event in a counter density matrix the packet index must be zero.<br>
    *           4th-Row: The number of photons (or counting events from a detector point of view).<br>
    *  <br>
    *  Example:
    *  @code
        in_term <<  0, 1,
                    H, H,
                    0, 0,
                    3, 3;

        prob=matrix->get_prob(in_term,example);
    * @endcode
    * @param qocircuit *qoc Circuit where we want to perform detection.
    * @return It is returned the probability of the event defined by mati def.
    * @see dmatrix *get_counter(int npack,qocircuit *qoc);
    * @ingroup Dens_basic
    */
    double get_result(mati def,qocircuit *qoc);
    /**
    *  Returns the diagonal elements of a density matrix as a set of probability bins.
    *  @return A set of probability bins with the values of the diaginal elements of the density matrix.
    *  @ingroup Dens_operation
    */
    p_bin *get_pbin();

    // Update matrix operations
    /** @defgroup Dens_update Operations to update the density matrix
    *   @ingroup Dens
    *   Update operations of the density matrix.
    */
    /**
    *  Adds a new state to the density matrix calculation according with detector definitions in qoc.
    *  @param state *newrun State to be added to the set of states used to calculate the density matrix entries.
    *  @param qocircuit *qoc Circuit where the detectors are defined.
    *  @ingroup Dens_update
    */
    void add_state(state *newrun, qocircuit *qoc);
    /**
    *  Adds a new state to the density matrix just summing the new state and increasing the amount of accounted states in one.
    *  This routine is the core of the process of adding states.
    *  @param state *newrun State to be added to the set of states used to calculate the density matrix entries.
    *  @ingroup Dens_update
    */
    void sum_state(state *newstate);
    /**
    *  Adds a new state only if the detection conditions in def are fulfilled. The degrees of freedom where the detection takes places are not stored like
    *  in a post-selection  process. This is and "advanced" form of post-selection where the characteristics
    *  of the detector to carry on the post-selection are considered. This is an auxiliary function, intended only for internal use of the library.
    *  @param mati def Definition of the conditional detection. "def" is a 2xn matrix where in the first
    *  row are defined the channels where the detection takes place while in the second row is the number
    *  of photons to be detected.
    *  @code
        dec_def << 1, 2,
                   1, 1;

        detector_matrix->add_state_cond(dec_def.col(),dec_def,a_state,example);
    *  @endcode
    *
    *  @param int ndec Number of detectors. It has to be the number of columns in def.
    *  @param mati def Detection condition table.
    *  @param state *newrun State to be added to the set of states used to calculate the density matrix entries.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @ingroup Dens_update
    */
    void add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc);
    /**
    *  Adds a new state only if the detection conditions in def are fulfilled. Furthermore it removes the extra loss channels
    *  in the appropriate way to obtain physical measurements. Internally it uses  dmatrix::add_state_cond.
    *  Auxiliary function, intended only for internal use of the library.
    *  @param int ndec Number of detectors. It has to be the number of columns in def.
    *  @param mati def Detection condition table.
    *  @param state *newrun State to be added to the set of states used to calculate the density matrix entries.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @see dmatrix::add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc);
    *  @ingroup Dens_update
    */
    void add_state_loss(int ndec, mati def,state *in_state, qocircuit *qoc);
    /**
    *  Calculate measure. It means to integrate/trace out the frequency degrees of freedom or both the frequency and time degrees of freedom in case there
    *  is no clock associated with the detectors.
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @return A new smaller matrix with the integrated degrees of freedom.
    *  @ingroup Dens_update
    */
    dmatrix *calc_measure(qocircuit *qoc);
    /**
    *  Integrates the frequency and time degrees of freedom to have only the total probability count by channel (and polarization)
    *  @param qocircuit *qoc Circuit to which the density matrix is related.
    *  @return A new smaller matrix with the integrated degrees of freedom.
    *  @ingroup Dens_update
    */
    dmatrix *get_counts(int npack,qocircuit *qoc);
    /**
    *  Calculates the partial trace of the elements of the density matrix. Usually is used to integrate the frequency degree of
    *  freedom to have only the total probability count by channel, polarization and time.
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
    /** @defgroup Dens_print Output information of the density matrix
    *   @ingroup Dens
    *   Methods to print the content of the density matrix.
    */

    /**
    *  Prints a density matrix. States are printed as a list of occupations by level.
    *  @see qocircuit::flag_prnt
    *  @see void state::prnt_ket()
    *  @ingroup Dens_print
    */
    void prnt_mtx();
    /**
    *  Prints a density matrix in human readable from. The base states corresponding with each row are printed in
    * the same format than state::prnt_ket depending on qocircuit::flag_prnt for the given circuit.
    *  @param qocircuit *qoc Circuit to which the density matrix is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void state::prnt_ket()
    *  @ingroup Dens_print
    */
    void prnt_mtx(qocircuit *qoc);
    /**
    *  Prints a relevant density submatrix in human readable from. Only the rows and columns with sum values
    *  greater than "thresh" are printed. This is a way to filter very small entries. The base states of each
    *  row are printed in the same format than state::prnt_ket depending on qocircuit::flag_prnt for the given circuit.
    *  @param double thresh Threshold value to print.
    *  @param qocircuit *qoc    Circuit to which the density matrix is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void state::prnt_ket()
    *  @ingroup Dens_print
    */
    void prnt_mtx(double thresh, qocircuit *qoc);
    /**
    *  Prints the diagonal elements of the density matrix.
    *  @see void state::prnt_ket()
    *  @ingroup Dens_print
    */
    void prnt_results();
    /**
    *  Prints the diagonal elements of the density matrix in human readable form.
    *  @param qocircuit *qoc    Circuit to which the density matrix is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void state::prnt_ket()
    *  @ingroup Dens_print
    */
    void prnt_results(qocircuit *qoc);

protected:
    // Auxiliary functions.
    /** @defgroup Dens_aux Auxiliary protected methods
    *   @ingroup Dens
    *   Auxiliary methods not intended for external use.
    */

    /**
    *  Auxiliary method to create a density matrix.
    *  @param int i_mem Row number of the density matrix.
    *  @ingroup Dens_aux
    */
    void create_dmtx(int i_mem);
    /**
    *  Auxiliary method to check if a density matrix element has an non-zero partial trace.
    *  @param state *A Row dictionary state.
    *  @param state *B Column dictionary state.
    *  @param qocircuit *qoc    Circuit to which the density matrix is related.
    *  @see dmatrix *partial_trace(mati pack_def,qocircuit *qoc);
    *  @ingroup Dens_aux
    */
    int ketcompatible(int A, int B,mati pack_idx,qocircuit *qoc);
};
