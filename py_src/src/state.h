/**************************************************************************
* @file state.h
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
* @title State definition library
* @brief In this library we define and manage the definition of a bosonic
*        quantum state
*
***************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************
class ket_list{
public:
    // Public methods
    // Management methods
    ket_list(int i_level);                                          //  Creates a ket list.The maximum number of kets is set by default.
    ket_list(int i_level,int i_maxket);                             //  Creates a ket list specifying the maximum number of kets.
    ket_list(int i_level, int i_maxket, int *i_vis);                //  Creates a ket list specifying the maximum number of kets and a vector of equivalence between state and circuit levels.
    ~ket_list();                                                    //  Destroys a ket_list
    void clear_kets();                                              //  Clear the list (without destroying it).

    // State manipulation methods.
    int add_ket(int *occ);                                          // Adds a new ket to a ket list
    int add_ket(hterm term, qocircuit *qoc);                        // Adds a new ket to a ket list ( human readable form)
    int find_ket(int *occ);                                         // Finds the position of a ket in the list given and occupation description
    int find_ket(mati def,qocircuit *qoc);                          // Finds the position of a ket described in human readable form
    ket_list *remove_time(qocircuit *qoc);                          // Remove levels with wave-packet indexes greater than zero. Needed for some specific operations

    //Print methods
    void  prnt_ket(int iket,qocircuit *qoc);                        // Prints a ket
    void  prnt_ket(int iket, bool loss,qocircuit *qoc);             // Prints a ket


protected:
    void create_ket_list(int i_level, int i_maxket);                // Create ket list auxiliary function
};


class state{
public:
    // Public methods
    // Management methods
    state(int i_level);                                             //  Creates a state.The maximum number of kets is set by default.
    state(int i_level,int i_maxket);                                //  Creates a state specifying the maximum number of kets.
    state(int i_level, int i_maxket, int *i_vis);                   //  Creates a state specifying the maximum number of kets and a vector of equivalence between state and circuit levels.

    ~state();                                                       //  Destroys a state
    state *clone();                                                 //  Copy a state
    void clear_state();                                             //  Clears a state

    // State manipulation methods.
    void add_term(cmplx i_ampl, int *occ);                          // Adds a new term to a state
    void add_term(cmplx i_ampl, hterm term, qocircuit *qoc);        // Adds a new term to a state ( human readable form)
    void dproduct(state *rhs);                                      // Direct product of states defined in non coincident channels).
    cmplx braket(state *bra);                                       // Calculates the braket <bra|state>
    void normalize();                                               // Normalizes the state
    state *post_selection(state *prj);                              // Post select a state using a "projector"
    state *prepare_reference(veci ch, int t,qocircuit *qoc);        // Remove all the levels corresponding to the specified channels provided that they are zero.

    //Print methods
    void  prnt_state();                                             // Prints a state ( various kets and amplitudes )
    void  prnt_state(qocircuit *qoc, int column);                   // Prints a state ( in human readable form )
    void  prnt_state(qocircuit *qoc, bool loss,int column);         // Prints a state ( in human readable form )

    // Emitters/Initial states
    void QD(mati ch, double k, double S, double tss, double thv, qocircuit *qoc);                                       // QD state generator model
    void QDPair(int i_ch0,int i_ch1, veci i_t, double dt, double k, double S, double tss, double thv, qocircuit *qoc);  // Single pair photon emission in a QD.
    void Bell(int i_ch0,int i_ch1, veci i_t, double phi,char kind, qocircuit *qoc);                                     // Non-ideal bell emitter with a phase e^(-i phi) in the second term.
    void Corr(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);                                                           // Correlated emitter
    void Rand(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);                                                           // Random emitter

protected:
    // Auxiliary methods
    void prnt_in_rows(qocircuit *qoc, bool loss);                   // Auxiliary method to print a state in rows
    void prnt_in_cols(qocircuit *qoc, bool loss);                   // Auxiliary method to print a state in columns
};

class projector : public state{
public:
    // Public methods
    // Management methods
    projector(int i_level);                                         //  Creates a projector.The maximum number of kets is set by default.
    projector(int i_level,int i_maxket);                            //  Creates a projector specifying the maximum number of kets.
    projector(int i_level, int i_maxket, int *i_vis);               //  Creates a projector specifying the maximum number of kets and a vector of equivalence between state and circuit levels.

    // Auxiliary methods
    void create_projector(int i_level, int i_maxket);               // Create projector auxiliary function
};


class p_bin: public ket_list{
public:
    // Public methods
    // Management methods
    p_bin(int i_level);                                             //  Creates a set of probability bins.The maximum number of outcomes is set by default
    p_bin(int i_level,int i_maxket);                                //  Creates a set of probability bins specifying the maximum number of outcomes
    p_bin(int i_level, int i_maxket, int *i_vis);                   //  Creates a set of probability bins specifying the maximum number of outcomes and a vector of equivalence between state and circuit levels.
    ~p_bin();                                                       //  Destroys a set of probability bins
    p_bin *clone();                                                 //  Copies a set of probability bins

    // Bin update methods
    int add_count(int *occ);                                        // Counts a new sample
    void add_bin(p_bin *input);                                     // Adds the statistics from another bin
    void add_state(state *input);                                   // Adds the statistics of a quantum state

    // Bin manipulation methods
    p_bin *calc_measure(qocircuit *qoc);                            // Calculates a measurement as described by the detectors in the circuit
    p_bin *post_selection(state *prj);                              // Perform post-selection over the bins given a projector
    p_bin *dark_counts(int S, qocircuit* qoc);                      // Computes dark counts effects
    p_bin *blink(int S, qocircuit* qoc);                            // Computes detector dead time effects
    p_bin *compute_loss(int ndec, mati def,qocircuit *qoc);         // Remove the additional channels to compute losses
    p_bin *compute_cond(int ndec, mati def,qocircuit *qoc);         // Applies the detection conditions
    p_bin *remove_freq(qocircuit* qoc);                             // Remove frequencies and leaves only times
    p_bin *perform_count(qocircuit* qoc);                           // Sum all the contributions to a channel independently of time or frequency
    p_bin *remove_time(qocircuit *qoc);                             // Rome the all the packet definitions that are not 0
    p_bin *white_noise(double stdev);                               // Adds Gaussian white noise

    //Bin consultation methods
    string tag(int index);                                          // Returns the occupation of the indexed bin in string format
    double prob(int index);                                         // Probability value of the indexed bin
    double prob(mati def,qocircuit *qoc);                           // Gets the probability of a bin defined by def.

    // Print bins
    void  prnt_bins();                                              // Prints the bin list ( various states and probabilities )
    void  prnt_bins(qocircuit *qoc, bool thresh);                   // Prints the bin list ( in human readable form )
    void  prnt_bins(qocircuit *qoc, double thresh, bool loss);      // Prints the bin list with a threshold limit
};

class ph_bunch{
public:
    // Public methods
    // Management methods
    ph_bunch(int i_level);                                          // Create a photon bunch
    ph_bunch(int i_level,int i_maxket);                             // Create a photon bunch
    ~ph_bunch();                                                    // Destroys a photon bunch object
    void clear();                                                   // Clears a photon bunch

    // Manipulation methods
    void add_photons(int N, int ch, qocircuit *qoc);                // Add photons
    void add_photons(int N, int ch, int P, double t, double f, double w, qocircuit *qoc);  // Add photons with physical description of the photon shape
    void send2circuit(qocircuit *qoc);                              // Send photons to the circuit
    void send2circuit(char ckind, int rand, qocircuit *qoc);        // Send photons to the circuit specifying the photon shape model
    void weight(cmplx A);                                           // Amplitude of the bunch if alternatices/superposition of state are considered
    void alternative();                                             // Superposition of states
    state *bunch_state();                                           // Return state of the photons

    // Print methods
    void  prnt_state();                                             // Prints a state ( various kets and amplitudes )
    void  prnt_state(qocircuit *qoc, int column);                   // Prints a state ( in human readable form )
    void  prnt_state(qocircuit *qoc, bool loss,int column);         // Prints a state ( in human readable form )
protected:
    // Auxiliary methods
    void create_ph_bunch(int i_level, int i_maxket);                // Auxiliary method to create photon bunches.
};

***********************************************************************************/


#include "qocircuit.h"
// #include <unordered_map>


// Alias
// typedef unordered_map<int,int> thash; ///< Simplified type name for a hash table
typedef MatrixXi hterm;               ///< Type "Human Term" (hterm). Matrix definition for a state in human readable terms.


// Constant defaults
const int DEFSTATEDIM=50;             ///< Default density matrix dimension
const double DEFTHOLDPRNT=0.0001;     ///< Default density matrix dimension

//Group definitions
/** @defgroup Ket_List Ket list
 *  Ket list states and methods
 */

/** \class ket_list
*   \brief Contains all the information
*          to create and manipulate a list of kets defined as levels and occupations
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup Ket_List
*
*/
class ket_list{
public:
    // Public variables
    int nket;              ///< Number of kets (C1|1>+C2|2>+...+Cn|nket>.
    int maxket;            ///< Maximum number of kets.
    int nlevel;            ///< Number of levels in each ket |0, 1, 2, ... nlevel>.

    // Ket list definition
    thash ketindex;        ///< Hash table of the dynamic base dictionary
    int **ket;             ///< Ket definitions. Level occupations of each ket/term
    int *vis;              ///< Print visibility vector.
                           ///< It stores to which level correspond each vector position.
                           ///< After post-selection it keeps track of the original level number instead of relying in order for printing.

    // Public methods
    // Management methods
    /** @defgroup Ket_management Ket List management
    *   @ingroup Ket_List
    *   Creation and management of ket list object.
    */

    /**
    *  Creates a ket_list object. The maximum number of kets is set by default.
    *  @param int i_level   Number of levels to describe a ket.
    *  @ingroup Ket_management
    */
    ket_list(int i_level);
    /**
    *  Creates a ket_list object specifying the maximum number of kets.
    *  @param int i_level   Number of levels to describe a ket.
    *  @param int i_maxket  Maximum number of different kets in the list.
    *  @ingroup Ket_management
    */
    ket_list(int i_level,int i_maxket);
    /**
    *  Creates a ket_list object specifying the maximum number of kets
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. Intended for internal use.
    *  @param int i_level   Number of levels to describe a ket.
    *  @param int i_maxket  Maximum number of different kets in the list.
    *  @param int i_vis  Vector that translates the ket levels to the circuit levels defined by the circuit indices.
    *  @ingroup Ket_management
    */
    ket_list(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys a ket_list object
    *  @ingroup Ket_management
    */
    ~ket_list();
    /**
    *   Creates a copy of this ket list.
    *   @returns    A copy of the present ket list.
    *   @ingroup Ket_management
    */
    ket_list *clone();
    /**
    *   Empties the ket list.
    *  @ingroup Ket_management
    */
    void clear_kets();

    // State manipulation methods.
    /** @defgroup Ket_operations Ket list operations
    *   @ingroup Ket_List
    *   List of operations of a ket list.
    */

    /**
    *   Adds a new ket to the list.
    *  @param int   *occ    List with the occupation of each level in the new ket.
    *  @ingroup Ket_operations
    */
    int add_ket(int *occ);
    /**
    *   Adds a new ket to the list (from a human readable configuration).
    *  @param hterm      term    Matrix that defines the new term. Each column defines the configuration of one level.<br>
    *                            There are four different possible kind of configurations depending on the number of rows
    *                            of the term matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavefunction/packet and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. The wavefunciont/time is assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavefinction are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in order.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore every configuration not mentioned is assumed to be zero.
    *
    *  @param qocircuit *qoc     Circuit to which the ket is related (to translate from human to numeric form).
    *  @ingroup Ket_operations
    */
    int add_ket(hterm term, qocircuit *qoc);
    /**
    *   Finds the position of a ket in the list.
    *  @param int   *occ    List with the occupation of each level in the ket we are searching.
    *  @return Position of the ket in the list.
    *  @ingroup Ket_operations
    */
    int find_ket(int *occ);
    /**
    *  Finds the position of a ket in the list (from a human readable configuration).
    *  @param hterm      term    Matrix that defines the searched term. Each column defines the configuration of one level.<br>
    *                            There are four different possible kind of configurations depending on the number of rows
    *                            of the term matrix:<br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavefunction/packet and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. The wavefunciont/time is assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavefinction are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in order.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore every configuration not mentioned is assumed to be zero.
    *
    *  @param qocircuit *qoc     Circuit to which the ket is related (to translate from human to numeric form).
    *  @return Position of the ket in the list.
    *  @ingroup Ket_operations
    */
    int find_ket(hterm def,qocircuit *qoc);
    /**
    *  Removes levels with a wave-packet index grater than zero.
    *  @param qocircuit *qoc     Circuit to which the ket list is related.
    *  @ingroup Ket_operations
    */
    ket_list *remove_time(qocircuit *qoc);

    //Print methods
    /** @defgroup Ket_output Ket list output
    *   @ingroup Ket_List
    *   Methods to print a ket from the list on screen.
    */

    /**
    *  Prints the ket "#iket" of the list.
    *  Note that the print format is different depending on the circuit print flag value.<br>
    *  <br>
    *  qocircuit::flag_prnt = 0 The ket is printed in numerical format: | (channel, polarization, wavefunction): occupation >. <br>
    *  Detailed format:
    *  @code
        |(ch1{, m1}{, t1}): o1,  ...,(chn{, mn}{, tn}): on >
    *  @endcode
    *  Example:
    *  @code
        | (0, 0, 0): 0, (0, 0, 1): 0, (0, 1, 0): 2, (0, 1, 1): 0 >
    *  @endcode
    *  <br>
    *  qocircuit::flag_prnt = 1 The polarization is referred by a letter H/V: | (channel, H/V, wavefunction): occupation >. <br>
    *| Detailed format:
    *  @code
        | (ch1{, H/V}{, t1}): o1, ...,(chn{, H/V}{, tn}): on >
    *  @endcode
    *  Example:
    *  @code
        | (0, H, 0): 0, (0, H, 1): 0, (0, V, 0): 2, (0, V, 1): 0 >
    *  @endcode
    *  <br>
    *  qocircuit::flag_prnt = 2  Human readable form: | [occupation]H/V(wavefunction)channel >  <br>
    *  Detailed format:
    *  @code
        {{ | {{[o1]}o1>1}{H/V}{(t1)}ch1 > } o1>0},...,{{ | {{[on]}o1>1}{{H/V}{(tn)}chn > } on>0}
    *  @endcode
    *  Example:
    *  @code
        | [2]V(0)0 >
    *  @endcode
    *  <br>
    *  If the pointer to the circuit "qocircuit* qoc" is Null then the occupation values are printed in ascending order. In this
    *  case is responsibility of the user to be aware to which circuit configuration corresponds with each level.<br>
    *  Example:<br>
    *  @code
        | 0, 0, 2, 0 >
    *  @endcode

    *  <br>
    *  Legend:<br>
    *  chi:Channel i, mi: polarization i, ti: wavefunction i. <br>
    *  {}   Optional values only to be printed if there is more than one possibility. <br>
    *  {{}v>i} Optional values only to be printed if the variable v>i.<br>
    *  @param int iket           Number of ket to be printed.
    *  @param qocircuit *qoc Circuit to which the ket is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void qocircuit::set_prnt_flag(int flag);
    *  @ingroup Ket_output
    */
    void  prnt_ket(int iket,qocircuit *qoc);
    /**
    *  Prints the ket "#iket" of the list.
    *  Note that the print format is different depending on the circuit print flag value.
    *  If loss=True loss channels are printed in blue.<br>
    *  <br>
    *  qocircuit::flag_prnt = 0 The ket is printed in numerical format: | (channel, polarization, wavefunction): occupation >. <br>
    *  Detailed format:
    *  @code
        |(ch1{, m1}{, t1}): o1,  ...,(chn{, mn}{, tn}): on >
    *  @endcode
    *  Example:
    *  @code
        | (0, 0, 0): 0, (0, 0, 1): 0, (0, 1, 0): 2, (0, 1, 1): 0 >
    *  @endcode
    *  <br>
    *  qocircuit::flag_prnt = 1 The polarization is referred by a letter H/V: | (channel, H/V, wavefunction): occupation >. <br>
    *| Detailed format:
    *  @code
        | (ch1{, H/V}{, t1}): o1, ...,(chn{, H/V}{, tn}): on >
    *  @endcode
    *  Example:
    *  @code
        | (0, H, 0): 0, (0, H, 1): 0, (0, V, 0): 2, (0, V, 1): 0 >
    *  @endcode
    *  <br>
    *  qocircuit::flag_prnt = 2  Human readable form: | [occupation]H/V(wavefunction)channel >  <br>
    *  Detailed format:
    *  @code
        {{ | {{[o1]}o1>1}{H/V}{(t1)}ch1 > } o1>0},...,{{ | {{[on]}o1>1}{{H/V}{(tn)}chn > } on>0}
    *  @endcode
    *  Example:
    *  @code
        | [2]V(0)0 >
    *  @endcode
    *  <br>
    *  If the pointer to the circuit "qocircuit* qoc" is Null then the occupation values are printed in ascending order. In this
    *  case is responsibility of the user to be aware to which circuit configuration corresponds with each level.<br>
    *  Example:<br>
    *  @code
        | 0, 0, 2, 0 >
    *  @endcode

    *  <br>
    *  Legend:<br>
    *  chi:Channel i, mi: polarization i, ti: wavefunction i. <br>
    *  {}   Optional values only to be printed if there is more than one possibility. <br>
    *  {{}v>i} Optional values only to be printed if the variable v>i.<br>
    *  @param int iket   Number of ket to be printed.
    *  @param bool loss  Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the ket is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void qocircuit::set_prnt_flag(int flag);
    *  @ingroup Ket_output
    */
    void  prnt_ket(int iket, bool loss,qocircuit *qoc);


protected:

    /**
    *  Auxiliary method to create a ket list.
    *  Each ket is described as level occupations.
    *  @param int i_level   Number of levels to describe a ket.
    *  @param int i_maxket  Maximum number of different kets in the list.
    */
    void create_ket_list(int i_level, int i_maxket);

};

/** @defgroup State
 *  State classes and methods
 */

/** \class state
*   \brief Contains all the information
*          to create and manipulate a state as a list of
*          levels and occupation.
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup State
*
*/
class state: public ket_list{
public:
    // Public variables
    // State definition
    // A state it is described as a summation of kets
    // multiplied by amplitudes.
    complex<double> *ampl; ///< Amplitudes of each ket/term

    // Public methods
    // Management methods
    /** @defgroup State_management State management
    *   @ingroup State
    *   Creation and management of state object.
    */

    /**
    *  Creates a state object. The maximum number of kets is set by default.
    *  @param int i_level   Number of levels to describe the state.
    *  @ingroup State_management
    */
    state(int i_level);
    /**
    *  Creates a state object specifying the maximum number of kets.
    *  @param int i_level   Number of levels to describe the state.
    *  @param int i_maxket  Maximum number of different kets in the summation.
    *  @ingroup State_management
    */
    state(int i_level,int i_maxket);
    /**
    *  Creates a state specifying the maximum number of kets
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. Intended for internal use.
    *  Useful to use after post-selection where some of the levels have been
    *  eliminated and the level ordering does not coincide anymore with the one
    *  defined by the circuit.
    *  @param int i_level   Number of levels to describe the state.
    *  @param int i_maxket  Maximum number of different kets in the summation.
    *  @param int i_vis  Vector that translates the state levels to the circuit levels defined by the circuit indices.
    *  @ingroup State_management
    */
    state(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys a state object.
    *  @ingroup State_management
    */
    ~state();
    /**
    *  Copies a state into a new one.
    *  @return Copy of the state.
    *  @ingroup State_management
    */
    state *clone();
    /**
    *   Re-initializes a state to zero terms.
    *  @ingroup State_management
    */
    void clear_state();


    // State manipulation methods
    /** @defgroup State_operations State operations
    *   @ingroup State
    *   List of state possible operations.
    */

    /**
    *   Adds a new term/ket to a state.
    *  @param cmplx i_ampl  Amplitude of the new term.
    *  @param int   *occ    List with the occupation of each level in the new term.
    *  @ingroup State_operations
    */
    int add_term(cmplx i_ampl, int *occ);
    /**
    *   Adds a new term/ket to a state (from a human readable configuration).
    *  @param cmplx      i_ampl  Amplitude of the new term.
    *  @param hterm      term    Matrix that defines the new term. Each column defines the configuration of one level.<br>
    *                            There are four different possible kind of configurations depending on the number of rows
    *                            of the term matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavefunction/packet and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. The wavefunciont/time is assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavefinction are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in order.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore every configuration not mentioned is assumed to be zero.
    *
    *  @param qocircuit *qoc     Circuit to which the term is related (to translate from human to numeric form).
    *  @ingroup State_operations
    */
    int add_term(cmplx i_ampl, hterm term, qocircuit *qoc);
    /**
    *  Direct product of states where occupations are different from zero at different levels "state x rhs".<br>
    *  @param state *rhs  State in the right hand side of the direct product.
    *  @ingroup State_operations
    */
    void dproduct(state *rhs);
    /**
    *  Performs the bra-ket operation <bra|state>.
    *  @return The complex number result of the projection.
    *  @ingroup State_operations
    */
    cmplx braket(state *bra);
    /**
    *  Normalizes the state.
    *  @ingroup State_operations
    */
    void normalize();
    /**
    *  Calculates a post-selected state using a "projector" where it is defined the occupation of the levels to be post-selected.
    *  @param state *prj   Projector with the description of the levels (and their occupations) to be post-selected.
    *  @return Returns the post selected state.
    *  @ingroup State_operations
    */
    state *post_selection(state *prj);
    /**
    *  Remove all the levels corresponding to the specified channels provided that they are zero.
    *  It also removes the levels corresponding with a non zero wave-packet definition if it==0.
    *  @param veci *ch   Channel list where the post-selection is performed.
    *  @param int it    Allowed "times" 0=Only 0/ 1=All.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @return Returns the post selected state.
    *  @ingroup State_operations
    */
    state *prepare_reference(veci ch, int it,qocircuit *qoc);

    //Print methods
    /** @defgroup State_output State output
    *   @ingroup State
    *   Methods to print a state on screen.
    */
    /**
    *  Prints a state. A state is the sum of various kets multiplied by their amplitudes.
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  Kets are printed in the same format that prnt_ket for the case where
    *  qocircuit::flag_prnt =0 .
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup State_output
    */
    void  prnt_state();
    /**
    *  Prints a state. A state is the sum of various kets multiplied by their amplitudes.
    *  If column is set to one the kets and their amplitudes are printed in column.
    *  @code
        |ket1>: A1
        |ket2>: A2
        .
        .
        .
        |ketn>: An
    *  @endcode
    *  Kets are printed in the same format that prnt_ket for the case where
    *  qocircuit::flag_prnt =0 .
    *  @param column Do we want each ket term written on a different line (useful when printing long states) 0=No/1=Yes.
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup State_output
    */
    void  prnt_state(int column);
    /**
    *  Prints a state in human readable form. A state is the sum of various kets multiplied by their amplitudes.
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  Kets are printed in the same format that prnt_ket depending on
    *  qocircuit::flag_prnt of the given circuit.
    *  @param int column Do we want each ket term written on a different line (useful when printing long states) 0=No/1=Yes.
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup State_output
    */
    void  prnt_state(qocircuit *qoc, int column);
    /**
    *  Prints a state in human readable form. A state is the sum of various kets multiplied by their amplitudes. If loss=True loss channels
    *  are printed in blue.
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  Kets are printed in the same format that prnt_ket depending on
    *  qocircuit::flag_prnt of the given circuit.
    *  @param int column Do we want each ket term written on a different line (useful when printing long states) 0=No/1=Yes.
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup State_output
    */
    void  prnt_state(qocircuit *qoc, bool loss,int column);

    // Emitters/Initial methods
    /**
    *  Quantum dot state generator. It creates a state compatible with the photon emission of the quantum dot specified by the next parameters.
    *  The parameter "ch" provides information about the configuration of the emitted photons.
    *  @param mati ch Photon configuration matrix.<br>
    *  For each emitted photon we must define (in column) three numbers. Channel and wave function numbers in case photons are horizontally H or vertically V polarized.
    *  There are as many columns as photons defined. Photons are created in pairs using the QDPairs function therefore "ch" has to be defined with an even number of columns.
    *  Distinguishability is defined in the configuration of the emitter circuit element.
    *  @param double k Fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multiphoton emission.
    *  @param double S S/hbar Fine Structure Splitting (FSS) constant.
    *  @param double tss Characteristic time of spin-scattering tss/tx.
    *  @param double thv Characteristic time of cross-dephasing thv/tx.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @see  void emitter(matc D)
    *  @see  void emitter(photon_mdl *mdl)
    *  @see  void QDPair(int i_ch0,int i_ch1, veci i_t, double dt, double k, double S, double tss, double thv, qocircuit *qoc)
    *  @ingroup Emitter_state
    */
    void QD(mati ch, double k, double S, double tss, double thv, qocircuit *qoc);
    /**
    *  Quantum dot photon pair generator. It creates a photon pair as it would be generated by a quantum dot excited by a laser pulse.
    *  Quantum dots emit Bell states, correlated states and random states with a certain probability. Therefore "QDPair" calls the functions "Bell", "Corr" , "Rand".
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t;  Vector with the wavepacket numbers to be assigned to each photon.
    *  @param double k Fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multiphoton emission.
    *  @param double S S/hbar Fine Structure Splitting (FSS) constant.
    *  @param double tss Characteristic time of spin-scattering tss/tx.
    *  @param double thv Characteristic time of cross-dephasing thv/tx.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @see  matc qocircuit::emitter(matc D)
    *  @see  void Bell(int i_ch0,int i_ch1, veci i_t, double phi,char kind, qocircuit *qoc
    *  @see  void Corr(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc)
    *  @see  void Rand(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc)
    *  @ingroup Emitter_state
    */
    void QDPair(int i_ch0,int i_ch1, veci i_t, double dt, double k, double S, double tss, double thv, qocircuit *qoc);
    /**
    *  It creates a photon pair in a Bell state.
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t;  Vector with the wavepacket numbers to be assigned to each photon.
    *  @param double phi   Phase difference between the first and second ket in the definition of the Bell state.
    *  @code
        | A, B > + exp(-i phi)| C, D >
    *  @endcode
    *  @param char   kind  Which of the four bell states is created. <br>
    *  @code
       '+'=|HH> + |VV>
       '-'=|HH> - |VV>
       'p'=|HV> + |VH>
       'm'=|HV> - |VH>
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Bell(int i_ch0,int i_ch1, veci i_t, double phi,char kind, qocircuit *qoc);
    /**
    *  It creates a photon pair correlated in polarization.
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t;  Vector with the wavepacket numbers to be assigned to each photon.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Corr(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);
    /**
    *  It creates a photon pair with random polarization.
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t;  Vector with the wavepacket numbers to be assigned to each photon.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Rand(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);

protected:
    // Auxiliary methods
    /** @defgroup Aux_operations Auxiliary state operations
    *   @ingroup State
    *   List of auxiliary state operations
    */
    /**
    *  Prints a state in rows.
    *  This method is written in support of @see void qocircuit::prnt_state(qocircuit *qoc, int column, bool loss);
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @param loss If true the last nch/2 channels are printed in blue.
    *  @ingroup Aux_operations
    *  @see void qocircuit::prnt_state(qocircuit *qoc, int column, bool loss);
    */
    void prnt_in_rows(qocircuit *qoc, bool loss);
    /**
    *  Prints a state in columns.
    *  This method is written in support of @see void qocircuit::prnt_state(qocircuit *qoc, int column, bool loss);
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @param loss If true the last nch/2 channels are printed in blue.
    *  @ingroup Aux_operations
    *  @see void qocircuit::prnt_state(qocircuit *qoc, int column, bool loss);
    */
    void prnt_in_cols(qocircuit *qoc, bool loss);
};


/** \class projector
*   \brief Is a variation of of the class state to be used
*          as the definition of the post-selection projectors.
*
*          It introduces levels of negative occupation that are
*          the ones ignored in the post-selection process.
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup State
*/
class projector : public state{
public:
    // Public methods
    // Management methods
    /**
    *  Creates a projector object. The maximum number of kets is set by default.
    *  @param int i_level   Number of levels to describe the projector.
    */
    projector(int i_level);
    /**
    *  Creates a projector specifying the maximum number of kets.
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different kets in the projector.
    */
    projector(int i_level,int i_maxket);
    /**
    *  Creates a projector specifying the maximum number of kets
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. Intended for internal use.
    *  Useful to use after a previous post-selection where some of the levels have been
    *  eliminated and the level ordering does not coincide anymore with the one
    *  defined by the circuit.
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different kets in the summation.
    *  @param int i_vis  Vector that translates the projector levels to the circuit levels defined by the circuit indices.
    */
    projector(int i_level, int i_maxket, int *i_vis);

    // Auxiliary methods
    /**
    *  Auxiliary method to create a projector.
    *  Same that to create a state but non-defined
    *  occupation levels are not assumed to be 0 occupied.
    *  They are initialized with a negative value instead.
    *  This way these levels are later on ignored/not
    *  considered in the post-selection operation.
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different kets in the projector.
    */
    void create_projector(int i_level, int i_maxket);
};


/** @defgroup P_Bin Probability bin
 *  Set of probability bins to contabilize samples.
 */

/** \class p_bin
*   \brief Contains all the information
*          to create and manipulate a set of probability bins defined as levels and occupations
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup P_Bin
*
*/
class p_bin: public ket_list{
public:
    double *p; ///< Probability of each samples state
    int N;     ///< Number of samples considered in the bn list


    // Public methods
    // Management methods
    /** @defgroup Bin_management Set of probability bins management
    *   @ingroup P_Bin
    *   Creation and management of a set of probability bins.
    */

    /**
    *  Creates a set of probability bins. The maximum number of bins is set by default.
    *  @param int i_level   Number of levels to describe the bins.
    *  @ingroup Bin_management
    */
    p_bin(int i_level);
    /**
    *  Creates a set of probability bins specifying the maximum number of bins.
    *  @param int i_level   Number of levels to describe the bins.
    *  @param int i_maxket  Maximum number of different bins in the set.
    *  @ingroup Bin_management
    */
    p_bin(int i_level,int i_maxket);
    /**
    *  Creates a set of probability bins specifying the maximum number of bins
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. Intended for internal use.
    *  Useful to use after a previous post-selection where some of the levels have been
    *  eliminated and the level ordering does not coincide anymore with the one
    *  defined by the circuit.
    *  @param int i_level   Number of levels to describe the bins.
    *  @param int i_maxket  Maximum number of different bins in the set.
    *  @param int i_vis  Vector that translates the p_bin levels to the circuit levels defined by the circuit indices.
    *  @ingroup Bin_management
    */
    p_bin(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys a set of probability bins object.
    *  @ingroup Bin_management
    */
    ~p_bin();
    /**
    *   Creates a copy of this set of probability bins.
    *   @returns    A copy of the present set of bins
    *   @ingroup Bin_management
    */
    p_bin *clone();


    // Bin update methods
    /** @defgroup Bin_update Set of probability bins update
    *   @ingroup P_Bin
    *   List of operations to add statistics to a set of probability bins
    */
    /**
    *  Counts a new sample.
    *  @param int   *occ    Occupation of the sampled state.
    *  @ingroup Bin_update
    */
    int add_count(int *occ);
    /**
    *  Adds the statistics from another bin.
    *  @param p_bin *input  Probability bins with additional statistics.
    *  @ingroup Bin_update
    */
    void add_bin(p_bin *input);
    /**
    *  Adds the statistics of a quantum state.
    *  @param state *input  Quantum state.
    *  @ingroup Bin_update
    */
    void add_state(state *input);

    // Bin manipulation methods
    /** @defgroup Bin_manipulation Set of probability bins operations
    *   @ingroup P_Bin
    *   List of operations to manipulate set of probability effects and to consider classical effects to measurement
    */
    /**
    *  Calculates the effects of the detectors over the outcome contained in this set of probability bins. This effects are described in the circuit.
    *  They are post-selection-conditions, dark counts, detector dead time effects, loss computations and circuit noise.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of physical detectors.
    *  @ingroup Bin_manipulation
    */
    p_bin *calc_measure(qocircuit *qoc);
    /**
    *  Calculates the labels of the bins after a post-selection operation using a "projector". If two resulting labels coincide their counts
    *  are summed if the post-selection is not possible the bin is erased.
    *  @param state *prj   Projector with the description of the levels (and their occupations) to be post-selected.
    *  @return Returns the post selected set of bins.
    *  @ingroup Bin_manipulation
    */
    p_bin *post_selection(state *prj);
    /**
    *  Computes dark counts effects. <br>
    *  <b> Intended for internal use of the library </b>
    *  @param int S The size of the sample for dark counts.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of dark counts.
    *  @ingroup Bin_manipulation
    */
    p_bin *dark_counts(int S, qocircuit* qoc);
    /**
    *  Computes detector dead time effects due other previous measurements. <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param int S The size of the sample for the calculation of detector dead time effects.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of detectors dead time.
    *  @ingroup Bin_manipulation
    */
    p_bin *blink(int S, qocircuit* qoc);
    /**
    *  It removes the extra channels used by the simulator to compute losses. This way hiding which channels are the responsible for the most losses.
    *  This is the usual situation in a real experiment with losses where the fate of the lost photons is unknown. <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param int ndecq Number of detector.
    *  @param mati def  Post-selection conditions if they are also present @see compute_cond .
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with los loss channels properly filtered.
    *  @see compute_cond(int ndec, mati def,qocircuit *qoc);
    *  @ingroup Bin_manipulation
    */
    p_bin *compute_loss(int ndec, mati def,qocircuit *qoc);
    /**
    *  It perform the post-selections compatible with the detector conditions defined in the circuit. <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param int ndecq Number of detector.
    *  @param mati def  Table with the  post-selection condition.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with the proper post-selections applied.
    *  @ingroup Bin_manipulation
    */
    p_bin *compute_cond(int ndec, mati def,qocircuit *qoc);
    /**
    *  Calculates the outcomes removing the frequency degrees of freedom and leaving only times. <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcomes indexed by times instead of packets with the frequency properly integrated out.
    *  @ingroup Bin_manipulation
    */
    p_bin *remove_freq(qocircuit* qoc);
    /**
    *  Calculates the outcomes if we ignore the packet degrees of freedom.<br>
    *  <b> Intended for internal use of the library </b>.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities independent of the packet index.
    *  @ingroup Bin_manipulation
    */
    p_bin *perform_count(qocircuit* qoc);
    /**
    *  It removes the packet of degrees of freedom (except one). <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with the proper post-selections applied.
    *  @ingroup Bin_manipulation
    */
    p_bin *remove_time(qocircuit *qoc);
    /**
    *  Computes gaussian white noise effects in the output. <br>
    *  <b> Intended for internal use of the library </b>.
    *  @param double stdev  Standard deviation of the Gaussian white noise.
    *  @return Returns a list of outcome probabilities with an added Gaussian white noise.
    *  @ingroup Bin_manipulation
    */
    p_bin *white_noise(double stdev);

    //Bin consultation methods
    /** @defgroup Bin_consult Set of probability bins consult statistics
    *   @ingroup P_Bin
    *   List of operations to read the statistics stored in the set of probability bins
    */
    /**
    *  Given a bin index returns the occupation of a bin in a string format.
    *  @param int index  Bin index.
    *  @return Returns a string with the bin name/occupation.
    *  @ingroup Bin_consult
    */
    string tag(int index);
    /**
    *  Given a bin index returns a the probability of that outcome.
    *  @param int index  Bin index.
    *  @return Returns the probability of the referred bin.
    *  @ingroup Bin_consult
    */
    double prob(int index);
    /**
    *  Given a bin definition returns a the probability of that outcome.
    *  @param mati def  Matrix that defines the bin. Each column defines the configuration of one level.<br>
    *                            There are four different possible kind of configurations depending on the number of rows
    *                            of the term matrix.<br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavefunction/packet and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. The wavefunciont/time is assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavefinction are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in order.<br>
    *                            <br>
    *  @return Returns the probability of the referred bin.
    *  @ingroup Bin_consult
    */
    double prob(mati def,qocircuit *qoc);

    // Print bins
    /** @defgroup Bin_output Set of bins output.
    *   @ingroup P_Bin
    *   Methods to print a set of bins on screen.
    */

    /**
    *  Prints a set of bins.
    *  Bin tags are printed in the same format that prnt_ket for the case where
    *  qocircuit::flag_prnt =0 .
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup Bin_output
    */
    void  prnt_bins();
    /**
    *  Prints a set of bins in human readable form.
    *  Bin tags are printed in the same format that prnt_ket depending on
    *  qocircuit::flag_prnt of the given circuit.
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @param double thresh Probabilities below this value are not printed.
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup Bin_output
    */
    void  prnt_bins(qocircuit *qoc, double thresh);
    /**
    *  Prints a set of bins in human readable form.
    *  Bin tags are printed in the same format that prnt_ket depending on
    *  qocircuit::flag_prnt of the given circuit.
    *  @param qocircuit *qoc Circuit to which the state is related (to translate from numeric to human form).
    *  @param double thresh Probabilities below this value are not printed.
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @see qocircuit::flag_prnt
    *  @see void prnt_ket()
    *  @ingroup Bin_output
    */
    void  prnt_bins(qocircuit *qoc, double thresh, bool loss);
};


/** @defgroup PBunch Photon bunch
 *  A bunch of photons to be emitted into a circuit
 */

/** \class ph_bunch
*   \brief Contains all the information
*          to create and manipulate a bunch of photons to be emitted into a circuit
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup PBunch
*
*/
class ph_bunch{
    int    npack;       ///< Number of wave-packets
    matd   pack_list;   ///< Wave-packet definitions
    state *bunch;       ///< Photon state
public:
    // Public methods
    // Management methods
    /** @defgroup PBunch_management Photons bunch management
    *   @ingroup PBunch
    *   Creation and management of photon bunches
    */

    /**
    *  Creates a photon bunch. The maximum number of alternatives is set by default.
    *  @param int i_level   Number of levels to describe the photon states.
    *  @ingroup PBunch_management
    */
    ph_bunch(int i_level);
    /**
    *  Creates a photon bunch specifying the maximum number of alternatives (kets in the state).
    *  @param int i_level   Number of levels to describe the photon states.
    *  @param int i_maxket  Maximum number of alternatives (kets in the state).
    *  @ingroup PBunch_management
    */
    ph_bunch(int i_level,int i_maxket);
    /**
    *  Destroys a photon bunch object.
    *  @ingroup PBunch_management
    */
    ~ph_bunch();
    /**
    *  Clears a photon bunch object.
    *  @ingroup PBunch_management
    */
    void clear();

    // Manipulation methods
    /** @defgroup PBunch_manipulation Photons bunch operations
    *   @ingroup PBunch
    *   Operations with photon bunches
    */
    /**
    *  Adds N photons to a photon bunch.
    *  @param int N Number of photons to be added.
    *  @param int ch Channel where the photons are added.
    *  @param qocircuit *qoc Circuit to which the channels are related.
    *  @ingroup PBunch_manipulation
    */
    void add_photons(int N, int ch, qocircuit *qoc);
    /**
    *  Adds N photons to a photon bunch with a defined wavepacket shape.
    *  @param int N Number of photons to be added.
    *  @param int ch Channel where the photons are added.
    *  @param int P Polarization of the photons.
    *  @param double t Emission time of the photons.
    *  @param double f Emission frequency of the photons.
    *  @param double w Emission width/decay time of the photons depending on the packet shape model.
    *  @param qocircuit *qoc Circuit to which the channels are related.
    *  @ingroup PBunch_manipulation
    */
    void add_photons(int N, int ch, int P, double t, double f, double w, qocircuit *qoc);
    /**
    *  Send photons to a circuit (Emit photons).
    *  @param qocircuit *qoc Circuit to which the photons are emitted.
    *  @ingroup PBunch_manipulation
    */
    void send2circuit(qocircuit *qoc);
    /**
    *  Send photons to a circuit specifying the packet shape model.
    *   @param char ckind  Kind of wave packets 'G'=Gaussian, 'E'=Exponential.
    *   @param int rand    How the phase times have been computed. <br>
    *       0= No compute. <br>
    *       1=Same than flight times. <br>
    *       2=Random times (within a certain distribution). <br>
    *  @param qocircuit *qoc Circuit to which the photons are emitted.
    *  @ingroup PBunch_manipulation
    */
    void send2circuit(char ckind, int rand, qocircuit *qoc);
    /**
    *  Amplitude of the current alternative/state ket. If we have a superposition of states
    *  this is the amplitude of the current one. With the ph_bunch::alternative method we change to
    *  a new superposition.
    *  @param cmplx A Amplitude of the current alternative/ket.
    *  @see void ph_bunch::alternative();
    *  @ingroup PBunch_manipulation
    */
    void weight(cmplx A);
    /**
    *  After alternative is called new photons are added to a new state that exists in superpostion
    *  with the previous ones. Useful to define entangled sets of photons.
    *  Example:
    *  @code
        bunch->add_photons(0, 0, circuit);
        bunch->add_photons(1, 0, circuit);
    *  @endcode
    *  means an state |11> while <br>
    *  @code
        bunch->add_photons(0, 0, circuit);
        bunch->alternative();
        bunch->add_photons(1, 0, circuit);
    *  @endcode
    *  means 1/sqrt(2)(|10>+|11>). A weight/amplitude one is assumed if nothing is specified and
    *  the state is normalized when sent to the circuit.
    *  @ingroup PBunch_manipulation
    */
    void alternative();
    /**
    *  Return the photon state stored in the bunch of photons.
    *  @return The state of the defined photons.
    *  @ingroup PBunch_manipulation
    */
    state *bunch_state();

    // Print methods
    /** @defgroup PBunch_output Bunch state output
    *   @ingroup PBunch
    *   Methods to print the bunch photon state on screen.
    */
    /**
    *  Prints the bunch photons state.
    *  Kets are printed in the same format than state::prnt_state when the print flag is zero.
    *  @param column Is each ket written on a different line (useful when printing long states) 0=No/1=Yes.
    *  @param qocircuit* qoc Circuit to which the bunch of photons are emitted.
    *  @see void qocircuit::set_prnt_flag(int flag);
    *  @see void state::prnt_state();
    *  @ingroup PBunch_output
    */
    void  prnt_state();
    /**
    *  Prints the bunch photons state.
    *  Kets are printed in the same format than state::prnt_state .
    *  @param column Is each ket written on a different line (useful when printing long states) 0=No/1=Yes.
    *  @param qocircuit* qoc Circuit to which the bunch of photons are emitted.
    *  @see void qocircuit::set_prnt_flag(int flag);
    *  @see void state::prnt_state();
    *  @ingroup PBunch_output
    */
    void  prnt_state(qocircuit *qoc, int column);
    /**
    *  Prints the bunch photons state. Loss channels are colored in blue if loss=true.
    *  Kets are printed in the same format than state::prnt_state
    *  @param column Is each ket written on a different line (useful when printing long states) 0=No/1=Yes
    *  @param loss Last nch/2 are printed in blue.
    *  @param qocircuit* qoc Circuit to wich the bunch of photons are emitted.
    *  @see void qocircuit::set_prnt_flag(int flag);
    *  @see void state::prnt_state();
    *  @ingroup PBunch_output
    */
    void  prnt_state(qocircuit *qoc, bool loss,int column);

protected:
    // Auxiliary methods
    /** @defgroup Aux_PBunch Auxiliary photon bunch operations
    *   @ingroup PBunch
    *   List of auxiliary  photon bunch operations
    */

    /**
    *  Auxiliary method to create a photon bunch.
    *  @param int i_level  Number of levels of the bunch state.
    *  @param int i_maxket Maximum number of alternatives/kets.
    *  @ingroup Aux_PBunch
    */

    void create_ph_bunch(int i_level, int i_maxket);
};
