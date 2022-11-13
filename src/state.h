/**************************************************************************
* @file state.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Bosonic state library
*
***************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************
class ket_list{
public:
    // Public methods
    // Management methods
    ket_list(int i_level);                                               //  Creates a ket list.The maximum number of kets is set by default.
    ket_list(int i_level,int i_maxket);                                  //  Creates a ket list specifying the maximum number of kets.
    ket_list(int i_level, int i_maxket, int *i_vis);                     //  Creates a ket list specifying the maximum number of kets and a vector of equivalence between state and circuit levels.
    ~ket_list();                                                         //  Destroys a ket_list
    void clone();                                                        //  Copies a ket list
    void clear_kets();                                                   //  Clear the list (without destroying it).

    // State manipulation methods.
    int add_ket(int *occ);                                               // Adds a new ket to a ket list
    int add_ket(hterm term, qocircuit *qoc);                             // Adds a new ket to a ket list ( human readable form)
    int find_ket(int *occ);                                              // Finds the position of a ket in the list given and occupation description
    int find_ket(mati def,qocircuit *qoc);                               // Finds the position of a ket described in human readable form
    ket_list *remove_time(qocircuit *qoc);                               // Remove levels with wave-packet indexes greater than zero. Needed for some specific operations

    //Print methods
    void  prnt_ket(int iket);                                            // Prints a ket
    void  prnt_ket(int iket, int format, qocircuit *qoc);                // Prints a ket
    void  prnt_ket(int iket, int format, bool loss,qocircuit *qoc);      // Prints a ket

protected:
    void create_ket_list(int i_level, int i_maxket);                     // Create ket list auxiliary function
};


class state{
public:
    // Public methods
    // Management methods
    state(int i_level);                                                  //  Creates a state.The maximum number of kets is set by default.
    state(int i_level,int i_maxket);                                     //  Creates a state specifying the maximum number of kets.
    state(int i_level, int i_maxket, int *i_vis);                        //  Creates a state specifying the maximum number of kets and a vector of equivalence between state and circuit levels.

    ~state();                                                            //  Destroys a state
    state *clone();                                                      //  Copy a state
    void clear();                                                        //  Clears a state

    // State manipulation methods.
    void add_term(cmplx i_ampl, int *occ);                               // Adds a new term to a state
    void add_term(cmplx i_ampl, hterm term, qocircuit *qoc);             // Adds a new term to a state ( human readable form)
    void dproduct(state *rhs);                                           // Direct product of states defined in non coincident channels).
    cmplx braket(state *bra);                                            // Calculates the braket <bra|state>
    void normalize();                                                    // Normalizes the state
    state *post_selection(state *prj);                                   // Post select a state using a "projector"
    state *remove_empty_channels(veci ch, int t,qocircuit *qoc);         // Remove all the levels corresponding to the specified channels provided that they are zero.
    state *convert(veci cnv, qocircuit *qoc);                            // Re-arranges the packet numeration

    //Print methods
    void  prnt_state();                                                  // Prints a state ( various kets and amplitudes )
    void  prnt_state(int column);                                        // Prints a state ( various kets and amplitudes )
    void  prnt_state(int format, int column, qocircuit *qoc);            // Prints a state ( in human readable form )
    void  prnt_state(int format, int column, bool loss, qocircuit *qoc); // Prints a state ( in human readable form )

    // Emitters/Initial states
    void QD(mati ch, double k, double S, double tss, double thv, qocircuit *qoc);                                       // QD state generator model
    void Bell (mati ch,char kind, double phi, qocircuit *qoc);                                                          // Non-ideal bell emitter with a phase e^(-i phi) in the second term. (Path encoding)
    void BellP(mati ch,char kind, double phi, qocircuit *qoc);                                                          // Non-ideal bell emitter with a phase e^(-i phi) in the second term. (Polarization encoding)
    void QDPair(int i_ch0,int i_ch1, veci i_t, double dt, double k, double S, double tss, double thv, qocircuit *qoc);  // Single pair photon emission in a QD.
    void Bell_Path(int i_ch0,int i_ch1, veci i_t, char kind, double phi, qocircuit *qoc);                               // Auxiliary method for non-ideal bell emission. (Path encoding)
    void Bell_Pol(int i_ch0,int i_ch1, veci i_t, char kind, double phi,qocircuit *qoc);                                 // Auxiliary method for non-ideal bell emission. (Polarization encoding)
    void Corr_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);                                                       // Auxiliary method for QD emission. Correlated emitter
    void Rand_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);                                                       // Auxiliary method for QD emission. Random emitter

protected:
    // Auxiliary methods
    void prnt_in_rows(qocircuit *qoc, bool loss);                         // Auxiliary method to print a state in rows
    void prnt_in_cols(qocircuit *qoc, bool loss);                         // Auxiliary method to print a state in columns
};

class projector : public state{
public:
    // Public methods
    // Management methods
    projector(int i_level);                                              //  Creates a projector.The maximum number of kets is set by default.
    projector(int i_level,int i_maxket);                                 //  Creates a projector specifying the maximum number of kets.
    projector(int i_level, int i_maxket, int *i_vis);                    //  Creates a projector specifying the maximum number of kets and a vector of equivalence between state and circuit levels.

    // Auxiliary methods
    void create_projector(int i_level, int i_maxket);                    // Create projector auxiliary function
};

***********************************************************************************/


#include "qocircuit.h"

// Alias
typedef MatrixXi hterm;               ///< Type "Human Term" (hterm). Matrix definition for a state in human readable terms.

// Constant defaults
const int    DEFFORMAT   = 3;         ///< Default print format.
const int    DEFSTATEDIM = 50;        ///< Default maximum number of kets.
const double DEFTHOLDPRNT= 0.0001;    ///< Default amplitude magnitude threshold for printing.


//Group definitions
/** @defgroup Ket_List Ket list
 *  List of kets
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
    thash ketindex;        ///< Hash table of the dynamic dictionary of kets
    int **ket;             ///< Ket definitions. Level occupations of each ket/term
    int *vis;              ///< Correspondence vector. Position to level index.
                           ///< It stores to which level correspond each vector position.
                           ///< After post-selection it keeps track of the original level number.

    // Public methods
    // Management methods
    /** @defgroup Ket_management List of kets management
    *   @ingroup Ket_List
    *   Creation and management of a list of kets.
    */

    /**
    *  Creates a list of kets. The maximum number of kets is set by default.
    *
    *  @param int i_level   Number of levels to describe a ket.
    *  @ingroup Ket_management
    */
    ket_list(int i_level);
    /**
    *  Creates a list of kets specifying the maximum number of kets.
    *
    *  @param int i_level   Number of levels to describe a ket.
    *  @param int i_maxket  Maximum number of different kets in the list.
    *  @ingroup Ket_management
    */
    ket_list(int i_level,int i_maxket);
    /**
    *  Creates a list of kets specifying the maximum number of them given the correspondence vector between state levels and
    *  circuit defined levels. <br>
    *  <b> Intended for internal use. </b>
    *
    *  @param int i_level   Number of levels to describe a ket.
    *  @param int i_maxket  Maximum number of different kets in the list.
    *  @param int i_vis  Vector that translates the ket levels to the circuit levels defined by the circuit indexes.
    *  @ingroup Ket_management
    */
    ket_list(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys this list of kets.
    *
    *  @ingroup Ket_management
    */
    ~ket_list();
    /**
    *   Creates a copy of this list of kets.
    *
    *   @returns    A copy of the present list of kets.
    *   @ingroup Ket_management
    */
    ket_list *clone();
    /**
    *  Empties this list of kets.
    *
    *  @ingroup Ket_management
    */
    void clear_kets();

    // State manipulation methods.
    /** @defgroup Ket_operations List of kets operations
    *   @ingroup Ket_List
    *   List of operations of a list of kets.
    */

    /**
    *   Adds a new ket to the list.
    *
    *  @param int   *occ    Array with the occupation of each level in the new ket.
    *  @ingroup Ket_operations
    */
    int add_ket(int *occ);
    /**
    *   Adds a new ket to the list using a description of the ket.
    *
    *  @param hterm      term    Matrix that defines the new ket. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create these configurations depending on the number of rows
    *                            of the term matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore, levels not configured are assumed to be initialized by zero.
    *
    *  @param qocircuit *qoc     Circuit to which the ket is related.
    *  @ingroup Ket_operations
    */
    int add_ket(hterm term, qocircuit *qoc);
    /**
    *  Finds the position of a ket in the list.
    *
    *  @param int   *occ    List with the occupation of each level of the ket we are searching.
    *  @return Position of the ket in the list.
    *  @ingroup Ket_operations
    */
    int find_ket(int *occ);
    /**
    *  Finds the position of a ket in the list using a description of the ket.
    *
    *  @param hterm      term    Matrix that defines the ket to be found. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create these configurations depending on the number of rows
    *                            of the term matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore, levels not configured are assumed to be zero.
    *
    *  @param qocircuit *qoc     Circuit to which the ket is related.
    *  @return Position of the ket in the list.
    *  @ingroup Ket_operations
    */
    int find_ket(hterm def,qocircuit *qoc);
    /**
    *  Removes levels with a wavepacket index grater than zero.
    *  @param qocircuit *qoc     Circuit to which the ket list is related.
    *  @ingroup Ket_operations
    */
    ket_list *remove_time(qocircuit *qoc);

    //Print methods
    /** @defgroup Ket_output List of kets output
    *   @ingroup Ket_List
    *   Methods to print on screren the content of a ket.
    */

    /**
    *  Prints a ket of the list given its position in the list.
    *
    *  @param int iket           Index of the ket to be printed.
    *  @ingroup Ket_output
    */
    void  prnt_ket(int iket);
    /**
    *  Prints a ket of the list given its position in the list.
    *
    *  @param int iket           Index of ket to be printed.
    *  @param int format         Format of the output.<br>
                        <b>0</b> The ket is printed in numerical format: | (channel, polarization, wavefunction): occupation >. <br>
    *                     Detailed format:
    *                     @code
                            |(ch1{, m1}{, t1}): o1,  ...,(chn{, mn}{, tn}): on >
    *                     @endcode
    *                     Example:
    *                     @code
                            | (0, 0, 0): 0, (0, 0, 1): 0, (0, 1, 0): 2, (0, 1, 1): 0 >
    *                    @endcode
    *                    <br>
    *                   <b>1</b> The polarization is referred by a letter H/V: | (channel, H/V, wavefunction): occupation >. <br>
    *                     Detailed format:
    *                     @code
                            | (ch1{, H/V}{, t1}): o1, ...,(chn{, H/V}{, tn}): on >
    *                     @endcode
    *                     Example:
    *                     @code
                            | (0, H, 0): 0, (0, H, 1): 0, (0, V, 0): 2, (0, V, 1): 0 >
    *                     @endcode
    *                     <br>
    *                  <b>2</b>  Human readable format: | [occupation]H/V(wavefunction)channel >  <br>
    *                     Detailed format:
    *                     @code
                            {{ | {{[o1]}o1>1}{H/V}{(t1)}ch1 > } o1>0},...,{{ | {{[on]}o1>1}{{H/V}{(tn)}chn > } on>0}
    *                     @endcode
    *                     Example:
    *                     @code
                            | [2]V(0)0 >
    *                     @endcode
    *                     <br>
                       <b>3</b>  Straightforward form: <br>
    *                     The occupation values are printed in ascending order. The method reverts to this format
    *                     independently of the configuration if the pointer to the circuit "qocircuit* qoc" is Null<br>
    *                     Example:<br>
    *                     @code
                            | 0, 0, 2, 0 >
    *                     @endcode
    *  <br>
    *  <b>Legend</b>:<br>
    *  chi:Channel i, mi: polarization i, ti: wavefunction i. <br>
    *  {}   Optional values only to be printed if there is more than one possibility. <br>
    *  {{}v>i} Optional values only to be printed if the variable v>i.<br>
    *  @param qocircuit *qoc Circuit to which the ket is related.
    *  @ingroup Ket_output
    */
    void  prnt_ket(int iket, int format, qocircuit *qoc);
    /**
    *  Prints a ket of the list given its position in the list. Loss channels may be printed in different color if explicit losses are considered.
    *
    *  @param int iket           Number of ket to be printed.
    *  @param int format         Format of the output.<br>
                        <b>0</b> The ket is printed in numerical format: | (channel, polarization, wavefunction): occupation >. <br>
    *                     Detailed format:
    *                     @code
                            |(ch1{, m1}{, t1}): o1,  ...,(chn{, mn}{, tn}): on >
    *                     @endcode
    *                     Example:
    *                     @code
                            | (0, 0, 0): 0, (0, 0, 1): 0, (0, 1, 0): 2, (0, 1, 1): 0 >
    *                    @endcode
    *                    <br>
    *                   <b>1</b> The polarization is referred by a letter H/V: | (channel, H/V, wavefunction): occupation >. <br>
    *                     Detailed format:
    *                     @code
                            | (ch1{, H/V}{, t1}): o1, ...,(chn{, H/V}{, tn}): on >
    *                     @endcode
    *                     Example:
    *                     @code
                            | (0, H, 0): 0, (0, H, 1): 0, (0, V, 0): 2, (0, V, 1): 0 >
    *                     @endcode
    *                     <br>
    *                  <b>2</b>  Human readable format: | [occupation]H/V(wavefunction)channel >  <br>
    *                     Detailed format:
    *                     @code
                            {{ | {{[o1]}o1>1}{H/V}{(t1)}ch1 > } o1>0},...,{{ | {{[on]}o1>1}{{H/V}{(tn)}chn > } on>0}
    *                     @endcode
    *                     Example:
    *                     @code
                            | [2]V(0)0 >
    *                     @endcode
    *                     <br>
                       <b>3</b>  Straightforward format: <br>
    *                     The occupation values are printed in ascending order. The method reverts to this format
    *                     independently of the configuration if the pointer to the circuit "qocircuit* qoc" is Null<br>
    *                     Example:<br>
    *                     @code
                            | 0, 0, 2, 0 >
    *                     @endcode
    *  <br>
    *  <b>Legend</b>:<br>
    *  chi:Channel i, mi: polarization i, ti: wavefunction i. <br>
    *  {}   Optional values only to be printed if there is more than one possibility. <br>
    *  {{}v>i} Optional values only to be printed if the variable v>i.<br>
    *  @param bool loss  Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the ket is related.
    *  @ingroup Ket_output
    */
    void  prnt_ket(int iket, int format, bool loss,qocircuit *qoc);

protected:

    /**
    *  Auxiliary method to create a ket list.
    *  Each ket is described as level occupations.
    *
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
    *
    *  @param int i_level   Number of levels to describe the state.
    *  @ingroup State_management
    */
    state(int i_level);
    /**
    *  Creates a state object specifying the maximum number of kets.
    *
    *  @param int i_level   Number of levels to describe the state.
    *  @param int i_maxket  Maximum number of different kets in the summation.
    *  @ingroup State_management
    */
    state(int i_level,int i_maxket);
    /**
    *  Creates a state specifying the maximum number of kets
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. <br>
    *  <b>Intended for internal use.</b>
    *
    *  @param int i_level   Number of levels to describe the state.
    *  @param int i_maxket  Maximum number of different kets in the summation.
    *  @param int i_vis  Vector that translates the state levels to the circuit levels defined by the circuit indexes.
    *  @ingroup State_management
    */
    state(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys a state object.
    *
    *  @ingroup State_management
    */
    ~state();
    /**
    *  Copies a state into a new one.
    *
    *  @return Copy of the state.
    *  @ingroup State_management
    */
    state *clone();
    /**
    *   Re-initializes a state to zero terms.
    *
    *  @ingroup State_management
    */
    void clear();


    // State manipulation methods
    /** @defgroup State_operations State operations
    *   @ingroup State
    *   List of possible operations of a state.
    */

    /**
    *   Adds a new term to the state.
    *
    *  @param cmplx i_ampl  Amplitude of the new term.
    *  @param int   *occ    List with the occupation of each level in the new term.
    *  @ingroup State_operations
    */
    int add_term(cmplx i_ampl, int *occ);
    /**
    *   Adds a new term to the state using a description of the term.
    *
    *  @param hterm      term    Matrix that defines the new tern. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create these configurations depending on the number of rows
    *                            of the term matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    *   Except in the 1-Row case the order is irrelevant. Furthermore, levels not configured are assumed to be initialized by zero.
    *
    *  @param qocircuit *qoc     Circuit to which the term is related.
    *  @return index of the new term in the list defining the state.
    *  @ingroup State_operations
    */
    int add_term(cmplx i_ampl, hterm term, qocircuit *qoc);
    /**
    *  Direct product like operation where state occupations between states are assumed to be different from zero at different levels.
    *  The kets of the resulting state are the result of the sum of the occupations of both the present and input state.
    *
    *  @param state *rhs  State in the right hand side of the direct product like operation.
    *  @return index of the new ket in the ket list defining the state.
    *  @ingroup State_operations
    */
    void dproduct(state *rhs);
    /**
    *  Performs the bra-ket operation <bra|state>.
    *
    *  @return The complex number result of the projection.
    *  @ingroup State_operations
    */
    cmplx braket(state *bra);
    /**
    *  Normalizes the state to one.
    *
    *  @ingroup State_operations
    */
    void normalize();
    /**
    *  Post-selection over states by the condition defined in the the "projector".
    *
    *  @param state *prj   Projector with the description of the levels (and their occupations) to be post-selected.
    *  @return Returns the post-selected state.
    *  @ingroup State_operations
    */
    state *post_selection(state *prj);
    /**
    *  Remove all the levels corresponding to the specified channels provided that they are zero.
    *  Additionally it also removes from the rest of the channels the levels corresponding with
    * a non zero wave-packet index if different packets are not allowed.
    *
    *  @param veci ch   List of channels to be removed.
    *  @param int it    Allowed "times" 0=Only 0/ 1=All.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @return Returns the post-selected state.
    *  @ingroup State_operations
    */
    state *remove_empty_channels(veci ch, int it,qocircuit *qoc);
    /**
    *  Reassign the packet numbers according with the new definition provided by cnv.<br>
    *  <b> Intended for internal use of the library </b>
    *
    *  @param veci cnv  Conversion vector where the index is the old packet number and the contents the new.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @return Returns the re-arranged state.
    *  @ingroup State_operations
    */
    state *convert(veci cnv, qocircuit *qoc);

    //Print methods
    /** @defgroup State_output State output
    *   @ingroup State
    *   Methods to print a state on screen.
    */
    /**
    *  Prints a state.
    *
    *  @ingroup State_output
    */
    void  prnt_state();
    /**
    *  Prints a state.
    *
    *  @param column Are the terms printed in column 0=No/1=Yes.<br>
    *  If column is set to zero the terms are printed in a row,
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  otherwise, if column is set to one the terms are printed in column.
    *  @code
        |ket1>: A1
        |ket2>: A2
        .
        .
        .
        |ketn>: An
    *  @endcode
    *  @ingroup State_output
    */
    void  prnt_state(int column);
    /**
    *  Prints a state in one of the specified the formats.
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket
    *  @param column Are the terms printed in column 0=No/1=Yes.<br>
    *  If column is set to zero the terms are printed in a row,
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  otherwise, if column is set to one the terms are printed in column.
    *  @code
        |ket1>: A1
        |ket2>: A2
        .
        .
        .
        |ketn>: An
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup State_output
    */
    void  prnt_state(int format, int column, qocircuit *qoc);
    /**
    *  Prints a state in one of the specified the formats. If loss=True loss channels are printed in blue.
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket
    *  @param column Are the terms printed in column 0=No/1=Yes.<br>
    *  If column is set to zero the terms are printed in a row,
    *  @code
        A1|ket1>+A2|ket2>+...+An|ketn>
    *  @endcode
    *  otherwise, if column is set to one the terms are printed in column.
    *  @code
        |ket1>: A1
        |ket2>: A2
        .
        .
        .
        |ketn>: An
    *  @endcode
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup State_output
    */
    void  prnt_state(int format, int column, bool loss, qocircuit *qoc);

    // Emitters/Initial methods
    /** @defgroup Emitter_state Initialization methods
    *   @ingroup State
    *   List of methods to create states compatible with an emitter specifications.
    */
    /**
    *  Quantum dot state generator. It creates a state compatible with the photon emission of a quantum dot.

    *  @param mati ch Photon configuration matrix.<br>
    *  For each emitted photon we must define (in column) three numbers.
    *  - 1st Channel.
    *  - 2nd Wavepacket number if the photon is emitted in a state of horizontal polarization.
    *  - 3rd Wavepacket number if the photon is emitted in a state of vertical polarization.
    *  @param double k Fraction of the emitted photon pairs which originate from a XX-X cascade over a background of noise.
    *  @param double S S/hbar Fine Structure Splitting (FSS) constant.
    *  @param double tx  Exciton lifetime
    *  @param double tss Characteristic time of spin-scattering tss.
    *  @param double thv Characteristic time of cross-dephasing thv.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void QD(mati ch, double k, double S, double tx, double tss, double thv, qocircuit *qoc);
    /**
    *  Method to create a photonic Bell state in path encoding.
    *
    *  @param mati ch Photon configuration matrix.<br>
    *  For each emitted photon we must define (in column) three numbers.
    *  - 1st Channel.
    *  - 2nd Wavepacket number if the photon is emitted in a state of horizontal polarization.
    *  - 3rd Wavepacket number if the photon is emitted in a state of vertical polarization.
    *  For each channel we must define (in column) two numbers. Channel and the corresponding wave function number for the photons emitted in that channel. There are as many columns as channels.
    *  @param char   kind  Which of the four bell states is created. <br>
    *  @code
       '+'=|00> + |11>
       '-'=|00> - |11>
       'p'=|01> + |10>
       'm'=|01> - |10>
    *  @endcode
    *  @param double phi   Phase difference between the first and second ket in the definition of the Bell state.
    *  @code
        | A, B > + exp(-i phi)| C, D >
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Bell (mati ch, char kind, double phi, qocircuit *qoc);
    /**
    *  Method to create a photonic Bell state in polarization encoding.
    *
    *  @param mati ch Photon configuration matrix.<br>
    *  For each channel we must define (in column) three numbers. Channel and the corresponding wave function numbers for the photons emitted in that channel. One packet number for H polarized
    *  photons and the other for P polarized photons. There are as many columns as channels.
    *  @param char   kind  Which of the four bell states is created. <br>
    *  @code
       '+'=|HH> + |VV>
       '-'=|HH> - |VV>
       'p'=|HV> + |VH>
       'm'=|HV> - |VH>
    *  @endcode
    *  @param double phi   Phase difference between the first and second ket in the definition of the Bell state.
    *  @code
        | A, B > + exp(-i phi)| C, D >
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void BellP(mati ch, char kind, double phi, qocircuit *qoc);
    /**
    *  Quantum dot photon pair generator. It creates a photon pair as it would be generated by a quantum dot excited by a single laser pulse.
    *
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t  Vector with the wavepacket numbers to be assigned to the photon depending on their polarization.
    *  First, the first channel horizontal and vertical wavepackets and then the second channel in the same order 1H1V2H2V.
    *  @param double k Fraction of the emitted photon pairs which originate from a XX-X cascade over a background of noise.
    *  @param double S S/hbar Fine Structure Splitting (FSS) constant.
    *  @param double tx  Exciton lifetime
    *  @param double tss Characteristic time of spin-scattering tss.
    *  @param double thv Characteristic time of cross-dephasing thv.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void QDPair(int i_ch0,int i_ch1, veci i_t, double dt, double k, double S, double tss, double thv, qocircuit *qoc);
    /**
    *  Auxiliary method to create a photonic Bell state in path encoding.
    *
    *  @param int i_ch0 Channel 0,
    *  @param int i_ch1 Channel 1.
    *  @param veci i_t  Vector with the wavepacket numbers to be assigned to the photon depending on their polarization.
    *  First, the first channel horizontal and vertical wavepackets and then the second channel in the same order 1H1V2H2V.
    *  @param char   kind  Which of the four bell states is created. <br>
    *  @code
       '+'=|00> + |11>
       '-'=|00> - |11>
       'p'=|01> + |10>
       'm'=|01> - |10>
    *  @endcode
    *  @param double phi   Phase difference between the first and second ket in the definition of the Bell state.
    *  @code
        | A, B > + exp(-i phi)| C, D >
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Bell_Path(int i_ch0,int i_ch1, veci i_t, char kind, double phi, qocircuit *qoc);
    /**
    *  Auxiliary method to create a photonic Bell state in polarization encoding.
    *
    *  @param int i_ch0 Channel 0.
    *  @param int i_ch1 Channel 1,
    *  @param veci i_t  Vector with the wavepacket numbers to be assigned to the photon depending on their polarization.
    *  First, the first channel horizontal and vertical wavepackets and then the second channel in the same order 1H1V2H2V.
    *  @param char   kind  Which of the four bell states is created. <br>
    *  @code
       '+'=|HH> + |VV>
       '-'=|HH> - |VV>
       'p'=|HV> + |VH>
       'm'=|HV> - |VH>
    *  @endcode
    *  @param double phi   Phase difference between the first and second ket in the definition of the Bell state.
    *  @code
        | A, B > + exp(-i phi)| C, D >
    *  @endcode
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Bell_Pol(int i_ch0,int i_ch1, veci i_t, char kind, double phi, qocircuit *qoc);
    /**
    *  It creates a photon pair correlated in polarization.
    *
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t  Vector with the wavepacket numbers to be assigned to the photon depending on their polarization.
    *  First, the first channel horizontal and vertical wavepackets and then the second channel in the same order 1H1V2H2V.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Corr_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);
    /**
    *  It creates a photon pair with random polarization.
    *
    *  @param int i_ch0 Channel were the first photon is created.
    *  @param int i_ch1 Channel were the second photon is created.
    *  @param veci i_t  Vector with the wavepacket numbers to be assigned to the photon depending on their polarization.
    *  First, the first channel horizontal and vertical wavepackets and then the second channel in the same order 1H1V2H2V.
    *  @param qocircuit *qoc Circuit to which the state to be created is related.
    *  @ingroup Emitter_state
    */
    void Rand_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc);

protected:
    // Auxiliary methods
    /** @defgroup Aux_operations State auxiliary methods
    *   @ingroup State
    *   List of auxiliary state operations
    */
    /**
    *  Prints a state in rows. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  This method is written in support of prnt_state(int format, int column, bool loss, qocircuit *qoc);
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket
    *  @param loss If true the last nch/2 channels are printed in blue.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @ingroup Aux_operations
    *  @see prnt_state(int format, int column, bool loss, qocircuit *qoc);
    */
    void prnt_in_rows(int format, bool loss, qocircuit *qoc);
    /**
    *  Prints a state in columns. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  This method is written in support of prnt_state(int format, int column, bool loss, qocircuit *qoc);
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket
    *  @param loss If true the last nch/2 channels are printed in blue.
    *  @param qocircuit *qoc Circuit to which the state is related.
    *  @ingroup Aux_operations
    *  @see prnt_state(int format, int column, bool loss, qocircuit *qoc);
    */
    void prnt_in_cols(int format, bool loss, qocircuit *qoc);
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
    *  Creates a projector object. The maximum number of terms is set by default.
    *
    *  @param int i_level   Number of levels to describe the projector.
    */
    projector(int i_level);
    /**
    *  Creates a projector specifying the maximum number of terms.
    *
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different kets in the projector.
    */
    projector(int i_level,int i_maxket);
    /**
    *  Creates a projector specifying the maximum number of terms and the vector of equivalence between state levels and
    *  circuit defined levels. <br>
    *  <b>Intended for internal use</b>.
    *
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different terms in the summation.
    *  @param int i_vis  Vector that translates the projector levels to the circuit levels defined by the circuit indices.
    */
    projector(int i_level, int i_maxket, int *i_vis);

    // Auxiliary methods
    /**
    *  Auxiliary method to create a projector.
    *  Same that to create a state but non-defined
    *  occupation levels are not assumed to be 0 but
    *  they are initialized with a negative value instead.
    *  This way these levels are later on ignored/not
    *  considered in the post-selection operation.
    *
    *  @param int i_level   Number of levels to describe the projector.
    *  @param int i_maxket  Maximum number of different terms in the projector.
    */
    void create_projector(int i_level, int i_maxket);
};
