/**************************************************************************
* @file pbin.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Probability bins library
*
***************************************************************************/

/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************
class p_bin: public ket_list{
public:
    // Public methods
    // Management methods
    p_bin(int i_level);                                                        //  Creates a set of probability bins.The maximum number of outcomes is set by default
    p_bin(int i_level,int i_maxket);                                           //  Creates a set of probability bins specifying the maximum number of outcomes
    p_bin(int i_level, int i_maxket, int *i_vis);                              //  Creates a set of probability bins specifying the maximum number of outcomes and a vector of equivalence between state and circuit levels.
    ~p_bin();                                                                  //  Destroys a set of probability bins
    p_bin *clone();                                                            //  Copies a set of probability bins
    void clear();                                                              //  Clears a set of probability bins

    // Bin basic operations
    double trace();                                                            // Obtain the total probability stored.
    void normalize();                                                          // Normalize the total probability to one.

    // Bin update methods
    int add_count(int *occ);                                                   // Counts a new sample
    void add_bin(p_bin *input);                                                // Adds the statistics from another bin
    void add_state(state *input);                                              // Adds the statistics of a quantum state

    // Bin manipulation methods
    p_bin *calc_measure(qocircuit *qoc);                                       // Calculates a measurement as described by the detectors in the circuit
    p_bin *post_selection(state *prj);                                         // Perform post-selection over the bins given a projector
    p_bin *dark_counts(int S, qocircuit* qoc);                                 // Computes dark counts effects
    p_bin *blink(int S, qocircuit* qoc);                                       // Computes detector dead time effects
    p_bin *compute_loss(qocircuit *qoc);                                       // Remove the additional channels to compute losses
    p_bin *compute_ignored(qocircuit *qoc);                                    // Remove channels flagged to be ignored
    p_bin *remove_channels(veci ch, qocircuit *qoc);                           // Remove channels given by a list
    p_bin *compute_cond(qocircuit *qoc);                                       // Applies the detection conditions from the circuit
    p_bin *post_select_cond(int ndec, mati def,qocircuit *qoc);                // Applies the detection conditions from a list
    p_bin *remove_freq(qocircuit* qoc);                                        // Remove frequencies and leaves only times
    p_bin *perform_count(qocircuit* qoc);                                      // Sum all the contributions to a channel independently of time or frequency
    p_bin *remove_time(qocircuit *qoc);                                        // Rome the all the packet definitions that are not 0
    p_bin *white_noise(double stdev);                                          // Adds Gaussian white noise

    //Bin consultation methods
    string tag(int index);                                                     // Returns the occupation of the indexed bin in string format
    double prob(int index);                                                    // Probability value of the indexed bin
    double prob(mati def,qocircuit *qoc);                                      // Gets the probability of a bin defined by def.
    double prob(mati def,qodev *dev);                                          // Gets the probability of a bin defined by def.

    // Print bins
    void prnt_bins();                                                          // Prints the bin list
    void prnt_bins(double thresh);                                             // Prints the bin list with a threshold limit
    void prnt_bins(int format,double thresh,qocircuit *qoc);                   // Prints the bin list ( in human readable form )
    void prnt_bins(int format, double thresh, qodev *dev);                     // Prints the bin list ( in human readable form )
    void prnt_bins(int format,double thresh,bool loss,qocircuit *qoc);         // Prints the bin list ( in human readable form )
    void prnt_bins( int format, double thresh, bool loss, qodev *dev);         // Prints the bin list ( in human readable form )
    void aux_prnt_bins( int format, double thresh, bool loss, qocircuit *qoc); // Auxiliary method to print the bin list

};

***********************************************************************************/


#include "qodev.h"

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
    double *p; ///< Probability of each sampled state
    int N;     ///< Number of samples considered in the bn list


    // Public methods
    // Management methods
    /** @defgroup Bin_management Set of probability bins management
    *   @ingroup P_Bin
    *   Creation and management of a set of probability bins.
    */

    /**
    *  Creates a set of probability bins. The maximum number of bins is set by default.
    *
    *  @param int i_level   Number of levels.
    *  @ingroup Bin_management
    */
    p_bin(int i_level);
    /**
    *  Creates a set of probability bins specifying the maximum number of bins.
    *
    *  @param int i_level   Number of levels.
    *  @param int i_maxket  Maximum number of different bins in the set.
    *  @ingroup Bin_management
    */
    p_bin(int i_level,int i_maxket);
    /**
    *  Creates a set of probability bins specifying the maximum number of bins
    *  and the vector of equivalence between state levels and
    *  circuit defined levels. <br>
    *  <b> Intended for internal use. </b>
    *
    *  @param int i_level   Number of levels.
    *  @param int i_maxket  Maximum number of different bins in the set.
    *  @param int i_vis  Vector that translates the p_bin levels to the circuit levels defined by the circuit indexes.
    *  @ingroup Bin_management
    */
    p_bin(int i_level, int i_maxket, int *i_vis);
    /**
    *  Destroys a set of probability bins object.
    *
    *  @ingroup Bin_management
    */
    ~p_bin();
    /**
    *   Creates a copy of this set of probability bins.
    *
    *   @returns    A copy of the present set of bins
    *   @ingroup Bin_management
    */
    p_bin *clone();
    /**
    *   Re-initializes a p_bin to zero outcomes.
    *
    *  @ingroup Bin_management
    */
    void clear();

    // Bin basic operations
    /** @defgroup Bin_basic Set of probability basic operations
    *   @ingroup P_Bin
    *   List of basic set of probability bins operations
    */
    /**
    *  Obtains the total probability of the set of bins.
    *  It may be inferior to one due post-selection or encoding procedures.
    *
    *  @ingroup Bin_basic
    */
    double trace();
    /**
    *  Normalizes the total probability of the set of bins to one.
    *
    *  @ingroup Bin_basic
    */
    void normalize();

    // Bin update methods
    /** @defgroup Bin_update Set of probability bins update operations
    *   @ingroup P_Bin
    *   List of operations to add statistics to a set of probability bins
    */
    /**
    *  Counts a new sample.
    *
    *  @param int   *occ    Occupation of the sampled state.
    *  @ingroup Bin_update
    */
    int add_count(int *occ);
    /**
    *  Adds the statistics from another bin.
    *
    *  @param p_bin *input  Probability bins with additional statistics.
    *  @ingroup Bin_update
    */
    void add_bin(p_bin *input);
    /**
    *  Adds the statistics of a quantum state.
    *
    *  @param state *input  Quantum state.
    *  @ingroup Bin_update
    */
    void add_state(state *input);


    // Bin manipulation methods
    /** @defgroup Bin_manipulation Set of probability bins measurement operations
    *   @ingroup P_Bin
    *   List of operations to manipulate a set of probability bins to compute post-selection conditions and effects present in physical detectors.
    */
    /**
    *  Calculates the effect of the detectors defined in a circuit over the outcome contained in this set of probability bins.
    *  This effect are post-selection-conditions, dark counts, detector dead time, loss computations and circuit noise.
    *
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of physical detectors.
    *  @see dark_counts(int S, qocircuit* qoc);
    *  @see blink(int S, qocircuit* qoc);
    *  @see compute_loss(qocircuit *qoc);
    *  @see compute_ignored(qocircuit *qoc);
    *  @see compute_cond(qocircuit *qoc);
    *  @see remove_freq(qocircuit* qoc);
    *  @see remove_time(qocircuit *qoc);
    *  @see perform_count(qocircuit* qoc);
    *  @see white_noise(double stdev);
    *  @ingroup Bin_manipulation
    */
    p_bin *calc_measure(qocircuit *qoc);
    /**
    *  Calculates the labels of the bins after a post-selection operation. If two resulting labels coincide their counts
    *  are summed if the post-selection is not possible the bin is erased.
    *
    *  @param state *prj   Projector with the description of the levels (and their occupations) to be post-selected.
    *  @return Returns the post selected set of bins.
    *  @ingroup Bin_manipulation
    */
    p_bin *post_selection(state *prj);
    /**
    *  Computes dark counts effects. <br>
    *  <b> Intended for internal use of the library </b>
    *
    *  @param int S The size of the sample for dark counts.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of dark counts.
    *  @ingroup Bin_manipulation
    */
    p_bin *dark_counts(int S, qocircuit* qoc);
    /**
    *  Computes detector dead time effects due other previous measurements. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param int S The size of the sample for the calculation of detector dead time effects.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities considering the effects of detectors dead time.
    *  @ingroup Bin_manipulation
    */
    p_bin *blink(int S, qocircuit* qoc);
    /**
    *  It removes the extra channels used by the simulator to compute losses explicitly. This way it is hidden which channels are the responsible for the losses.
    *  This is the usual situation in a real experiment where the fate of the lost photons is unknown. The definition of the loss channels
    *  is within the circuit object.<br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit to which the channels are referred.
    *  @return Returns a list of outcome probabilities with the loss channels properly filtered.
    *  @see remove_channels(veci ch, qocircuit *qoc);
    *  @ingroup Bin_manipulation
    */
    p_bin *compute_loss(qocircuit *qoc);
    /**
    *  It removes channels flagged to be ignored by the user. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit to which the channels are referred.
    *  @return Returns a list of outcome probabilities with the ignored channels properly filtered.
    *  @see remove_channels(veci ch, qocircuit *qoc);
    *  @ingroup Bin_manipulation
    */
    p_bin *compute_ignored(qocircuit *qoc);
    /**
    *  General method to remove channels from a list. It perform post-selection by all the possible outcomes of those channels (the channels are traced out from the result). <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param veci ch List of channels to be traced out.
    *  @param qocircuit *qoc  Circuit to which the channels are referred.
    *  @return Returns a list of outcome probabilities with the indicated channels properly filtered.
    *  @see post_select_cond(int ndec, mati def,qocircuit *qoc);
    *  @ingroup Bin_manipulation
    */
    p_bin *remove_channels(veci ch, qocircuit *qoc);

    /**
    *  It perform post-selections compatible with the detector conditions defined in the circuit. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with the proper post-selections applied.
    *  @see post_select_cond(int ndec, mati def,qocircuit *qoc);
    *  @ingroup Bin_manipulation
    */
    p_bin *compute_cond(qocircuit *qoc);
    /**
    *  General method to perform post-selections compatible with a list of detector conditions. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param int ndec Number of detector.
    *  @param mati def  Table with the  post-selection condition.
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with the proper post-selections applied.
    *  @ingroup Bin_manipulation
    */
    p_bin *post_select_cond(int ndec, mati def,qocircuit *qoc);
    /**
    *  Calculates the outcomes removing the frequency degrees of freedom and leaving only times. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcomes indexed by times instead of packets with the frequency properly integrated out.
    *  @ingroup Bin_manipulation
    */
    p_bin *remove_freq(qocircuit* qoc);
    /**
    *  Calculates the outcomes if we ignore the packet degrees of freedom.<br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities independent of the packet index.
    *  @ingroup Bin_manipulation
    */
    p_bin *perform_count(qocircuit* qoc);
    /**
    *  It removes the packet of degrees of freedom (except one). <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param qocircuit *qoc  Circuit where the detectors are defined.
    *  @return Returns a list of outcome probabilities with the proper post-selections applied.
    *  @ingroup Bin_manipulation
    */
    p_bin *remove_time(qocircuit *qoc);
    /**
    *  Computes Gaussian white noise effects in the output. <br>
    *  <b> Intended for internal use of the library </b>.
    *
    *  @param double stdev  Standard deviation of the Gaussian white noise.
    *  @return Returns a list of outcome probabilities with an added Gaussian white noise.
    *  @ingroup Bin_manipulation
    */
    p_bin *white_noise(double stdev);

    //Bin consultation methods
    /** @defgroup Bin_consult Set of probability bins consult statistics
    *   @ingroup P_Bin
    *   List of operations to read the statistics stored in a set of probability bins
    */
    /**
    *  Returns the occupation of the bin referred by the index in a string format
    *
    *  @param int index  Bin index.
    *  @return Returns a string with the bin occupation.
    *  @ingroup Bin_consult
    */
    string tag(int index);
    /**
    *  Returns the probability of the outcome stored in the bin referred by the index.
    *
    *  @param int index  Bin index.
    *  @return Returns the probability of the referred bin.
    *  @ingroup Bin_consult
    */
    double prob(int index);
    /**
    *  Returns the probability of the outcome stored in the bin described by the provided definition. (Circuit version)
    *
    *  @param mati def  Matrix that defines the bin. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create a definition depending on the number of rows
    *                            of the matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    *  @param qocircuit *qoc    Circuit to which the outcomes are referred.
    *  @return Returns the probability of the referred bin.
    *  @ingroup Bin_consult
    */
    double prob(mati def,qocircuit *qoc);
    /**
    *  Returns the probability of the outcome stored in the bin described by the provided definition. (Device version)
    *
    *  @param mati def  Matrix that defines the bin. Each column defines the configuration of one level.<br>
    *                            There are four different ways to create a definition depending on the number of rows
    *                            of the matrix: <br>
    *                            <br>
    *                            4-Row: Channels, polarization, wavepacket and occupation in this order.<br>
    *                            3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases.<br>
    *                            2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases.<br>
    *                            1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering.<br>
    *                            <br>
    *  @param qodev *dev    Device to which the outcomes are referred.
    *  @return Returns the probability of the referred bin.
    *  @ingroup Bin_consult
    */
    double prob(mati def,qodev *dev);

    // Print bins
    /** @defgroup Bin_output Set of probability bins output.
    *   @ingroup P_Bin
    *   Methods to print a set of probability bins on screen.
    */

    /**
    *  Prints a set of probability bins.<br>
    *  Bin tags are printed in the same format that prnt_ket for the case where format = 3
    *
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void  prnt_bins();
    /**
    *  Prints a set of probability bins.<br>
    *  Bin tags are printed in the same format that prnt_ket for the case where format = 3.
    *
    *  @param double thresh Probabilities below this value are not printed.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void  prnt_bins(double thresh);
    /**
    *  Prints a set of probability bins. (Circuit version).<br>
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Probabilities below this value are not printed.
    *  @param qocircuit *qoc Circuit to which the set of probability bins is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void prnt_bins(int format, double thresh, qocircuit *qoc);
    /**
    *  Prints a set of probability bins. (Device version).<br>
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Probabilities below this value are not printed.
    *  @param qodev *dev Device to which the set of probability bins is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void prnt_bins(int format, double thresh, qodev *dev);
    /**
    *  Prints a set of probability bins. Loss channels may be colored. (Circuit version).<br>
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Probabilities below this value are not printed.
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the set of probability bins is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void prnt_bins( int format, double thresh, bool loss, qocircuit *qoc);
    /**
    *  Prints a set of probability bins. Loss channels may be colored. (Device version).<br>
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Probabilities below this value are not printed.
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @param qodev *dev Device to which the set of probability bins is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void prnt_bins( int format, double thresh, bool loss, qodev *dev);
    /**
    *  Auxiliary method to print set of bins.
    *
    *  @param format Specifies the format in which the ket is printed. The formats are the same than in prnt_ket.
    *  @param double thresh Probabilities below this value are not printed.
    *  @param bool loss Print loss channels in different color? False=No/True=Yes.
    *  @param qocircuit *qoc Circuit to which the set of probability bins is related.
    *  @see prnt_ket(int iket, int format, bool loss,qocircuit *qoc);
    *  @ingroup Bin_output
    */
    void aux_prnt_bins( int format, double thresh, bool loss, qocircuit *qoc);

};
