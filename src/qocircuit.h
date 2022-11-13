/***********************************************************************************
* @file qocircuit.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Optical circuit library
* @brief A circuit is defined from its optical elements.
*        In this library are implemented the local definitions of these elements
*        and is managed the way they are interconnected.
*
***********************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************
class qocircuit{
public:
    // Public functions
    // Management functionss
    qocircuit(int i_nch);                                                    // Create circuit
    qocircuit(int i_nch, int i_nm, int i_ns);                                // Create circuit
    qocircuit(int i_nch, int i_nm, int i_ns, int clock, char i_ckind);       // Create circuit
    qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);        // Create circuit
    ~qocircuit();                                                            // Destroy circuit
    qocircuit *clone();                                                      // Copy a circuit
    void reset();                                                            // Reset circuit
    int concatenate(qocircuit *qoc);                                         // Concatenates two circuits
    int num_levels();                                                        // Returns the number of levels of this circuit


    // Circuit elements
    //      Basic elements
    void random_circuit();                                                   // Adds a random circuit
    void NSX(int i_ch1, int i_ch2, int i_ch3);                               // Adds a built-in NSX circuit
    void beamsplitter(int i_ch1, int i_ch2, double theta, double phi);       // Adds a ideal beamsplitter to the circuit
    void dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r);                 // Adds a physical dieletric beamsplitter
    void MMI2(int i_ch1, int ch2);                                           // Adds a 2x2 ideal MMI
    void rewire(int i_ch1,int i_ch2);                                        // Adds a swap gate between two channels
    void phase_shifter(int i_ch, double phi);                                // Adds a phase shifter to the circuit
    void dispersion(int i_ch, double dt);                                    // Adds a phase shift that depends on the optical path (in time units)
    void dispersion(int i_ch, int P, double dt);                             // Adds a phase shift that depends on the optical path provided polarization P is fulfilled.
    void phase_shifter(int i_ch, cmplx t);                                   // Adds a General phase shifter with losses
    void loss(int i_ch, double l);                                           // Adds a lossy medium
    void custom_gate(mati iodef, matc mtx);                                  // Adds a custom gate

    //      Polarization elements
    void rotator(int i_ch, double d_theta, double d_phi);                    // Adds a rotator to the circuit.
    void polbeamsplitter(int i_ch1, int i_ch2, int P);                       // Adds a polarizing beamsplitter
    void waveplate(int i_ch, double d_alpha, double d_gamma;                 // Adds a waveplate
    void half(int i_ch, double alpha);                                       // Adds a half-waveplate
    void quarter(int i_ch, double alpha);                                    // Adds a quarter-waveplate


    //      Detection elements
    void ignore(int i_ch);                                                   // Flags a channel to be ignored
    void detector(int i_ch);                                                 // Adds a detector
    void detector(int i_ch, int cond);                                       // Adds a conditional detection
    int detector(int i_ch, int cond, double eff, double blnk, double gamma); // Adds a general physical detector
    int remdec();                                                            // Returns the remaining number of not defined detectors
    void noise(double stdev2);                                               // Adds noise to the output


    //      Emitter and distinguishability model.
    void def_packet(int n, double t, double f, double w);                    // Creates a new packet definition
    cfg_rnd def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand); // Creates a new packet definition
    double emitted_vis(int i,int j);                                         // Probability of two wave packets given by def_packet to overlap
    void emitter ();                                                         // Adds a emitter using the packet definition given by def_packet
    veci emitter (int npack, matd packets, char ckind, int rand);            // Adds a emitter with packets defined in a matrix
    int  emitter(photon_mdl *mdl);                                           // Adds a emitter using a photon model
    void emitter(matc D);                                                    // Adds a emitter using a coupling matrix
    int  delay(int ch, double dt);                                           // Adds a delay using the packet definition given by def_packet
    photon_mdl *delay(int ch, double dt, photon_mdl *mdl);                   // Adds a delay using a given photon model
    void delay(int i_ch, matc T, int nT);                                    // Adds a using an extended coupling matrix


    // Print functions
    void prnt(int format);                                                   // Print circuit matrix


    // Auxiliary functions
    void create_circuit(int i_nch, int i_nm, int i_ns, bool loss);           // Auxiliary function to create circuits
    void compute_losses();                                                   // Auxiliary function to compute the states with less photons at the input than the output.
};


class photon_mdl{
public:
    // Public functions
    // Management functions
    photon_mdl();                                                           // Creates an empty photon model.
    photon_mdl(mati pack_def, matd times, matd freq, char ckind);           // Creates a photon model from parameters
    ~photon_mdl();                                                          // Destroys photon model.
    photon_mdl *clone();                                                    // Creates a copy of a photon model.
    void update(double dt);                                                 // Updates the photon model to include definitions of packets displaced a quantity dt.

    // Utility functions
    matd create_packet_mtx();                                               // Creates a packet matrix model from the packet definitions/specifications
    double visibility(int i,int j);                                         // Probability of two wave packets to overlap

    // Print functions
    void prnt();                                                            //  Prints the information related with a photon model.
    void prnt_times();                                                      //  Prints the the table of times of the photon model.
    void prnt_freqs();                                                      //  Prints the the table of frequencies of the photon model.
    void prnt_packets();                                                    //  Prints the the table of packet index of the photon model.
};

// Auxiliary methods
mati create_packet_idx(mati pack_def);                                      // Converts an unordered packet definition into an ordered one
***********************************************************************************/


#include "util.h"

// Alphabet for a simplified printing.
// Note that the alphabet is designed for binary choices.
// Larger spaces are possible but in this case the user hast to resort
// to a numerical output (also available) or extend the alphabet for
// better human reading.
const char pl[2]={'H','V'};  ///< Alphabet to print polarizations.

// Constant.
const int H=0;               ///< Numerical value for horizontal polarization
const int V=1;               ///< Numerical value for vertical polarization


// Phase configuration in time.
/**
*  \struct cfg_rnd
*  \brief  Phase configuration given in time units.
*/
struct cfg_rnd{
    double dt1;              ///< Phase time 1. Meaning depends on the emitter configuration. See def_packet.
    double dt2;              ///< Phase time 2. Meaning depends on the emitter configuration. See def_packet.
};

//Group definitions
/** @defgroup Photonmdl Packet model
 *  Photon wave packet model class and methods.
 */

/** @defgroup Circuit
 *  Circuit class and methods
 */

/** \class photon_mdl Photon packet model
*   \brief  Contains all the information and methods to create and manipulate
*    a set of photon wave packets
*
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup Photonmdl
*/
class photon_mdl{
public:
    int kind;           ///< Kind of the wave packet 0='Gaussian', 1='Exponential'.
    mati pack_def;      ///< Ordered definition of a photon packet. Each column defines a packet where the first row is an index to the times matrix and the second an index to the freq matrix.
    matd times;         ///< Time definitions. Each column contains one definition. The first row are the relative times to the leading packet while the second the phase times ( The difference is convenient for calculation purposes).
    matd freq;          ///< Frequency definitions. Each column contains one definition. The first row are frequencies while the second are widths (Gaussian) or characteristic decay times (Exponentials).


    // Public functions
    // Management functions
    /** @defgroup Photonmdl_management Packet model management
    *   @ingroup Photonmdl
    *   Creation and management of the photon packet model object.
    */
    /**
    *  Creates an empty packet model.
    *
    *   @ingroup Photonmdl_management
    */
    photon_mdl();
    /**
    *  Creates a packet model given a set of parameters. These are, a packet definition and lists of available times and frequencies.
    *
    *  @param mati pack_def  Packet definition. It is an integer matrix of three rows where each column contains,
    *   a wave packet. <br>
    *                   1st row = an index that labels the packet uniquely.<br>
    *                   2nd row = an index that refers to the central time of the packet as specified in the matrix times.<br>
    *                   3rd row = an index that refers to the central frequency and width of the packet as specified in the matrix freq.<br>
    *   @param  matfd times Two row matrix that contains the discrete sets of times in which a packet may exist in the current simulation
    *                   1st row = is the packet central time.<br>
    *                   2nd row = is the packet initial phase time. If the phase time equals the central time the initial phase is zero.<br>
    *   @param  matd freq Two row matrix that contains the discrete set of frequency configurations that a packet may have in the current simulation.<br>
    *                   1st row = is the packet central frequency.<br>
    *                   2nd row = is the packet width or characteristic time.<br>
    *   @param char ckind  Kind of the wave packet 'G'=Gaussian, 'E'=Exponential.
    *   @ingroup Photonmdl_management
    */
    photon_mdl(mati pack_def, matd times, matd freq, char ckind);
    /**
    *  Destroys the packet model.
    *
    *   @ingroup Photonmdl_management
    */
    ~photon_mdl();
    /**
    *   Creates a copy of this packet model.
    *
    *   @returns    A copy of the present photon model.
    *   @ingroup Photonmdl_management
    */
    photon_mdl *clone();
    /**
    *  Updates the packet model to include definitions of packets displaced by quantity dt with respect the ones already defined by the model.
    *  The update is performed doubling the number of possible wave packets. The new set of wavepackets are displaced in time a quantity dt with respect the originals.
    *
    *  @param double dt Quantity of time the new set of packets are delayed with respect the old ones.
    *  @return An updated photon packet model.
    *  @ingroup Photonmdl_management
    */
    void update(double dt);

    // Utility functions
    /** @defgroup Photonmdl_utility Packet model utilities
    *   @ingroup Photonmdl
    *   Utilities to process the information of a packet model.
    */
    /**
    *  Coverts the packet model into a matrix description useful to calculate overlapping coefficients between packets.
    *  Each column is indexed by the packet number and the rows contain real numbers storing time and frequency values. <br>
    *  <b>This routine is intended for internal use.</b>
    *
    *  @return A photon packet matrix that contains a summary of the packet model information in a convenient order to perform calculations of overlapping matrix.
    *  @ingroup Photonmdl_utility
    */
    matd create_packet_mtx();
    /**
    *  Probability of two wave packets to overlap.
    *
    *  @param int i  Packet i number.
    *  @param int j  Packet j number.
    *  @return Probability of overlapping.
    *  @ingroup Photonmdl_utility
    */
    double visibility(int i,int j);
    /**
    *  Returns the packet definition in a three row format.
    *  This is:
    *       - 1st-Row: Packet number
    *       - 2nd-Row: Index to the time matrix.
    *       - 3rd-Row: Index to the frequency matrix.
    *
    *  @return Three row packet definition.
    *  @ingroup Photonmdl_utility
    */
    mati return_packet_def();

// Print functions
    /** @defgroup Photonmdl_print Packet model print functions
    *   @ingroup Photonmdl
    *   Functions to print on screen and debug the information in a packet model.
    *
    */
    /**
    * Prints on screen the information stored in a packet model.
    *
    *  @ingroup Photonmdl_print
    */
    void prnt();
    /**
    * Prints on screen the table of times of the packet model.
    *
    *  @ingroup Photonmdl_print
    */
    void prnt_times();
    /**
    * Prints on screen the table of frequencies of the packet model.
    *
    *  @ingroup Photonmdl_print
    */
    void prnt_freqs();
    /**
    * Prints on screen the table of packet indexes of the packet model.
    *
    *  @ingroup Photonmdl_print
    */
    void prnt_packets();
};

// Auxiliary methods
/**
*  Converts an unordered packet definition into an ordered one where columns are indexed by the packet number. <br>
*  <b> This routine is intended for internal use of the library. </b>
*
*  @param mati pack_def  A three row matrix packet definition. Where the first one is the packet index, the second one the time index and the third one the frequency index.
*  @return A photon definition where the columns are indexed by the packet index.
*  @ingroup Wavepackets
*/
mati create_packet_idx(mati pack_def);



/** \class qocircuit
*   \brief  Contains all the information to create and
*           manipulate a quantum optical circuit
*
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup Circuit
*/
class qocircuit{
public:
    //Type definitions
    /**
    *  \struct level
    *  \brief  Definition of a level
    */
    struct level{
        int ch;            ///< Channel number.
        int m;             ///< Mode (Polarization).
        int s;             ///< Wavepacket index.
    };

    // Public variables
    int    nlevel;          ///< Number of levels.
    int    nch;             ///< Number of channels.
    int    nm;              ///< Number of modes.
    int    ns;              ///< Number of different wavepackets possible.

    // Dictionary
    level *idx;             ///< Level index: index->level.
    int ***i_idx;           ///< Inverse level index: level->index
    matc   circmtx;         ///< Matrix of the circuit as a function of the level number.

    // Emitter model
    char   ckind;           ///< Kind of emitter 'G'=Gaussian/'E'=Exponential
    int    emiss;           ///< Has been defined an emitter model 0=No/1=Yes
    int    npack;           ///< Number of packets
    matd   pack_list;       ///< Packet list
    matc   init_dmat;       ///< Initial matrix due to emitter
    photon_mdl *emitted;    ///< Packet model of the emitter
    double confidence;      ///< Confidence in the Cholesky approximation.

    // Detector model
    int    losses;          ///< Do we consider a loss model 0=No/1=Yes
    int    ncond;           ///< Post selection condition length
    int    ndetc;           ///< Number of detectors
    int    nignored;        ///< Number of channels ignored (which turn off detectors).
    int    timed;           ///< There is a clock related with the detectors 0=No/1=Yes
    mati   det_def;         ///< Post-selection condition definition
    matd   det_par;         ///< Detector physical parameters
    veci   ch_ignored;      ///< Channels with no detectors or ignored channels.
    int    R;               ///< Number of iterations to calculate detector dead time and dark-counts
    double dev;             ///< Detector standard deviation squared of the Gaussian noise



    // Public functions
    // Management functions
    /** @defgroup Circuit_management Circuit management
    *   @ingroup Circuit
    *   Creation and management of the quantum optical circuit object.
    */

    /**
    *  Creates an optical circuit object.
    *
    *  @param int i_nch  Number of channels.
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch);
    /**
    *  Creates an optical circuit object.
    *
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of packets.
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns);
    /**
    *  Creates an optical circuit object.
    *
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of packets.
    *  @param int clock  Detector behavior configuration.<br>
    *           - 0: No clock. The detectors behave as counters. <br>
    *           - 1: Clock. The detectors are able distinguish arrival times but they are blind to frequency. <br>
    *           - 2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. <br>
    *           - 3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. <br>
    *                 <b>Warning!</b> This is an advanced configuration. It is not recommended unless the user has knowledge on how SOQCS internally works. <br>
    *  @param char ckind Packet model 'G': Gaussian/ 'E': Exponential
    *
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns, int clock, char i_ckind);
    /**
    *  Creates an optical circuit object.
    *
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of packets.
    *  @param int clock  Detector behavior configuration.<br>
    *           - 0: No clock. The detectors behave as counters. <br>
    *           - 1: Clock. The detectors are able distinguish arrival times but they are blind to frequency. <br>
    *           - 2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. <br>
    *           - 3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. <br>
    *                 <b>Warning!</b> This is an advanced configuration. It is not recommended unless the user has knowledge on how SOQCS internally works. <br>
    *
    *
    *  @param int i_R    Number of iterations in the calculation of detector dead time and dark counts effects.
    *  @param bool loss  Are the losses calculated explicitly? True/False.
    *  @param char ckind Packet model 'G': Gaussian/ 'E': Exponential
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    /**
    *  Destroys an optical circuit.
    *
    *  @ingroup Circuit_management
    */
    ~qocircuit();
    /**
    *  Copies a circuit object into a new one.
    *
    *  @return Copy of the ciruit.
    *  @ingroup Circuit_management
    */
    qocircuit *clone();
    /**
    *  Resets circuit.
    *
    *  It resets the circuit back to a blank state except for the definitions in the number of degrees of freedom.
    *  @ingroup Circuit_management
    */
    void reset();
    /**
    *  Concatenates two circuits (with some limitations).<br>
    *  All the packets have to be defined in the first circuit and all the detector in the last.
    *  Circuit with delays can not be concatenated unless the delay is in the first one.
    *
    *  @param qocircuit* qoc Circuit to be concatenated.
    *  @return 0 if success -1 if an error happened.
    *  @ingroup Circuit_management
    */
    int concatenate(qocircuit *qoc);
    /**
    *  Returns the total number of levels of the circuit.
    *
    *  @return Number of levels of the circuit.
    *  @ingroup Circuit_management
    */
    int num_levels();



    // Circuit elements
    /** @defgroup Circuit_elements Circuit elements
    *   @ingroup Circuit
    *   List of optical elements that can be added to a quantum optical circuit.
    */
    //     Basic elements
    /** @defgroup Circuit_basic Basic circuit elements
    *   @ingroup Circuit_elements
    *   List of basic circuit elements
    */

    /**
    *  Adds a circuit defined by a random unitary matrix.
    *
    *  @ingroup Circuit_basic
    */
    void random_circuit();
    /**
    *  Adds a NSX circuit element. Post-selection still has to be carried out to obtain the proper functionality.
    *
    *  @param int i_ch1     NSX input channel 1.
    *  @param int i_ch2     NSX input channel 2.
    *  @param int i_ch3     NSX input channel 3.
    *  @ingroup Circuit_basic
    */
    void NSX(int i_ch1, int i_ch2, int i_ch3);
    /**
    *  Adds a beamsplitter to the circuit attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1     Beamsplitter input channel 1.
    *  @param int i_ch2     Beamsplitter input channel 2.
    *  @param double theta  Angle theta in degrees.
    *  @param double phi    Angle phi in degrees.
    *  @ingroup Circuit_basic
    */
    void beamsplitter(int i_ch1, int i_ch2, double theta, double phi);
    /**
    *  Adds a dieletric film attached to channels i_ch1 and i_ch2
    *  It may also work as a dieletric beamsplitter.
    *
    *  @param int i_ch1     Dielectric input channel 1.
    *  @param int i_ch2     Dielectric input channel 2.
    *  @param cmplx t       Transmission amplitude of probability.
    *  @param cmplx r       Reflection amplitude of probability.
    *  @ingroup Circuit_basic
    */
    void dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r);
    /**
    *  Adds a 2x2 MMI device to the circuit attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1     MMI input channel 1.
    *  @param int i_ch2     MMI input channel 2.
    *  @ingroup Circuit_basic
    */
    void MMI2(int i_ch1, int i_ch2);
    /**
    *  Adds a swap gate between two channels.
    *
    *  @param int i_ch1     Channel 1.
    *  @param int i_ch2     Channel 2.
    *  @ingroup Circuit_basic
    */
    void rewire(int i_ch1,int i_ch2);
    /**
    *  Adds a phase shifter to the circuit in channel i_ch.
    *
    *  @param int i_ch      Phase shifter input channel.
    *  @param double phi    Angle phi in degrees.
    *  @ingroup Circuit_basic
    */
    void phase_shifter(int i_ch, double phi);
    /**
    *  Adds a phase shift to the photons in channel i_ch depending on
    *  their frequency and optical path dt.
    *
    *  @param int i_ch      Phase shifter input channel.
    *  @param double dt     Optical path given in time units.
    *  @ingroup Circuit_basic
    */
    void dispersion(int i_ch, double dt);
    /**
    *  Adds a phase shift to the photons in channel i_ch depending on
    *  their frequency and optical path dt provided the photons match the given polarization P.
    *
    *  @param int i_ch      Phase shifter input channel.
    *  @param int P         Polarization required to the photons to apply the phase shift.
    *  @param double dt     Optical path given in time units.
    *  @ingroup Circuit_basic
    */
    void dispersion(int i_ch, int P, double dt);
    /**
    *  Adds a general phase shifter with losses to the circuit in channel i_ch.
    *
    *  @param int i_ch      Phase shifter input channel.
    *  @param cmplx t       Transmission amplitude.
    *  @ingroup Circuit_basic
    */
    void phase_shifter(int i_ch, cmplx t);
    /**
    *  Adds a lossy medium with loss probability l to the circuit in channel i_ch.
    *
    *  @param int i_ch      Loss medium input channel.
    *  @param double l      Loss probability.
    *  @ingroup Circuit_basic
    */
    void loss(int i_ch, double l);
    /**
    *  Adds a custom gate to the circuit.
    *
    *  @param mati iodef  List of channels and polarizations that define (by column) the input to the gate.
    *  @param matc U      Custom gate nxn matrix definitions.
    *  @ingroup Circuit_basic
    */
    void custom_gate(mati iodef, matc U);

    //     Polarization elements elements
    /** @defgroup Circuit_polar Circuit polarization elements
    *   @ingroup Circuit_elements
    *   List of circuit elements that act on photon polarization
    */
    /**
    *  Adds a polarization rotation device to the circuit attached to the channel i_ch.
    *
    *  @param int i_ch      Rotator input channel.
    *  @param double theta  Angle theta in degrees.
    *  @param double phi    Angle phi in degrees.
    *  @ingroup Circuit_polar
    */
    void rotator(int i_ch, double d_theta, double d_phi);
    /**
    *  Adds a polarized beamsplitter attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1     Beamsplitter input channel 1.
    *  @param int i_ch2     Beamsplitter input channel 2.
    *  @param int P         Polarization to which the beamsplitter is sensitive.
    *  @ingroup Circuit_polar
    */
    void polbeamsplitter(int i_ch1, int i_ch2, int pol);
    /**
    *  Adds a general waveplate attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @param double d_gamma Wavenumber delay of one component with respect the other in degrees.
    *  @ingroup Circuit_polar
    */
    void waveplate(int i_ch, double d_alpha, double d_gamma);
    /**
    *  Adds a half waveplate attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @ingroup Circuit_polar
    */
    void half(int i_ch, double alpha);
    /**
    *  Adds a quarter waveplate attached to channels i_ch1 and i_ch2.
    *
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @ingroup Circuit_polar
    */
    void quarter(int i_ch, double alpha);

    //      Emitter and distinguishability model.
    /** @defgroup Circuit_distin Emission and distinguishability model
    *   @ingroup Circuit_elements
    *   List of advanced circuit elements to configure the photon emission model and the distinguishability between photons in the circuit.
    */

    /**
    *  Adds to the circuit a new definition of a wave packet. <br>
    *  <b> Warning! </b> Note that when the detectors are not configured as counters the packet numeration may change on emission.
    *
    *  @param int n Suggested wavepacket number.
    *  @param double t Wavepacket characteristic emission time.
    *  @param double f Wavepacket characteristic frequency.
    *  @param double w Width or decay length depending if the packet shape model is Gaussian or Exponential.
    *  @ingroup Circuit_distin
    *  @see qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    *  @see emitter();
    */
    void def_packet(int n, double t, double f, double w);
    /**
    *  Adds to the circuit a new definition of a wave packet. <br>
    *  <b> Warning! </b> Note that when the detectors are not configured as counters the packet numeration may change on emission.
    *
    *
    *  @param int n Suggested wavepacket number.
    *  @param double t Wavepacket characteristic emission time.
    *  @param double f Wavepacket characteristic frequency.
    *  @param double w Width or decay length depending if the packet shape model is Gaussian or Exponential.
    *  @param double tp Initial phase in time units. If tp=t then the initial phase is zero at the central (Gaussian)/front(exponential) time of the packet.
    *  @param double val1 Value 1, depends on the rand configuration
    *  @param double val2 Value 2 ,depends on the rand configuration
    *  @param int rand Phase computation. <br>
    *               - 0 = Physical manual.<br>
    *                   The packet is configured with the given parameters.<br>
    *                   val1=irrelevant. <br>
    *                   val2=irrelevant. <br>
    *                   <br>
    *               - 1 = Cascade fully automatic.<br>
    *                   An extra random (exponentially distributed) emission time is calculated to automatically compute for the probabilistic delay between an exciton X
    *                   and a bi-exciton XX in a XX-X-G cascade of a QD.<br>
    *                   val1= Characteristic time txx of the previous photon in the cascade used to calculate a random additional emission time. <br>
    *                   val2= irrelevant. <br>
    *                   The random generated value is returned in the dt1 and dt2 fields of the cfg_rnd structure.<br>
    *                   <br>
    *               - 2 = Cascade automatic.<br>
    *                   Extra random emission and phase times are calculated separately to automatically compute for the probabilistic delay between an exciton X
    *                   and a bi-exciton XX in a XX-X-G cascade of a QD. In this case the phase time is calculated with an exponential distribution
    *                   that depends on T*. <br>
    *                   val1= Characteristic time txx of the previous photon in the cascade used to calculate a random additional emission time. <br>
    *                   val2= T* pure dephasing characteristic time used to calculate an extra random phase.<br>
    *                   The extra randomly generated times for delay dt1 and phase dt2 are returned in the cfg_rnd structure. <br>
    *                   <br>
    *               - 3 = Cascade manual.<br>
    *                   Same as 2 but in this case we take those extra random times from val1 and val2 assuming they are properly generated externally.<br>
    *                   val1=Extra emission time. <br>
    *                   val2=Extra emission phase time.<br>
    *                   In this case val1 and val2 are returned in the cfg_rnd structure as dt1 and dt2 fields respectively.<br>
    *                   <br>
    *               - 4 = Cascade semi-manual.<br>
    *                   Same as 2 but only for the phase time. This is useful when we want to maintain fixed a particular overlapping between two photons. <br>
    *                   val1= Irrelevant. <br>
    *                   val2=T* pure dephasing characteristic time used to calculate an extra random phase.<br>
    *                   In this case val1=0 and val2 are returned in the cfg_rnd structure as dt1 and dt2 respectively. <br>
    *                   Note that both tx and T* are needed to calculate the phase time. <br>
    *
    *
    *  @return conf_rand structure with two random times generated by this method whose meaning depends on the configuration above. They can be used to configure manually other packets with identical values.
    *  @see qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    *  @see emitter();
    *  @ingroup Circuit_distin
    */
    cfg_rnd def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand);
    /**
    *  Overlapping probability between two wave packets.
    *
    *  @param int i  Packet i number.
    *  @param int j  Packet j number.
    *  @return Probability of overlapping.
    *  @see def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand);
    *  @ingroup Circuit_distin
    */
    double emitted_vis(int i,int j);
    /**
    *  Adds an emitter to the circuit using the packet definition given by def_packet.
    *
    *  @return  Vector with the new packet numbers indexed by the old ones ( when the detectors are not configured as counters the packet numeration may change on emission).
    *  @see def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand);
    *  @see qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    *  @ingroup Circuit_distin
    */
    veci emitter();
    /**
    *  Adds an emitter to the circuit using a matrix with the parameters needed to create a packet model of the emitted photons.
    *
    *  @param int npack Number of packets.
    *  @param matd packets Matrix that contains the definition of the packets. Each columns contains the information of a single packet with its parameters in rows:
    *   - 1st row: Packet numbers.
    *   - 2nd row: Packet central times.
    *   - 3rd row: Packet central frequencies.
    *   - 4th row: Packet width or decay time.
    *   - 5th row: Packet initial phase in time units.
    +
    *  @see def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand);
    *  @see qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    *  @ingroup Circuit_distin
    */
    veci emitter(int npack, matd packets);
    /**
    *  Adds an emitter to the circuit using a packet model of the emitted photonss.
    *
    *  @param photon_mdl *mdl Packet model.
    *  @return 0 if success -1 if an error happened.
    *  @see  photon_mdl(mati pack_def, matd times, matd freq, char ckind);
    *  @ingroup Circuit_distin
    */
    int emitter(photon_mdl *mdl);
    /**
    *  Adds an emitter to the circuit using a matrix of overlapping coefficients between packets.
    *
    *  @param matc D Overlapping coefficient matrix for the wavepackets of the initial states.
    *  @ingroup Circuit_distin
    * \xrefitem know "KnowIss" "Known Issues"  When two identical full overlapping packets are present the library does not provide a correct result. A workaround is to define the overlapping as
    *  0.999999 with as many 9's as possible. This bug is inherited from the mathematical library Eigen 3.
    */
    void emitter(matc D);
    /**
    *  Increases the optical path of a channel by a quantity dt.
    *
    *  @param int  i_ch  Channel where the delay is introduced.
    *  @param double dt  Delay time introduced in the channel.
    *  @return 0 if success -1 if an error happened.
    *  @see def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand);
    *  @see qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char i_ckind);
    *  @ingroup Circuit_distin
    * \xrefitem know "KnowIss" "Known Issues" Delay can not be called if the circuit is configured to calculate explicit losses (flag loss=true in qocircuit constructor).
    */
    int delay(int ch, double dt);
    /**
    *  Increases the optical path of a channel by a quantity dt using a packet model of the emitted wavepackets.
    *
    *  @param int  i_ch  Channel where the delay is introduced.
    *  @param double dt  Delay time introduced in the channel.
    *  @param photon_mdl *mdl Packet model.
    *  @return Returns an updated packet model with the new delayed packets.
    *  @see emitter(photon_mdl *mdl);
    *  @ingroup Circuit_distin
    */
    photon_mdl *delay(int ch, double dt, photon_mdl *mdl);
    /**
    *  Creates the delay gate that increases the optical path of a channel from the Gram-Schmidt coefficients of the orthonormalized packets present in the circuit.
    *  The new delayed packets are orthonormal between themselves and the previously defined packets.
    *
    *  @param int i_ch  Channel where the delay is introduced.
    *  @param matc   T  Complex matrix with Gram-Schmidt coefficients that  contains the full information about the indistinguishability of the packets in the circuit.
    *  @param int   nT  Number of wavefunctions that will be added.
    *  @see matc emitter(matc D).
    *  @ingroup Circuit_distin
    */
    void delay(int i_ch, matc T, int nT);

    //      Detection elements
    /** @defgroup Circuit_detector Detectors
    *   @ingroup Circuit_elements
    *   List of circuit elements to define the circuit detectors.
    */
    /**
    *  Flags the channel to be ignored by outcome calculations in probability bins and density matrices.
    *
    *  @param int i_ch Detector channel.
    *  @return 0 if success -1 if an error happened.
    *  @ingroup Circuit_detector
    */
    int ignore(int i_ch);
    /**
    *  Adds a detector to a channel of the circuit.
    *
    *  @param int i_ch Detector channel.
    *  @return 0 if success -1 if an error happened.
    *  @ingroup Circuit_detector
    */
    int detector(int i_ch);
    /**
    *  Adds a detector with a post-selection condition to a channel of the circuit.
    *  Outcomes in the remaining channels are ignored by calculations in probability bins and density matrices if the number of photons in this channel is different to the one provided in cond.
    *
    *  @param int i_ch      Detector channel.
    *  @param int cond      Detection condition.<br>
    *       cond>=0:  Readings in the remaining channels are considered only by calculations in probability bins and density matrices if the number of photons in this channel are equal to cond.<br>
    *       cond=-1:  There is no condition and works as a normal detector.<br>
    *       cond=-2:  The channel is ignored by outcome calculations in probability bins and density matrices.<br>
    *  @return 0 if success -1 if an error happened.
    *  @ingroup Circuit_detector
    */
    int detector(int i_ch, int cond);
    /**
    *  Adds a general physical detector to a channel of the circuit.
    *
    *  @param int i_ch      Detector channel.
    *  @param int cond      Detection condition.<br>
    *       cond>=0:  Readings in the remaining channels are considered only by calculations in probability bins and density matrices if the number of photons in this channel are equal to cond.<br>
    *       cond=-1:  There is no condition and works as a normal detector.<br>
    *       cond=-2:  The channel is ignored by outcome calculations in probability bins and density matrices.<br>
    *  @param double eff    Efficiency of the detector.
    *  @param duble  blnk   Ratio of time in which the detector is inactive due other detections.
    *  @param double gamma  Average rate of dark counts in this channel.
    *  @return 0 if success -1 if an error happened.
    *  @ingroup Circuit_detector
    */
    int detector(int i_ch, int cond, double eff, double blnk, double gamma);
    /**
    *  Number of channels that remain without a detector.
    *
    *  @return Remaining channels without a detector.
    *  @ingroup Circuit_detector
    */
    int remdec();
    /**
    *  Adds Gaussian white noise to the output.
    *
    *  @param stdev2 Dispersion of the Gaussian distribution.
    *  @ingroup Circuit_detector
    */
    void noise(double stdev2);

    // Print functions
    /** @defgroup Circuit_print Circuit output information
    *   @ingroup Circuit
    *   Methods to print the circuit state
    */

    /**
    *  Prints the circuit matrix.
    *
    *  @param int format  Flag that controls the print style.<br>
    *                       0 = Prints numerically.<br>
    *                       N = Prints the polarization using the alphabet (H/V).<br>
    *  @ingroup Circuit_print
    */
    void prnt(int format);

    /** @defgroup Circuit_aux Circuit auxiliary methods
    *   @ingroup Circuit
    *   Auxiliary methods. Intended for internal use of the library
    */

    /**
    *  Auxiliary method to create a circuit.<br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of spatial wavepackets to be considered.
    *  @param int clock  Detector behavior configuration.<br>
    *           - 0: No clock. The detectors behave as counters. <br>
    *           - 1: Clock. The detectors are able distinguish arrival times but they are blind to frequency. <br>
    *           - 2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. <br>
    *           - 3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. <br>
    *                 <b>Warning!</b> This is an advanced configuration. It is not recommended unless the user has knowledge on how SOQCS internally works. <br>
    *
    *
    *  @param int i_R    Number of iterations in the calculation of detector dead time and dark counts effects.
    *  @param bool loss  Are the losses calculated explicitly? True/False.
    *  @param char ckind Packet model 'G': Gaussian/ 'E': Exponential
    *  @ingroup Circuit_aux
    */
    void create_circuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char ckind);
    /**
    *  Explicit computation of losses. It creates extra channels to compute the losses of each physical channel and
    *  calculates the proper coupling between the physical and the extra loss channels. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @ingroup Circuit_aux
    */
    void compute_losses();
};

