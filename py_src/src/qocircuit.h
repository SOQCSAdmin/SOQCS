/***********************************************************************************
* @file qocircuit.h
* @version 3.7.1
* @date 18/06/2022
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
    qocircuit(int i_nch, int i_nm, int i_ns, int clock);                     // Create circuit
    qocircuit(int i_nch, int i_nm, int i_ns, int clock,int i_R, bool loss);  // Create circuit
    ~qocircuit();                                                            // Destroy circuit
    qocircuit *clone();                                                      // Copy a circuit
    void reset();                                                            // Reset circuit
    int num_levels();                                                        // Returns the number of levels of this circuit


    // Circuit elements
    //      Basic elements
    void random_circuit();                                                   // Adds a random circuit
    void NSX(int i_ch1, int i_ch2, int i_ch3);                               // Adds a built-in NSX circuit
    void beamsplitter(int i_ch1, int i_ch2, double theta, double phi);       // Adds a ideal beamsplitter to the circuit
    void dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r);                 // Adds a physical dieletric beamsplitter
    void MMI2(int i_ch1, int ch2);                                           // Adds a 2x2 ideal MMI
    void phase_shifter(int i_ch, double phi);                                // Adds a phase_shifter to the circuit
    void phase_shifter(int i_ch, cmplx t);                                   // Adds a General phase shifter with losses
    void loss(int i_ch, double l);                                           // Adds a lossy medium
    void custom_gate(mati iodef, matc mtx);                                  // Adds a custom gate

    //      Polarization elements
    void rotator(int i_ch, double d_theta, double d_phi);                    // Adds a rotator to the circuit.
    void polbeamsplitter(int i_ch1, int i_ch2, int P);                       // Adds a polarizing beamsplitter
    void waveplate(int i_ch, int P, double alplha, double gamma);            // Adds a waveplate
    void half(int i_ch, int P, double alplha, double gamma);                 // Adds a half-waveplate
    void quarter(int i_ch, int P, double alplha, double gamma);              // Adds a quarter-waveplate


    //      Detection elements
    void detector(int i_ch);                                                 // Adds a detector
    void detector(int i_ch, int cond);                                       // Adds a conditional detection
    void detector(int i_ch, int cond, double eff, double blnk, double gamma);// Adds a general physical detector
    void noise(double stdev2);                                               // Adds noise to the output


    //      Emitter and distinguishability model.
    void def_packet(int n, double t, double f, double w);                    // Creates a new packet definition
    double emitted_vis(int i,int j);                                         // Probability of two wave packets given by def_packet to overlap
    void emitter (char ckind, int rand);                                     // Adds a emitter using the packet definition given by def_packet
    void emitter (int npack, matd packets, char ckind, int rand);            // Adds a emitter with packets defined in a matrix
    void emitter(photon_mdl *mdl);                                           // Adds a emitter using a photon model
    void emitter(matc D);                                                    // Adds a emitter using a coupling matrix
    void delay(int ch, double dt);                                           // Adds a delay using the packet definition given by def_packet
    photon_mdl *delay(int ch, double dt, photon_mdl *mdl);                   // Adds a delay using a given photon model
    void add_delay(int i_ch, matc T, int nT);                                // Adds a using an extended coupling matrix


    // Print functions
    void prnt();                                                             // Print circuit matrix
    void set_prnt_flag(int flag);                                            // Sets the flag that controls the print style.


    // Auxiliary functions
    void create_circuit(int i_nch, int i_nm, int i_ns, bool loss);           // Auxiliary function to create circuits

    void compute_losses();                                                   // Auxiliary function to compute the states with less photons at the input than the output.
};


class photon_mdl{
public:
    // Public functions
    // Management functions
    photon_mdl();                                                           // Creates an empty photon model.
    photon_mdl(mati pack_def, vecd times, matd freq, char ckind, int rand); // Creates a photon model from parameters
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

// Constants
const int DF=3;              ///< Number of degrees of freedom.
const int PD=4;              ///< Packet matrix number of fields.
const int PT=2;              ///< Packet definition number of fields.
const int PTE=PT+1;          ///< Extended packet definition number of fields.
const int TD=2;              ///< Times definition number of fields.
// Alphabet for a simplified printing.
// Note that the alphabet is designed for binary choices.
// Larger spaces are possible but in this case the user hast to resort
// to a numerical output (also available )or extend the alphabet for
// better human reading.
const char pl[2]={'H','V'};  ///< Alphabet to print polarizations.

// Constant.
const int H=0;               ///< Numerical value for horizontal polarization
const int V=1;               ///< Numerical value for vertical polarization

//Group definitions
/** @defgroup Photonmdl Photon model
 *  Photon wave packet model classes and methods.
 */

/** @defgroup Circuit
 *  Circuit classes and methods
 */

/** \class photon_mdl Photon model
*   \brief  Contains all the information and methods to create and manipulate
*    the definition of a set of photon wave packets
*
*   \version 1.1
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
    mati pack_def;      ///< Ordered definition of a photon packet. Each column defines a packet where the first row is an index to the times matrix and the second an index to the freq matrix.
    matd times;         ///< Time definitions. Each column contains one definition. The first row are the relative times to the leading packet while the second the phase times (Separation convenient for calculation).
    matd freq;          ///< Frequency definitions. Each column contains one definition. The first row are frequencies while the second are widths (Gaussian) or characteristic decay times (Exponentials).
    int kind;           ///< Kind of the wave packet 0='Gaussian', 1='Exponential'.
    int rand;           ///< How the phase times have been computed. 0= No compute, 1=Same than flight times, 2=Random times (within a certain distribution)


    // Public functions
    // Management functions
    /** @defgroup Photonmdl_management Photon model management
    *   @ingroup Photonmdl
    *   Creation and management of the photon model object.
    */
    /**
    *  Creates an empty photon model.
    *   @ingroup Photonmdl_management
    */
    photon_mdl();
    /**
    *  Creates a photon model from parameters. These are, a packet definition and lists of available times and frequencies.
    *  @param mati pack_def  Packet definition. It is an integer matrix of two rows where each column contains,
    *   a wave packet. <br>
    *                   1st row = the index that refers to the central time of the packet as specified in the vector times.<br>
    *                   2nd row = the index that refers to the central frequency and width of the packet as specified in the matrix freq.<br>
    *   @param  vecd times Vector that specifies the discrete sets of times in which a packet may exist in the current simulation
    *   @param  matd freq Two row matrix with the discrete set of frequency configurations that a packet may have in the current simulation where,<br>
    *                   1st row = is the packet central frequency.<br>
    *                   2nd row = is the packet width.<br>
    *   @param char ckind  Kind of the wave packet 'G'=Gaussian, 'E'=Exponential.
    *   @param int rand    How the phase times have been computed. <br>
    *       0= No compute. <br>
    *       1=Same than flight times. <br>
    *       2=Random times (within a certain distribution). <br>
    *   @ingroup Photonmdl_management
    */
    photon_mdl(mati pack_def, vecd times, matd freq, char ckind, int rand);
    /**
    *  Destroys the object photon model.
    *   @ingroup Photonmdl_management
    */
    ~photon_mdl();
    /**
    *   Creates a copy of this photon model.
    *   @returns    A copy of the present photon model.
    *   @ingroup Photonmdl_management
    */
    photon_mdl *clone();
    /**
    *  Updates the photon model to include definitions of packets displaced a quantity dt.
    *  The update is performed doubling the number of possible wave packets. The new set of
    *  wavepackets are displaced in time a quantity dt with respect the originals. This method
    *  is intended for internal use.
    *  @param double dt Quantity of time the new set of packets are delayed with respect the old ones.
    *  @return An updated photon model.
    *  @ingroup Photonmdl_management
    */
    void update(double dt);

    // Utility functions
    /** @defgroup Photonmdl_utility Photon model utilities
    *   @ingroup Photonmdl
    *   Utilities to process the information of a photon model.
    */
    /**
    *  Coverts the photon model into a matrix description that many subroutines can use.
    *  It is a transformation in the way information is stored. This routine is intended for internal use.
    *  @return A photon packet matrix that contains the same information in a different a convenient order
    *  to do certain calculations.
    *  @ingroup Photonmdl_utility
    */
    matd create_packet_mtx();
    /**
    * Probability of two wave packets to overlap.
    *  @param int i  Packet i number.
    *  @param int j  Packet j number.
    *  @return Probability of overlapping.
    *  @see class:qocircuit void emitter(photon_mdl mdl)
    *  @ingroup Photonmdl_utility
    */
    double visibility(int i,int j);
    /**
    *  Returns the packet definition in a three row format.
    *  This is:
    *       - 1st-Row: Packet number
    *       - 2nd-Row: Index to the time matrix.
    *       - 3rd-Row: Index to the frequency matrix.
    *  @return Three row packet definition.
    *  @ingroup Photonmdl_utility
    */
    mati return_packet_def();

// Print functions
    /** @defgroup Photonmdl_print Photon model print functions
    *   @ingroup Photonmdl
    *   Functions to print on screen and debug the information in a photon model.
    */
    /**
    * Prints on screen the information stored in a photon model.
    *  @ingroup Photonmdl_print
    */
    void prnt();
    /**
    * Prints on screen the table of times of the photon model.
    *  @ingroup Photonmdl_print
    */
    void prnt_times();
    /**
    * Prints on screen the table of frequencies of the photon model.
    *  @ingroup Photonmdl_print
    */
    void prnt_freqs();
    /**
    * Prints on screen the table of packet index of the photon model.
    *  @ingroup Photonmdl_print
    */
    void prnt_packets();
};

// Auxiliary methods
/**
*  Converts an unordered packet definition into an ordered one where columns are indexed by the packet number.
*  This routine is intended for internal use of the library.
*  @param mati pack_def  A three column matrix packet definition. Where the first one is the packet index, the second one the time index and the third one the frequency index.
*  @return A photon definition where the columns are indexed by the packet index.
*  @ingroup Wavepackets
*/
mati create_packet_idx(mati pack_def);



/** \class qocircuit
*   \brief  Contains all the information to create and
*           manipulate a quantum optical circuit
*
*   \version 3.7
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
        int s;             ///< Spatial part of the wavepacket. (Integer from a list of availables).
    };

    // Public variables
    int    nlevel;          ///< Number of levels.
    int    nch;             ///< Number of channels.
    int    nm;              ///< Number of modes.
    int    ns;              ///< Number of different wavepackets possible.


    // Print variables
    int flag_prnt=3;        ///< Print flag.<br>
                            ///< 0=Numerical print.<br>
                            ///< 1=Alphabetical print (more human readable).<br>

    // Dictionary
    level *idx;             ///< Level index: index->level.
    int ***i_idx;           ///< Inverse level index: level->index
    matc   circmtx;         ///< Matrix of the circuit as a function of the level number.

    // Emitter model
    int    emiss;           ///< There is an emitter model 0=No/1=Yes
    int    npack;           ///< Number of packets
    matd   pack_list;       ///< Packet list
    matc   init_dmat;       ///< Initial matrix due to emitter
    photon_mdl *emitted;    ///< Photon model of the emitter

    // Detector model
    int    losses;          ///< Do we consider a loss model 0=No/1=Yes
    int    ncond;           ///< Post selection condition length
    int    ndetc;           ///< Number of detectors
    int    timed;           ///< There is a clock related with the detectors 0=No/1=Yes
    mati   det_def;         ///< Post-selection condition definition
    matd   det_par;         ///< Detector physical parameters
    int    R;               ///< Number of iterations to calculate detector dead time and dark-counts
    double dev;             ///< Detector standard deviation squared of gaussian noise
    double confidence;      ///< Confidence in the Cholesky approximation.


    // Public functions
    // Management functions
    /** @defgroup Circuit_management Circuit management
    *   @ingroup Circuit
    *   Creation and management of the circuit object.
    */

    /**
    *  Creates an optical circuit object.
    *  @param int i_nch  Number of channels.
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch);
    /**
    *  Creates an optical circuit object.
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of spatial wavepackets to be considered.
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns);
    /**
    *  Creates an optical circuit object.
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of spatial wavepackets to be considered.
    *  @param int clock  Number that determines if there is a clock in the circuit. This configures how the detectors behave.<br>
    *           - 0: No clock. The detectors behave as counters. <br>
    *           - 1: Clock. The detectors are able distinguish arrival times but they are blind to frequency. <br>
    *           - 2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. <br>
    *           - 3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. <br>
    *                 <b>Warning!</b> This is an advanced configuration. It is not recommended unless the user has knowledge on how SOQCS internally works. <br>
    *
    *
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns, int clock);
    /**
    *  Creates an optical circuit object.
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of spatial wavepackets to be considered.
    *  @param int clock  Number that determines if there is a clock in the circuit. This configures how the detectors behave.<br>
    *           - 0: No clock. The detectors behave as counters. <br>
    *           - 1: Clock. The detectors are able distinguish arrival times but they are blind to frequency. <br>
    *           - 2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. <br>
    *           - 3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. <br>
    *                 <b>Warning!</b> This is an advanced configuration. It is not recommended unless the user has knowledge on how SOQCS internally works. <br>
    *
    *
    *  @param int i_R    Number of iterations in the calculation of detector dead time and dark counts effects.
    *  @param bool loss  Do we use a loss models to calculate the states with less number of photons than the input True/False.
    *  @ingroup Circuit_management
    */
    qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss);
    /**
    *  Destroys an optical circuit.
    *  @ingroup Circuit_management
    */
    ~qocircuit();
    /**
    *  Copies a circuit object into a new one.
    *  @return Copy of the ciruit.
    *  @ingroup Circuit_management
    */
    qocircuit *clone();
    /**
    *  Resets circuit.
    *  It resets the circuit back to a blank state except for the already
    *  defined degrees of freedom.
    *  @ingroup Circuit_management
    */
    void reset();

    /**
    *  Returns the total number of levels of the circuit.
    *  @return Number of levels of the circuit.
    *  @ingroup Circuit_management
    */
    int num_levels();



    // Circuit elements
    /** @defgroup Circuit_elements Circuit elements
    *   @ingroup Circuit
    *   List of elements that can be added into a circuit.
    */
    //     Basic elements
    /** @defgroup Circuit_basic Basic circuit elements
    *   @ingroup Circuit_elements
    *   List of basic circuit elements
    */

    /**
    *  Adds a circuit defined by a random unitary matrix.
    *  @ingroup Circuit_basic
    */
    void random_circuit();
    /**
    *  Adds a NSX circuit element. This is the built-in version of the NSX circuit.
    *  Post-selection still has to be carried out to obtain the proper functionality.
    *  @param int i_ch1     NSX input channel 1.
    *  @param int i_ch2     NSX input channel 2.
    *  @param int i_ch3     NSX input channel 3.
    *  @ingroup Circuit_basic
    */
    void NSX(int i_ch1, int i_ch2, int i_ch3);
    /**
    *  Adds a beamsplitter to the circuit attached to channels i_ch1 and i_ch2.
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
    *  @param int i_ch1     MMI input channel 1.
    *  @param int i_ch2     MMI input channel 2.
    *  @ingroup Circuit_basic
    */
    void MMI2(int i_ch1, int ch2);
    /**
    *  Adds a phase shifter to the circuit in channel i_ch.
    *  @param int i_ch      Phase shifter input channel.
    *  @param double phi    Angle phi in degrees..
    *  @ingroup Circuit_basic
    */
    void phase_shifter(int i_ch, double phi);
    /**
    *  Adds a general phase shifter with losses to the circuit in channel i_ch.
    *  @param int i_ch      Phase shifter input channel.
    *  @param cmplx t       Transmission amplitude.
    *  @ingroup Circuit_basic
    */
    void phase_shifter(int i_ch, cmplx t);
    /**
    *  Adds a lossy medium with loss probability l to the circuit in channel i_ch.
    *  @param int i_ch      Loss medium input channel.
    *  @param double l      Loss probability.
    *  @ingroup Circuit_basic
    */
    void loss(int i_ch, double l);
    /**
    *  Adds a custom gate to the circuit.
    * @param mati iodef  List of channels and polarization that define the input to the gate.
    * @param matc U      Custom gate nxn matrix definitions.
    *  @ingroup Circuit_basic
    */
    void custom_gate(mati iodef, matc U);

    //     Polarization elements elements
    /** @defgroup Circuit_polar Polarization circuit elements
    *   @ingroup Circuit_elements
    *   List of circuit elements that act on photon polarization
    */
    /**
    *  Adds a polarization rotation device to the circuit attached to channels i_ch.
    *  @param int i_ch     Rotator input channel.
    *  @param double theta  Angle theta in degrees.
    *  @param double phi    Angle phi in degrees.
    *  @ingroup Circuit_polar
    */
    void rotator(int i_ch, double d_theta, double d_phi);
    /**
    *  Adds a polarized beamsplitter attached to channels i_ch1 and i_ch2.
    *  @param int i_ch1     Beamsplitter input channel 1.
    *  @param int i_ch2     Beamsplitter input channel 2.
    *  @param int P         Polarization to which the beamsplitter is sensitive.
    *  @ingroup Circuit_polar
    */
    void polbeamsplitter(int i_ch1, int i_ch2, int pol);
    /**
    *  Adds a general waveplate attached to channels i_ch1 and i_ch2.
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @param double d_gamma Wavenumber delay of one component with respect the other in degrees.
    *  @ingroup Circuit_polar
    */
    void waveplate(int i_ch, double d_alpha, double d_gamma);
    /**
    *  Adds a half waveplate attached to channels i_ch1 and i_ch2.
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @ingroup Circuit_polar
    */
    void half(int i_ch, double alpha);
    /**
    *  Adds a quarter waveplate attached to channels i_ch1 and i_ch2.
    *  @param int i_ch1      Waveplate input channel 1.
    *  @param int i_ch2      Waveplate input channel 2.
    *  @param double d_alpha Rotation angle in degrees.
    *  @ingroup Circuit_polar
    */
    void quarter(int i_ch, double alpha);

    //      Emitter and distinguishability model.
    /** @defgroup Circuit_distin Emission and distinguishability model
    *   @ingroup Circuit_elements
    *   List of advanced circuit elements to define photon emission and the distinguishability in the circuit.
    */

    /**
    *  Adds to the circuit a new defifinition of a wave packet. <b> Warning! </b> Note that if the circuit is timed the numeration is not guaranteed.
    *
    *  @param int n Suggested wavepacket number.
    *  @param double t Wavepacket emission time.
    *  @param double f Wavepacket frequency.
    *  @param double w Width or decay length depending if the packet shape model is Gaussian or Exponential.
    *  @ingroup Circuit_distin
    *  @see void qocircuit:: clock();
    */
    void def_packet(int n, double t, double f, double w);
    /**
    *  Overlapping probability of two wave packets defined with def_packet.
    *  @param int i  Packet i number.
    *  @param int j  Packet j number.
    *  @return Probability of overlapping.
    *  @see void def_packet(int n, double t, double f, double w);
    *  @ingroup Circuit_distin
    */
    double emitted_vis(int i,int j);
    /**
    *  Adds an emitter to the circuit using the packet definition given by def_packets.
    *
    *  @param char ckind Packet model 'G': Gaussian/ 'E': Exponential.
    *  @param int rand Initial flight phase computation. <br>
    *               0= No compute.<br>
    *               1=Same than flight times.<br>
    *               2=Random times (within a certain distribution).<br>
    *  @see void def_packet(int n, double t, double f, double w);
    *  @ingroup Circuit_distin
    */
    void emitter(char ckind, int rand);
    /**
    *  Adds an emitter to the circuit using a physical photon model of the emitted wavepackets
    *  parametrized by the elements in the matrix "packets". <br>
    *  Each columns contains the information of a single packet with its parameters in rows:
    *
    *   - 1st row: Packet numbers.
    *   - 2nd row: Packet central times.
    *   - 3rd row: Packet central frequencies
    *   - 4th row: Packet width or decay time
    *
    *  @param int npack Number of packets.
    *  @param matd packets Matrix that contains the definition of the packets.
    *  @param char ckind Packet model 'G': Gaussian/ 'E': Exponential
    *  @param int rand Initial flight phase computation. <br>
    *               0= No compute.<br>
    *               1=Same than flight times.<br>
    *               2=Random times (within a certain distribution).<br>
    *  @ingroup Circuit_distin
    */
    void emitter(int npack, matd packets, char ckind, int rand);
    /**
    *  Adds an emitter to the circuit using a physical photon model of the emitted wavepackets.
    *
    *  @param photon_mdl *mdl Photon mdl/definition.
    *  @see photon_mdl create_photon_mdl(mati pack_def, vecd times, matd freq, char ckind, int rand).
    *  @ingroup Circuit_distin
    */
    void emitter(photon_mdl *mdl);
    /**
    *  Adds an emitter to the circuit using a coupling model.
    *  The coupling model are the overlapping coefficients of a non-orthornormal base
    *  of spatial wavefunctions provided in the matrix D.
    *
    *  @param matc D Overlapping coefficient matrix for the wavefunctions of the initial states.
    *  @ingroup Circuit_distin
    * \xrefitem know "KnowIss" "Known Issues"  If we have two identical full overlapping packets the library does not provide a correct result. The overlapping in this case
    *  has to be 0.999999 with as many 9's as possible. This bug is inherited from Eigen 3.
    */
    void emitter(matc D);
    /**
    *  Adds a delay in a channel using the packet definition given by def_packets.
    *  @param int  i_ch  Channel where the delay is introduced.
    *  @param double dt  Delay time introduced in the channel (because of larger path length for example).
    *  @see void def_packet(int n, double t, double f, double w);
    *  @ingroup Circuit_distin
    * \xrefitem know "KnowIss" "Known Issues" Delay can not be called for explicit loss calculations (flag loss=true)
    *  because delay of a packet is a classical operation in nature but it still works for implicit loss calculations.
    */
    void delay(int ch, double dt);
    /**
    *  Adds a delay in a channel using a physical photon model of the emitted wavepackets.
    *  @param int  i_ch  Channel where the delay is introduced.
    *  @param double dt  Delay time introduced in the channel (because of larger path length for example).
    *  @param photon_mdl *mdl Photon mdl/definition.
    *  @return Returns a updated photon model definition.
    *  @see void emitter(photon_mdl *mdl).
    *  @ingroup Circuit_distin
    */
    photon_mdl *delay(int ch, double dt, photon_mdl *mdl);
    /**
    * Adds a delay in a channel using a visibility model.
    * We provide to the routine the Gram-Schmidt basis coefficients of the new wavefunctions
    * orthonormalized with the previous ones and themselves.
    *  @param int i_ch  Channel where the delay is introduced.
    *  @param matc   T  Complex matrix wich contains the full information about the indistinguishability model.
    *  @param int   nT  Number of wavefunctions that have been added (that are new in this update).
    *  @see matc emitter(matc D).
    *  @see matc delay(int i_ch, matc D, matc T);
    *  @ingroup Circuit_distin
    */
    void delay(int i_ch, matc T, int nT);

    //      Detection elements
    /** @defgroup Circuit_detector Detectors
    *   @ingroup Circuit_elements
    *   List of circuit elements to define the circuit detectors.
    */
    /**
    *  Adds a detector to a channel of the circuit.
    *  @param int i_ch      Detector channel.
    *  @ingroup Circuit_detector
    */
    void detector(int i_ch);
    /**
    *  Adds a detector with a detection condition to a channel of the circuit.
    *  The readings in the rest of the channels are accepted only in the number of photons in this channel is equal to cond.
    *  @param int i_ch      Detector channel.
    *  @param int cond      Detection condition. If cond -1 there is no condition and works as a normal detector.
    *  @ingroup Circuit_detector
    */
    void detector(int i_ch, int cond);
    /**
    *  Adds a general physical detector to a channel of the circuit.
    *  @param int i_ch      Detector channel.
    *  @param int cond      Detection condition.<br>
    *       cond>=0: The readings in the rest of the channels are accepted only in the number of photons in this channel is equal to cond.<br>
    *       cond<0:  There is no condition and works as a normal detector.<br>
    *  @param double eff    Efficiency of the detector.
    *  @param duble  blnk   Ratio of time in which the detector is inactive due other detections.
    *  @param double gamma  Average rate of dark counts in this channel.
    *  @ingroup Circuit_detector
    */
    void detector(int i_ch, int cond, double eff, double blnk, double gamma);
    /**
    *  Adds Gaussian white noise to the output.
    *  @ingroup Circuit_detector
    */
    void noise(double stdev2);


    // Print functions
    /** @defgroup Circuit_print Circuit output information
    *   @ingroup Circuit
    *   Methods to print the circuit state
    */

    /**
    *  Prints the circuit matrix
    *  @ingroup Circuit_print
    */
    void prnt();
    /**
    *  Sets the flag that controls the print style.
    *  @param int flag  Flag that controls the print style.<br>
    *                   0 = Prints numerically.<br>
    *                   N = Prints the polarization using the alphabet (H/V).<br>
    *  @ingroup Circuit_print
    */
    void set_prnt_flag(int flag);


    /** @defgroup Circuit_aux Circuit auxiliary methods
    *   @ingroup Circuit
    *   Auxiliary methods. Intended for internal use of the library
    */

    /**
    *  Auxiliary method to create a circuit.
    *
    *  @param int i_nch  Number of channels.
    *  @param int i_nm   Number of modes.
    *  @param int i_ns   Number of spatial wavepackets to be considered.
    *  @param bool loss  Are losses going to be calculated explicitly? True=Yes/False=No.
    *  @ingroup Circuit_aux
    */
    void create_circuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss);
    /**
    *  Adds a virtual circuit element to compute the states with less photons at the input than the output.
    *  It creates extra channels to compute the losses of each physical channel and calculates the proper coupling
    *  between the physical and the extra loss channels.
    *
    *  @ingroup Circuit_aux
    */
    void compute_losses();
};

