/**************************************************************************
* @file sim.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Simulator library
*
***************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************

class simulator{
//  Private variables

public:
    // Public functions
    // Management functions
    simulator();                                                                  // Create a circuit simulator. Memory quantity set by default.
    simulator(int i_mem);                                                         // Create a circuit simulator. Memory quantity explicitly.
    ~simulator();                                                                 // Destroy circuit simulator

    // Simulation execution functions
    p_bin *run(qodev *circuit, int method);                                       // Calculate output of a device
    p_bin *run(qodev *circuit, int method, int nthreads);                         // Calculate output of a device with multi-threading support
    state *run(state *istate,qocircuit *qoc, int method );                        // Calculate output state as function of the input state
    state *run(state *istate,qocircuit *qoc, int method, int nthreads );          // Calculate output state as function of the input state with multi-threading support
    state *run( state *istate, ket_list *olist, qocircuit *qoc, int method );     // Calculates the output amplitudes of the kets specified
    state *run( state *istate, ket_list *olist, qocircuit *qoc, int method, int nthreads );                    // Calculates the output amplitudes of the kets specified with multi-threading support
    p_bin *sample( qodev *input, int N );                                         // Calculate output sample of a device ( Clifford A )
    p_bin *sample( state *istate,qocircuit *qoc, int N );                         // Calculate output sample as function of the input state ( Clifford A )
    tuple<p_bin*, double> metropolis( qodev *input ,int method, int N, int Nburn, int Nthin);                  // Calculate output sample of a device ( Metropolis )
    tuple<p_bin*, double> metropolis( state *istate, qocircuit *qoc ,int method, int N, int Nburn, int Nthin); // Calculate output sample as function of the input state  ( Metropolis )

protected:
    state *DirectF(state *istate,qocircuit *qoc );                                // Direct  full distribution
    state *DirectR(state *istate,qocircuit *qoc );                                // DirectR restricted distribution
    state *GlynnF (state *istate,qocircuit *qoc );                                // Glynn  full distribution
    state *GlynnR (state *istate,qocircuit *qoc );                                // GlynnR restricted distribution
    state *RyserF( state *istate, qocircuit *qoc, int nthreads);                  // Ryser full distribution with multi-threading support
    state *RyserR( state *istate, qocircuit *qoc, int nthreads);                  // RyserR restricted distribution with multi-threading support
    state *Fast_Ryser(state *istate,qocircuit *qoc, bool F, int nthreads);        // Ryser with the distribution restricted by the post-selection condition. This method supports multi-threading.
    void aux_RyserF( state *istate, state* ostate, qocircuit *qoc, int c_nph, veci constraint, int nthreads);   // Auxiliary method to calculate RyserF
    void aux_RyserR( state *istate, state* ostate, qocircuit *qoc, int c_nph, veci constraint, int nthreads);   // Auxiliary method to calculate RyserR

    state *DirectS( state *istate, ket_list *olist, qocircuit *qoc );             // Direct single set of kets
    state *GlynnS( state *istate, ket_list *olist, qocircuit *qoc );              // Glynn single set of kets
    state *RyserS( state *istate, ket_list *olist, qocircuit *qoc, int nthreads );// Ryser single set of kets. This method supports multi-threading..

    tuple<int*, double> classical_sample(int *ilist, int nph, bool gral, bool uniform, qocircuit *qoc);         // Generate a classically calculated sample
    tuple<int*, double> uniform_general(int nph, qocircuit *qoc);                                               // Generate a uniformly distributed sample ( General    )
    tuple<int*, double>uniform_restricted(int nph, qocircuit *qoc);                                             // Generate a uniformly distributed sample ( Restricted )
};
***********************************************************************************/


#include "dmat.h"
#include <thread>
#include <future>

// Constant defaults
const int DEFSIMMEM=1000;          ///< Default simulator reserved memory for output (in bytes).


/** @defgroup Simulator
 *  Simulator classes and methods
 */

/** \class simulator
*   \brief Contains all the information the perform simulations of devices and circuits.
*
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*   @ingroup Simulator
*/


class simulator{
public:
    // Public variables
    int mem;                       ///< Memory reserved for operations


    // Public functions
    // Management functions
    /** @defgroup Simulation_management Simulator management
    *   @ingroup Simulator
    *   Creation and management of the simulation object.
    */

    /**
    *  Creates a simulator object.
    *
    *  @ingroup Simulation_management
    */
    simulator();
    /**
    *  Creates a simulator object.
    *
    *  @param int i_mem Reserved memory for the output expressed as a maximum number of terms.
    *  @ingroup Simulation_management
    * \xrefitem know "KnowIss" "Known Issues" If the quantity of memory reserved is too large it is possible to see certain slow down due the time to reserve and initialize memory.
    */
    simulator(int i_mem);
    /**
    *  Destroys the simulator object.
    *
    *  @ingroup Simulation_management
    */
    ~simulator();


    // Simulation execution functions
    /** @defgroup Simulation_execution Simulator execution
    *   @ingroup Simulator
    *   Methods to run a simulation.
    */

    /**
    *  Calculates an output outcome from a device using the selected core method and the physical detectors definitions established in that device description. ( Default single thread version ).
    *
    *  @param qodev  *circuit  Device to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">1</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">3</b> = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <b style="color:blue;">5</b> = <b>Ryser restricted</b>: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">6</b> = <b>Fast Ryser method</b>:  Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output.<br>
    *                            <b style="color:blue;">7</b> = <b>Fast Ryser restricted</b>: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. <br>
    *                            <br>
    *  @return Returns the final outcomes and their probabilities.
    *  @ingroup Simulation_execution
    */
    p_bin *run(qodev *circuit, int method);
    /**
    *  Calculates an output outcome from a device using the selected core method and the physical detectors definitions established in that device description. ( Multi-thread version ). <br>
    *  Multi-threading available for Ryser based methods. The number of threads parameter is ignored for the rest of the methods.
    *
    *  @param qodev  *circuit  Device to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">1</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">3</b> = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <b style="color:blue;">5</b> = <b>Ryser restricted</b>: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">6</b> = <b>Fast Ryser method</b>:  Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output.<br>
    *                            <b style="color:blue;">7</b> = <b>Fast Ryser restricted</b>: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. <br>
    *                            <br>
    *  @param int nthreads Number of threads.
    *  @return Returns the final outcomes and their probabilities.
    *  @ingroup Simulation_execution
    */
    p_bin *run(qodev *circuit, int method, int nthreads);
    /**
    *  Calculates an output state as a function of an input initial state using the selected core method according to the rules established by a quantum circuit. ( Default single thread version ).
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">1</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">3</b> = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <b style="color:blue;">5</b> = <b>Ryser restricted</b>: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">6</b> = <b>Fast Ryser method</b>:  Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output.<br>
    *                            <b style="color:blue;">7</b> = <b>Fast Ryser restricted</b>: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. <br>
    *                            <br>
    *  @param int nthreads Number of threads.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run(state *istate,qocircuit *qoc, int method);
    /**
    *  Calculates an output state as a function of an input initial state using the selected core method according to the rules established by a quantum circuit. ( Multi-thread version ). <br>
    *  Multi-threading available for Ryser based methods. The number of threads parameter is ignored for the rest of the methods.
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">1</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">3</b> = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <b style="color:blue;">5</b> = <b>Ryser restricted</b>: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <b style="color:blue;">6</b> = <b>Fast Ryser method</b>:  Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output.<br>
    *                            <b style="color:blue;">7</b> = <b>Fast Ryser restricted</b>: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. <br>
    *                            <br>
    *  @param int nthreads Number of threads.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run(state *istate,qocircuit *qoc, int method, int nthreads );

    /**
    *  Calculates the output amplitudes for a given list of kets as a function of an input initial state using the selected core method according to the rules established by a quantum circuit. ( Default single thread version ).<br>
    *  In this case the methods available are "Direct", Glynn" and "Ryser.
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <br>
    *  @return Returns a final state that only contains terms with kets of the provided list with the amplitudes obtained from the application of the circuit rules to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run( state *istate, ket_list *olist, qocircuit *qoc, int method );
    /**
    *  Calculates the output amplitudes for a given list of kets as a function of an input initial state using the selected core method according to the rules established by a quantum circuit. ( Multi-thread version ).<br>
    *  In this case the methods available are "Direct", Glynn" and "Ryser. <br>
    *  Multi-threading available for the Ryser method. The number of threads parameter is ignored for the rest of the methods.
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int method  Core method. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">4</b> = <b>Ryser method</b>:  The calculation is performed using permanents. We use the Ryser formula to calculate the permanents.<br>
    *                            <br>
    *  @return Returns a final state that only contains terms with kets of the provided list with the amplitudes obtained from the application of the circuit rules to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run( state *istate, ket_list *olist, qocircuit *qoc, int method, int nthreads );
    /**
    *  Sampling of a device using Clifford A algorithm. <br>
    *  <b> Proceedings of the 2018 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA). Page 146-155. SIAM Publications Library (2018). </b>  <br>
    *  <b>Warning!</b> Clifford A is defined to be used with a single input ket. Therefore neither Bell or QD initializations are recommended.
    *
    *  @param qodev  *input  Device to be sampled.
    *  @param int N Number of samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    p_bin *sample( qodev *input, int N );
    /**
    *  Sampling of a circuit using Clifford A algorithm.<br>
    *  <b> Proceedings of the 2018 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA). Page 146-155. SIAM Publications Library (2018). </b> <br>
    *  <b>Warning!</b> Clifford A is defined to be used with a single input ket. Input states with multiple kets are not recommended.
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be sampled.
    *  @param int N Number of samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    p_bin *sample( state *istate,qocircuit *qoc, int N );
    /**
    *  Sampling of a device using a metropolis algorithm. <br>
    *  <b> Nature Physics 13, 1153-1157 (2017). </b> <br>
    *  <b>Warning!</b> Metropolis defined to be used with a single input ket. Therefore neither Bell or QD initializations are recommended.
    *
    *  @param qodev  *input  Device to be sampled.
    *  @param int method Sampling method. <br>
    *                    <b style="color:blue;">0</b> = <b> f classical. </b><br>
    *                    <b style="color:blue;">1</b> = <b> g Uniform. </b><br>
    *                    <b style="color:blue;">2</b> = <b> g Classical. </b><br>
    *                    <b style="color:blue;">3</b> = <b> f classical ( Restricted ).</b> <br>
    *                    <b style="color:blue;">4</b> = <b> g Uniform  ( Restricted ). </b> <br>
    *                    <b style="color:blue;">5</b> = <b> g Classical  ( Restricted ). </b> <br>
    *  @param int N Number of samples.
    *  @param int Nburn Number of initial samples to be skipped.
    *  @param int Nthin Number of thinning samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    tuple<p_bin*, double> metropolis( qodev *input ,int method, int N, int Nburn, int Nthin);
    /**
    *  Sampling of a circuit using a metropolis algorithm. <br>
    *  <b> Nature Physics 13, 1153-1157 (2017). </b> <br>
    *  <b>Warning!</b> Metropolis is defined to be used with a single input ket. Input states with multiple kets are not recommended.
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be sampled.
    *  @param int method Sampling method. <br>
    *                    <b style="color:blue;">0</b> = <b> f classical. </b><br>
    *                    <b style="color:blue;">1</b> = <b> g Uniform. </b><br>
    *                    <b style="color:blue;">2</b> = <b> g Classical. </b><br>
    *                    <b style="color:blue;">3</b> = <b> f classical ( Restricted ).</b> <br>
    *                    <b style="color:blue;">4</b> = <b> g Uniform  ( Restricted ). </b> <br>
    *                    <b style="color:blue;">5</b> = <b> g Classical  ( Restricted ). </b> <br>
    *  @param int N Number of samples.
    *  @param int Nburn Number of initial samples to be skipped.
    *  @param int Nthin Number of thinning samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    tuple<p_bin*, double> metropolis( state *istate, qocircuit *qoc ,int method, int N, int Nburn, int Nthin);


protected:
    /** @defgroup Simulation_auxiliary Simulator auxiliary methods
    *   @ingroup Simulator
    *   Auxiliary methods to run a simulation.
    */

    /**
    *  Calculates an output state as a function of an input initial state using the Direct method for a full output distribution.
    *  This is, calculating the amplitude of probability of every possible ket with the same number of photons that the input.
    *  In the direct method the output is calculated in the same way to how it is done analytically. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    */
    state *DirectF(state *istate,qocircuit *qoc );
    /**
    *  Calculates an output state as a function of an input initial state using the Direct method for a restricted output distribution.
    *  In this case the amplitudes of probability are calculated only for kets with 0 or 1 number of photons in each level. The rest
    *  is the same than in the DirectF method. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @see DirectF(state *istate,qocircuit *qoc );
    *  @ingroup Simulation_auxiliary
    */
    state *DirectR( state *istate,qocircuit *qoc );
    /**
    *  Calculates an output state as a function of an input initial state using a permanent calculation method for a full output distribution.
    *  This is, calculating the amplitude of probability of every possible ket with the same number of photons that the input.
    *  We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    */
    state *GlynnF( state *istate,qocircuit *qoc );
    /**
    *  Calculates an output state as a function of an input initial state using a permanent calculation method for a restricted output distribution.
    *  In this case the amplitudes of probability are calculated only for kets with 0 or 1 number of photons in each level. The rest
    *  is the same than in the GlynnF method. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanent. <br>
    *  <b> Intended for internal use of the library. </b>s
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see GlynnF( state *istate,qocircuit *qoc );
    */
    state *GlynnR( state *istate,qocircuit *qoc );

    /**
    *  Calculates an output state as a function of an input initial state using a permanent calculation method for a full output distribution.
    *  This is, calculating the amplitude of probability of every possible ket with the same number of photons that the input.
    *  We use the Ryser formula implemented in gray code to calculate the permanents. This method supports multi-threading.<br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    */
    state *RyserF( state *istate, qocircuit *qoc, int nthreads);
    /**
    *  Calculates an output state as a function of an input initial state using a permanent calculation method for a restricted output distribution.
    *  In this case the amplitudes of probability are calculated only for kets with 0 or 1 number of photons in each level. The rest
    *  is the same than in the RyserF method. We use the Ryser formula implemented in gray code to calculate the permanent. This method supports multi-threading. <br>
    *  <b> Intended for internal use of the library. </b>s
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int nthreads Number of threads.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see RyserF( state *istate,qocircuit *qoc );
    */
    state *RyserR( state *istate, qocircuit *qoc, int nthreads);
    /**
    *  Calculates an output state as a function of an input initial state using a permanent calculation method for a restricted output distribution.
    *  In this case the amplitudes of probability are calculated only for kets with non-zero contribution to post-selection. We use the Ryser formula
    *  implemented in gray code to calculate the permanent. This method supports multi-threading. <br>
    *  <b> Intended for internal use of the library. </b>s
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param bool F. True: The full distribution for the non-post selected channels is calculated. False: Only those kets with at most one photon by mode.
    *  @param int nthreads Number of threads.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see RyserF( state *istate,qocircuit *qoc );
    */
    state *Fast_Ryser(state *istate,qocircuit *qoc, bool F, int nthreads);
    /**
    *  Auxiliary method to calculate the full output distribution using the Ryser formula.
    *  <b> Intended for internal use of the library. </b>s
    *
    *  @param state     *istate Initial state.
    *  @param state     *ostate Output state. <b> Warning! this is an output variable </b>
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int c_nph  Constrained number of photons. This is, the number of photons used to describe the post-selection condition.
    *  @param veci c_nph  List of occupation by mode representing a post-selection condition. If the occupation of a mode is a negative value then the mode is not being post-selected.
    *  @param int nthreads Number of threads.
    *  @ingroup Simulation_auxiliary
    *  @see state *Fast_Ryser(state *istate,qocircuit *qoc, bool F, int nthreads);
    */
    void aux_RyserF( state *istate, state* ostate, qocircuit *qoc, int c_nph, veci constraint, int nthreads);
    /**
    *  Auxiliary method to calculate the restricted output distribution using the Ryser formula.
    *  <b> Intended for internal use of the library. </b>s
    *
    *  @param state     *istate Initial state.
    *  @param state     *ostate Output state. <b> Warning! this is an output variable </b>
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int c_nph  Constrained number of photons. This is, the number of photons used to describe the post-selection condition.
    *  @param veci c_nph  List of occupation by mode representing a post-selection condition. If the occupation of a mode is a negative value then the mode is not being post-selected.
    *  @param int nthreads Number of threads.
    *  @ingroup Simulation_auxiliary
    *  @see state *Fast_Ryser(state *istate,qocircuit *qoc, bool F, int nthreads);
    */
    void aux_RyserR( state *istate, state* ostate, qocircuit *qoc, int c_nph, veci constraint, int nthreads);

    /**
    *  Calculates the output amplitudes of the specified list of kets as a function of an input initial state according to the rules established by a quantum circuit using the Direct method.
    *  In the direct method we calculate the output in the same way we would do it analytically. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see DirectF(state *istate,qocircuit *qoc );
    */
    state *DirectS( state *istate, ket_list *olist, qocircuit *qoc );
    /**
    *  Calculates the output amplitudes of the specified list of kets as a function of an input initial state according to the rules established by a quantum circuit using a permanent calculation.
    *  We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanent. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see GlynnF( state *istate,qocircuit *qoc );
    */
    state *GlynnS( state *istate, ket_list *olist, qocircuit *qoc );
    /**
    *  Calculates the output amplitudes of the specified list of kets as a function of an input initial state according to the rules established by a quantum circuit using a permanent calculation.
    *  We use the Ryser formula implemented in gray code to calculate the permanent. This method supports multi-threading. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int nthreads Number of threads.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_auxiliary
    *  @see GlynnF( state *istate,qocircuit *qoc );
    */
    state *RyserS( state *istate, ket_list *olist, qocircuit *qoc, int nthreads );

    /**
    *  Calculates a sample for a circuit given an initial state assuming al photons are distinguishable. <br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param int *ilist Initial photon configuration
    *  @param int nph Total number of photons/ilist length
    *  @param bool gral Flag that configures the generator to use a general or a restricted Hilbert space true='General'/false='Restricted'
    *  @param bool uniform Flag that configures the generator to use an uniform distribution instead.  true='Uniform'/false='Classical distribution of a circuit'
    *  @param qocircuit *qoc    Circuit to be sampled classically.
    *  @return a sample and the success probability of the procedure.
    *  @ingroup Simulation_auxiliary
    *  @see metropolis( state *istate, qocircuit *qoc ,int method, int N, int Nburn, int Nthin);
    */
    tuple<int*, double> classical_sample(int *ilist, int nph, bool gral, bool uniform, qocircuit *qoc);
    /**
    *  Obtains a sample from a general Hilbert space with a uniform distribution. (Same probability for every element of the base)<br>
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param int nph Total number of photons.
    *  @param qocircuit *qoc    Circuit to which the sample is referred.
    *  @return a sample uniformly distributed and the success probability of the procedure.
    *  @ingroup Simulation_auxiliary
    *  @see classical_sample(int *ilist, int nph, bool gral, bool uniform, qocircuit *qoc);
    */
    tuple<int*, double> uniform_general(int nph, qocircuit *qoc);
    /**
    *  Obtains a sample from a restricted Hilbert space with a uniform distribution. (Same probability for every element of the base)<br>
    *  In a restricted Hilbert space we only allow maximum one photon by level
    *  <b> Intended for internal use of the library. </b>
    *
    *  @param int nph Total number of photons.
    *  @param qocircuit *qoc    Circuit to which the sample is referred.
    *  @return a sample uniformly distributed and the success probability of the procedure.
    *  @ingroup Simulation_auxiliary
    *  @see classical_sample(int *ilist, int nph, bool gral, bool uniform, qocircuit *qoc);
    */
    tuple<int*, double>uniform_restricted(int nph, qocircuit *qoc);
};





