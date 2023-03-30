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
    simulator();                                    // Create a circuit simulator. Memory quantity set by default.
    simulator((int i_backend,int i_mem));           // Create a circuit simulator. Memory quantity and backend set explicitly.
    simulator(const char* i_back,int i_mem);        // Create a circuit simulator. Memory quantity and backend by name set explicitly.
    ~simulator();                                   // Destroy circuit simulator

    // Simulation execution functions
    p_bin *run(ph_bunch *input,qocircuit *qoc );                      // Calculate output outcome as obtained from physical detector from an input photon bunch
    state *run(state *istate,qocircuit *qoc );                        // Calculate output state as function of the input state
    state *run( state *istate, ket_list *olist, qocircuit *qoc );     // Calculates the output amplitudes of the kets specified
    p_bin *sample( state *istate,qocircuit *qoc, int N );             // Calculate output sample as function of the input state
    p_bin *sample( ph_bunch *input, qocircuit *qoc, int N );          // Calculate output sample as obtained from physical detector from an input photon bunch

    // SHOES. Assembler methods
    matc assemble(mati in_qdef, mati out_qdef,state* anzilla,qocircuit *qoc);          // Assemble circuit matrix in q-bit encoding
    matc init_mtx(mati qdef, state* input,qocircuit *qoc);            // Assemble initialization matrix in q-bit encoding
protected:
    state *DirectF(state *istate,qocircuit *qoc );                    // Direct  full distribution
    state *DirectR( state *istate,qocircuit *qoc );                   // DirectR restricted distribution
    state *GlynnF( state *istate,qocircuit *qoc );                    // Glynn  full distribution
    state *GlynnR( state *istate,qocircuit *qoc );                    // GlynnR restricted distribution

    state *DirectS( state *istate, ket_list *olist, qocircuit *qoc ); // Direct single set of kets
    state *GlynnS( state *istate, ket_list *olist, qocircuit *qoc );  // Glynn single set of keys
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

    // Print variables
    int backend=0;                 ///< Backend/Core selected


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
    *  @param int i_backend  Backend/core. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">0</b> = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">1</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">2</b> = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">3</b> = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <br>
    *  @param int i_mem Reserved memory for the output expressed as a maximum number of terms.
    *  @ingroup Simulation_management
    * \xrefitem know "KnowIss" "Known Issues" If the quantity of memory reserved is too large it is possible to see certain slow down due the time to reserve and initialize memory.
    */
    simulator(int i_backend,int i_mem);
    /**
    *  Creates a simulator object.
    *
    *  @param const char* i_back Name of the selected backend/core. There are four to choose:
    *                            <br>
    *                            <b style="color:blue;">"Direct"</b>  = <b>Direct method</b>: The calculation is performed similarly on how it is done analytically.<br>
    *                            <b style="color:blue;">"DirectR"</b> = <b>Direct restricted</b>: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output.<br>
    *                            <b style="color:blue;">"Glynn"</b>   = <b>Glynn method</b>:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents.<br>
    *                            <b style="color:blue;">"GlynnR"</b>  = <b>Glynn restricted</b>: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. <br>
    *                            <br>
    *  @param int i_mem Reserved memory for the output expressed as a maximum number of terms.
    *  @ingroup Simulation_management
    * \xrefitem know "KnowIss" "Known Issues" If the quantity of memory reserved is too large it is possible to see certain slow down due the time to reserve and initialize memory.
    */
    simulator(const char* i_back,int i_mem);
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
    *  Calculates an output outcome from a device using the selected backend/core and the physical detectors definitions established in that device description.
    *
    *  @param qodev  *circuit  Device to be simulated.
    *  @return Returns the final outcomes and their probabilities.
    *  @ingroup Simulation_execution
    */
    p_bin *run(qodev *circuit);
    /**
    *  Calculates an output state as a function of an input initial state using the selected backend/core according to the rules established by a quantum circuit .
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns the final state that correspond to an application of the circuit to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run(state *istate,qocircuit *qoc );

    /**
    *  Calculates the output amplitudes for a given list of kets as a function of an input initial state using the selected backend/core according to the rules established by a quantum circuit .
    *  In this case the backends available are "Direct" and "Glynn".
    *
    *  @param state     *istate Initial state.
    *  @param ket_list  *olist List of the kets whose output amplitude we intend to calculate.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @return Returns a final state that only contains terms with kets of the provided list with the amplitudes obtained from the application of the circuit rules to the
    *  initial state.
    *  @ingroup Simulation_execution
    */
    state *run( state *istate, ket_list *olist, qocircuit *qoc );
    /**
    *  Sampling of a device using Clifford A algorithm. <br>
    *  <b>Warning!</b> Clifford A is defined to be used with a single input ket. Therefore neither Bell or QD initializations are recommended.
    *
    *  @param qodev  *input  Device to be simulated.
    *  @param int N Number of samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    p_bin *sample( qodev *input, int N );
    /**
    *  Sampling of a circuit using Clifford A algorithm.<br>
    *  <b>Warning!</b> Clifford A is defined to be used with a single input ket. Input states with multiple kets are not recommended.
    *
    *  @param state     *istate Initial state.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @param int N Number of samples.
    *  @return Returns a set of probability bins with the number of samples for each state.
    *  @ingroup Simulation_execution
    */
    p_bin *sample( state *istate,qocircuit *qoc, int N );

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
};





