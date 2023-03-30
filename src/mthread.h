/**************************************************************************
* @file mthread.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Multi thread server library
*
***************************************************************************/


/***********************************************************************************
 QUICK GLANCE INDEX
************************************************************************************

class mthread{
    vector<thread> threads;         /// Vector of active threads
    vector<future<qelem>> fqueue;   /// Vector of future results
public:
    // Public functions
    // Management functions
    mthread();                                                     //  Create a multi-thread "server" for works
    ~mthread();                                                    //  Destroy a multi-thread server

    // Server work handling methods
    void send_work(state *input, simulator *sim, qocircuit *qoc);  //  Send a work to the server
    state *receive_work();                                         //  Receive a work from the server.

};

void new_thread(state *input, simulator *sim, qocircuit *qoc, promise<qelem> p);  // Work to be carried by a thread.
***********************************************************************************/


#include "sim.h"

/** @defgroup Mt_sim Multi-thread server
 *  Multi-thread server
 */

//Type definitions
/**
*  \struct qelem
*  \brief  Definition of a work
*/
struct qelem{
    state* input;       ///< Input state.
    state* output;      ///< Output state or result.
    qocircuit* qoc;     ///< Circuit to which input and output are referred.
};


/** \class mthread
*   \brief This object is a server of works. Different works can be executed in parallel if an input state, a simulator and a circuit are provided to the works server.
*
*   \author Javier Osca
*   \author Jiri Vala
*
*   \copyright Copyright &copy; 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved. <br>
*              The contents and use of this document and the related code are subject to the licence terms detailed in <a  href="../assets/LICENCE.TXT"> LICENCE.txt </a>.
*
*  @ingroup Mt_sim
*/

class mthread{
    // Private variables
    vector<thread> threads;         ///< Vector of active threads.
    vector<future<qelem>> fqueue;   ///< Vector of future results.

public:
    // Public functions
    // Management functions
    /** @defgroup Serv_management Server management
    *   @ingroup Mt_sim
    *   Creation and management of the server object.
    */

    /**
    *  Creates a server object.
    *
    *  @ingroup Serv_management
    */
    mthread();
    /**
    *  Destroys a server object.
    *
    *  @ingroup Serv_management
    */
    ~mthread();


    // Server work handling methods
    /** @defgroup Serv_handling Server handling
    *   @ingroup Mt_sim
    *   Methods to send and receive works from the multi-thread "server".
    */

    /**
    *  Sends a work to the "server".
    *
    *  @param state     *istate Initial state.
    *  @param simulator *sim    Simulator employed to perform the work.
    *  @param qocircuit *qoc    Circuit to be simulated.
    *  @ingroup Serv_handling
    */
    void send_work(state *input, simulator *sim, qocircuit *qoc);
    /**
    *  Receives a work form the "server" and returns the output state.
    *  If no work has finished then the main process is
    *  suspended until a work ends.
    *
    *  @return Returns the final state that corresponds to an application of the circuit to the
    *  initial state using a simulator.
    *  @ingroup Serv_handling
    */
    state *receive_work();

};

/**
*   Work to be carried by a thread. <br>
*   <b>Not intended to be used outside the library.</b>
*
*  @param state     *istate Initial state.
*  @param simulator *sim    Simulator employed to perform the work.
*  @param qocircuit *qoc    Circuit to be simulated.
*  @param promise<qelem> p  Index of the promised value with the task finalization.
*  @ingroup Serv_handling
*/
void new_thread(state *input, simulator *sim, qocircuit *qoc, promise<qelem> p);

