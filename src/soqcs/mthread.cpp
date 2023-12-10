//======================================================================================================
// File mthread.cpp
//
// MULTI THREAD SERVER LIBRARY
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "mthread.h"


//----------------------------------------
//
//  Create a multi-thread "server" for works
//
//----------------------------------------
mthread::mthread(){


    sim=new simulator();
}


//----------------------------------------
//
//  Create a multi-thread "server" for works
//
//----------------------------------------
mthread::mthread(int mem){
//  int i_mem            // Number of memory positions reserved.


    sim=new simulator(mem);
}



//----------------------------------------
//
//  Destroy a multi-thread server
//
//----------------------------------------
mthread::~mthread(){


    delete sim;
}


//----------------------------------------
//
//  Send a work to the server
//
//----------------------------------------
void mthread::send_work(state *input, qocircuit *qoc, int method){
//  state     *input;       // Input state to be run
//  simulator *sim;         // Simulator employed to calculate the output from input.
//  qocircuit *qoc;         // Circuit employed to run the simulation
//  int        method;      // Simulation method
//  Variables
    promise<qelem> paux;    // Promise of a future result.
    state    *copyinput;    // Internal copy of the input state
    qocircuit*copyqoc;      // Internal copy of the circuit


    // Make internal copies of variables than can be modified between runs
    // to maintain a record of their state when this routine is called.
    // Note that the simulator usually is not modified between runs.
    copyinput=input->clone();
    copyqoc=qoc->clone();

    // Pus into the queue the promise of future result.
    fqueue.push_back(paux.get_future());
    // Create a new thread and execute it.
    threads.push_back(thread(new_thread, copyinput, copyqoc, sim, method, move(paux)));
}


//----------------------------------------
//
//  Receive a work from the server.
//  Work results are tried to be read in launch order.
//  If a work is not finished the main process is
//  suspended to wait for the task to end.
//
//----------------------------------------
state *mthread::receive_work(){
//  Values
    qelem receive;            // Element to be received as a result of a work


    // Get the first value from the vector of future results once it is ready.
    receive=fqueue.front().get();
    auto &th=threads.front();

    // Wait to the thread to end
    th.join();

    // Erase the value and the thread information from their corresponding vectors.
    fqueue.erase(fqueue.begin());
    threads.erase(threads.begin());

    // Delete the unneeded internal copies of the input state and circuit
    // once the calculation is finished.
    delete receive.input;
    delete receive.qoc;

    // Return the output state
    return receive.output;
}


//----------------------------------------
//
//  Work to be carried by a thread.
//
//----------------------------------------
void new_thread(state *input, qocircuit *qoc, simulator *sim, int method, promise<qelem> p){
//  state     *input;       // Input state to be run
//  qocircuit *qoc;         // Circuit employed to run the simulation
//  simulator *sim;         // Simulator employed to calculate the output from input.
//  int        method;      // Simulation method
//  promise<qelem> p        // Promise of future result variable
//  Variables
    qelem send;             // Element to be sent as a result


    // Calculate and compose the result to be sent once the task is finished
    send.input=input;
    send.output=sim->run(input,qoc,method);
    send.qoc=qoc;

    // Set the value of the promise to its definitive result.
    p.set_value(send);
}

