/** \example live2.cpp
*   \brief <b>EXAMPLE 2</b>: CNOT circuit.<br>
*   Example of the calculation of the output of a CNOT. <br>
*   <br>
*   <b>Description</b>:<br>
*   We use SOQCS to simulate a CNOT circuit as described in ref. [1]. <br>
*   <br>
*   [1] J L O'Brien, G J Pryde, A G White, T C Ralph, D Branning <i>Demonstration of an all-optical quantum controlled-NOT gate.</i> <b> Nature 426:264</b> (2003) <br>
*   <br>
*
*   <b>Output example</b>:<br>
*   | ./live2.x |
*   | :----------- |
*   | \image html live2.png |
*   <br>
*/


#include "soqcs.h"
#include <chrono>


int main()
{
    // Introduction message
    cout << "* Example 2: CNOT circuit." << endl;
    cout << endl;


    // Input configuration
    veci qinit;
    qinit.resize(2);
    qinit << 1, 0;

    // Qubit - channel modes mapping
    mati qmap;
    qmap.resize(2,2);
    qmap << 1 , 3,
            2 , 4;

    // Create circuit and photon bunches
    auto cnot = new qodev(2,6);
    // Add photons. Qubit initialization (Path encoding)
    cnot->qubits(qinit,qmap);
    // Build circuit
    //  First column of beamsplitters
    cnot->beamsplitter(3,4, -45.0,0.0);
    //  Second column of beamsplitters
    cnot->beamsplitter(0,1,180*acos(1.0/sqrt(3.0))/pi,0.0);
    cnot->beamsplitter(2,3,180*acos(1.0/sqrt(3.0))/pi,0.0);
    cnot->beamsplitter(4,5,180*acos(1.0/sqrt(3.0))/pi,0.0);
    //  Third column of beamsplitters
    cnot->beamsplitter(3,4, -45.0,0.0);
    //  Final column of phase shifters
    cnot->phase_shifter(1, 180);
    cnot->phase_shifter(3, 180);
    // Detectors
    cnot->detector(0,0);
    cnot->detector(1);
    cnot->detector(2);
    cnot->detector(3);
    cnot->detector(4);
    cnot->detector(5,0);

    // Create a simulator
    simulator *sim= new simulator();
    // Simulate
    auto outcome=sim->run(cnot);
    auto encoded=outcome->translate(qmap, cnot);

    // Print measures
    cout << "Input: "   << endl << endl;
    cout << "| 1, 0 >"  << endl << endl;
    cout << "Outcome: " << endl << endl;
    encoded->prnt_bins();
    cout << endl;
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
