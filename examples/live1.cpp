/** \example live1.cpp
*   \brief <b>EXAMPLE 1</b>: Elementary example program.<br>
*   Elementary example of SOQCS for a simple circuit with only one beamsplitter. <br>
*   <br>
*   <b>Structure</b>:<br>
*      - Configure SOQCS (maximum number of photons).<br>
*      - Create a circuit.<br>
*      - Create bunch of photons for that circuit.<br>
*      - Attach those photons as the input of the circuit.<br>
*      - Build the circuit.<br>
*      - Create a simulator to run it.<br>
*      - Run the simulation.<br>
*      - Print the output probabilities.<br>
*
*   <b>Output example</b>:<br>
*   | ./live1.x |
*   | :----------- |
*   | \image html live1.png |
*   <br>
*/


#include "soqcs.h"
#include <chrono>


int main()
{
    // Introduction message
    cout << "* Example 1: Elementary example program" << endl;
    cout << endl;

    // Create circuit and photon bunches
    auto example = new qodev(2,2);

    // Build circuit
    example->add_photons(2, 1);
    example->beamsplitter(0,1,45.0,0.0);
    example->detector(0);
    example->detector(1);

    // Create a simulator
    simulator *sim= new simulator();
    // Simulate
    auto measured=sim->run(example);

    // Print measures
    cout << "Probability outcome: " << endl << endl;
    measured->prnt_bins();
    cout << endl;
}

//*********************************************************************************************************
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
