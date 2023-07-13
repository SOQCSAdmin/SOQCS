/** \example live3.cpp
*   \brief <b>EXAMPLE 3</b>: CSign circuit.<br>
*   Example of the calculation of the output of a CSign. <br>
*   <br>
*   <b>Description</b>:<br>
*   We use SOQCS to simulate a CSign circuit as described in ref. [1]. <br>
*   <br>
*   [1] E. Knill, R. Laflamme, G. J. Milburn, <i>A scheme for efficient quantum computation with linear </i>, <b> Nature 409 46-52</b> (2001) <br>
*   <br>
*
*   <b>Output example</b>:<br>
*   | ./live3.x |
*   | :----------- |
*   | \image html live3.png |
*   <br>
*/


#include "soqcs.h"
#include <chrono>


int main()
{
    // Introduction message
    cout << "* Example 3: CSign circuit." << endl;
    cout << endl;


    // Define the initial state directly
    auto qubit= new state(2);
    int occ1[2]={0,0};
    qubit->add_term(0.5,occ1);
    int occ2[2]={0,1};
    qubit->add_term(0.5,occ2);
    int occ3[2]={1,0};
    qubit->add_term(0.5,occ3);
    int occ4[2]={1,1};
    qubit->add_term(0.5,occ4);
    cout << "Input: "   << endl << endl;
    qubit->prnt_state(1);

    // Qubit - channel modes mapping
    mati qmap;
    qmap.resize(2,2);
    qmap << 0 , 2,
            1 , 3;

    // Define the ancilla values
    veci ancilla;
    ancilla.resize(4);
    ancilla << 1, 0, 1, 0;


    // Create circuit and photon bunches
    auto csign = new qodev(3,8);
    // Build circuit
    csign->beamsplitter(0,2,45.0,0.0);
    csign->NSX(0, 4, 5);
    csign->NSX(2, 6, 7);
    csign->beamsplitter(0,2,-45.0,0.0);
    csign->detector(0);
    csign->detector(1);
    csign->detector(2);
    csign->detector(3);
    csign->detector(4,1);
    csign->detector(5,0);
    csign->detector(6,1);
    csign->detector(7,0);

    // Create a simulator
    simulator *sim= new simulator();

    auto decoded=qubit->decode(qmap,ancilla,csign->circ);
    auto raw_state=sim->run(decoded,csign->circ,0);
    auto output=csign->apply_condition(raw_state);
    auto encoded=output->encode(qmap,csign->circ);
    encoded->normalize();

    // Print measures
    cout << "Outcome: " << endl << endl;
    encoded->prnt_state(1);
    cout << endl;
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
