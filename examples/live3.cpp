/** \example live3.cpp
*   \brief <b>EXAMPLE 3</b>: Simulation of a delay in the middle of a circuit.<br>
*   Demonstration of the use of a delay as a circuit element.<br>
*   <br>
*
*   <b>Description</b>:<br>
*   We consider a circuit made of two ideal balanced beamsplitters with two photons of exponential shape in each of the input channels as a theoretical representation of the two photons interference
*   experiment reported in ref.[1]. We consider a delay dt in one of the channels between the two beamsplitters and we print the probability of these two photons to
*   be measured at different times in the circuit output. The result is printed in the file Results.txt int two columns {dt,Probability}. In this case we configure an ideal
*   detector and circuit therefore the result only depends on the photon distinguishability. <br>
*   <br>
*   <br>
*   [1] Santori, C., Fattal, D., Vučković, J. et al. <i>Indistinguishable photons from a single-photon device</i>. <b>Nature 419, 594:597</b> (2002) <br>
*   <br>
*
*
*   <b>Experimental device</b>:<br>
*   | Experiment device as described in fig. 3a of ref. [1] |
*   | :---------------------------------------------------- |
*   | \image html live3_device.png width=500px              |
*   <br>
*
*   <b>Simulated circuit</b>:<br>
*   | Circuit      |
*   | :----------- |
*   | \image html live3_circuit.png width=700px |
*   <br>
*
*
*   @param const int N Number of points in the plot.
*   @param const double dtm The plot is printed in the range {-dtm,dtm} whit N points.
*
*   <b>Output example</b>:<br>
*   | ./example11b.x |
*   | :----------- |
*   | \image html live3.png |
* <br>
*/


#include "soqcs.h"
#include <fstream>

const int    N      =  101;    // Number of points
const double dtm    =  4;      // Max delay

int main(){

    cout << "* Example 3: Simulation of a delay in the middle of a circuit." << endl;
    cout << endl;

    // Configure SOQCS
    cfg_soqcs(2);

    // Create objects
    hterm in_term;
    in_term.resize(4,2);
    auto example = new qocircuit(2,1,2,1);
    auto photons=new ph_bunch(example->num_levels(),1);
    auto sim= new simulator();

    // Create output file
    ofstream file;
    file.open ("Results.txt");

    // Loop dephasing times
    auto dt=-dtm;
    for(int i=0;i<N;i++){
        // Create circuit
        example->reset();
        photons->clear();
        photons->add_photons(1,0, H, 0.0, 1.0, 1.0,example);
        photons->add_photons(1,1, H, 0.0, 1.0, 1.0,example);
        photons->send2circuit('E',0,example);
        example->beamsplitter(0,1,45.0,0.0);
        example->delay(1,dt);
        example->beamsplitter(0,1,45.0,0.0);
        example->detector(0);
        example->detector(1);

        // Run
        auto measured=sim->run(photons,example);

        // Print result
        in_term << 0, 1,
                   H, H,
                   0, 1,
                   1, 1;
        auto prob=measured->prob(in_term,example);
        file << setw(10) <<real(dt) << " "  << setw(10) << prob << endl;

        // Advance
        dt=dt+2.0*dtm/((double) (N-1));

        // Free memory
        delete measured;
    }
    cout << "Finished. Output printed in Results.txt" << endl;
}

//*********************************************************************************************************
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
