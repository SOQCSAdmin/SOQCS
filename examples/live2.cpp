/** \example live2.cpp
*   \brief <b>EXAMPLE 2</b>: HOM Visibility simulation of a 2x2 MMI beamsplitter.<br>
*   Example of a HOM visibility calculation using a physical beamsplitter and physical detectors. Losses in in the photon propagation are also considered.<br>
*   <br>
*
*   <b>Description</b>:<br>
*   We simulate a circuit made of a 2x2 MMI beamsplitter with two photons of Gaussian shape in each of the input channels and a set of detectors that work as photon counters. We consider the time, frequency and width
*   given in random adimensional units. At  the output we print the probability of having two photons in two different channels depending on the delay time between them. For delay dt=0 both photons are
*   indistinguishable and the probability at the output is zero in ideal conditions. We consider time dependent losses in one of the channels and physical detectors that incorporate effects of efficiency, detector dead time, and dark counts. Furthermore
*   we also include the effect of the presence of a white Gaussian noise over the output. The result is printed in the file Results.txt in two columns {dt,Probability}. The output is similar to the experimental
*   measurements that can be done for real MMI devices like the ones in ref. [1]. <br>
*   <br>
*   <br>
*   [1] Alberto Peruzzo, et Al. <i>Multimode quantum interference of photons in multiport integrated devices.</i> <b>Nature Communications, 2:224</b> (2011)<br>
*   <br>
*
*   @param const int N Number of points in the plot.
*   @param const double dtm The plot is printed in the range {-dtm,dtm} whit N points.
*
*   <b>Output example</b>:<br>
*   | Output       |
*   | :----------- |
*   | \image html live2.png |
*   <br>
*/


#include "soqcs.h"
#include <fstream>


//  Configuration constants
const int     N     = 100;   // Number of points
const double dtm    =   4;   // Max delay

int main(){

    cout << "* Example 2: HOM Visibility simulation of a 2x2 MMI beamsplitter." << endl;
    cout << endl;
    // Create objects
    hterm in_term;
    in_term.resize(1,2);
    auto example = new qodev(2,2,1,2,0,10000,true,'G',1);
    auto sim= new simulator();

    // Create output file
    ofstream file;
    file.open ("Results.txt");

    // Loop dephasing times
    cout << "Calculating output" << endl;
    auto dt=-dtm;
    for(int i=0;i<N;i++){
        // Create circuit
        example->reset();
        example->add_photons(1,0, H, 0.0, 1.0, 1.0);
        example->add_photons(1,1, H,  dt, 1.0, 1.0);
        example->loss(1,0.5*(dtm+dt)/(2*dtm));
        example->MMI2(0,1);
        example->detector(0, -1, 0.85, 0.1, 0.4);
        example->detector(1, -1, 0.85, 0.1, 0.4);
        example->noise(0.001);

        // Run
        auto measured=sim->run(example);

        // Print result
        in_term << 1, 1;
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
