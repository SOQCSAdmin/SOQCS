/** \example live5.cpp
*   \brief <b>EXAMPLE 5</b>: Partial distinguishability in a beamsplitter.<br>
*   Partial distinguishability example of SOQCS for a simple circuit with only one beamsplitter using a non-trivial input state. <br>

*   <b>Description</b>:<br>
*   We simulate a circuit made of a single balanced beamsplitter with three photons of Gaussian shape in each of the input channels.
*   This is a | 3, 3 > input. We consider the time, frequency and width given in random adimensional units. At the output we print the
*   probability of having the outcomes | 3, 3 >, | 4, 2 >, | 5, 1> and |6, 0 >  in the two different channels depending on the delay time
*   between them. We reproduce here numerically the analytical result presented in figure 12a of ref. [1].<br>
*   <br>
*   <br>
*   [1] Malte C Tichy. <i>Interference of identical particles from entanglement to boson-sampling </i>. <b>Journal of Physics B: Atomic, Molecular and Optical Physics, 47(10):103001 </b> (2014)<br>
*   <br>
*
*   @param const int N Number of points of the plot ( N=100 is a recommended value for good resolution)
*   @param const int occi Photon occupation of the input channels (occi,occi).
*   @param const double dtm The plot is printed in the range {-dtm,dtm} whit N points.
*   <br>
*   <br>
*
*   <b>Output example</b>:<br>
*   | Plotted output |
*   | :----------- |
*   | \image html test3c.png |
*   <br>
*   Numerical reproduction of fig. 12a of ref. [1].<br>
*   <br>
*/

#include "soqcs.h"
#include <fstream>


int main()
{
//------------------------------------------------------------------------
//  Configuration constants
    const int    N    =100;    // Number of iterations
    const int    occi =  3;     // Initial state occupation (3,3).
    const double dtm  =  4;     // Max delay
//------------------------------------------------------------------------

    cout << "* Example 5: Partial distinguishability" << endl;
    cout << endl;

    // Create objects
    auto example = new qodev(6,2,1,2,0,'G');
    auto sim= new simulator();

    // Create output file
    ofstream file;
    file.open ("Results.txt");

    //Loop occupations to be measured at output
    cout << "Calculating output" << endl;
    for(int j=0;j<=occi;j++){
        double dt=-dtm;
        int och0=occi+j;
        int och1=occi-j;
        // Set up projector
        hterm in_term;
        in_term.resize(4,2);
        in_term << 0, 1,
                   H, H,
                   0, 0,
                   och0, och1;

        file << "Output: | " << och0 << ", " << och1 << ">" << endl;
        // Loop dephasing times
        for(int i=0;i<N;i++){
            // Create circuit
            example->reset();
            example->add_photons(occi,0, H, 0.0, 1.0, 1.0);
            example->add_photons(occi,1, H,  dt, 1.0, 1.0);
            example->dispersion(1,dt);
            example->beamsplitter(0,1,45.0,0.0);
            example->detector(0);
            example->detector(1);

            // Run
            auto measured=sim->run(example);

            // Print result
            auto prob=measured->prob(in_term,example);
            file << real(dt) << " "  << prob << endl;

            // Advance
            dt=dt+2.0*dtm/((double) (N-1));

            // Free memory
            delete measured;
        }
        file << endl << endl;;
    }
    cout << "Finished. Output printed in Results.txt" << endl;
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
