/** \example live6.cpp
*   \brief <b>EXAMPLE 6</b>: Simulation of a delay in the middle of a circuit.<br>
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
*   | \image html live6_device.png width=500px              |
*   <br>
*
*   <b>Simulated circuit</b>:<br>
*   | Circuit      |
*   | :----------- |
*   | \image html live6_circuit.png width=700px |
*   <br>
*
*
*   @param const int N Number of points in the plot.
*   @param const double dtm The plot is printed in the range {-dtm,dtm} whit N points.
*
*   <b>Output example</b>:<br>
*   | ./live6.x |
*   | :----------- |
*   | \image html live6.png |
* <br>
*/


#include "soqcs.h"
#include <fstream>

const int    N   =  50;    // Number of iterations
const double dtm = 9.0;    // Sweep time

int main(){

    cout << "* Example 6: Simulation of a delay in the middle of a circuit." << endl;
    cout << endl;

    // Create objects
    hterm in_term;
    in_term.resize(4,2);
    vecd prob;
    prob.setZero(N);
    auto example = new qodev(2,2,1,4,4,3.0,3,0,false,'E',1);
    auto sim= new simulator();

    //Set up variables
    double delta=dtm/((double) (N-1));
    in_term << 0, 1,
               H, H,
               0, 2,
               1, 1;

    // Perform main loop calculation
    cout << "Calculating output" << endl;
    double t1=0.0002;
    for(int i=0;i<N;i++){
        double t2=0;
        for(int j=0;j<i;j++){
            // Create circuit
            example->reset();
            in_term(2,1)=example->add_photons(0,0, H,   t2, 1.0, 0.01);
            in_term(2,0)=example->add_photons(0,0, H,   t1, 1.0, 0.01);
            example->add_photons(1,0, H,  0.001, 1.0, 0.3);
            example->add_photons(1,1, H,  3.101, 1.0, 0.3);
            example->beamsplitter(0,1,45.0,0.0);
            example->delay(1);
            example->beamsplitter(0,1,45.0,0.0);
            example->detector(0);
            example->detector(1);

            // Run
            auto measured=sim->run(example);

            // Store
            double dt=t1-t2;
            int k=floor(dt/delta);
            if(k>0) prob(k)=prob(k)+measured->prob(in_term,example);

            // Advance
            t2=t2+delta;
            // Free memory
            delete measured;
        }
        // Advance
        t1=t1+delta;
    }

    // Create output file
    ofstream file;
    file.open ("Results.txt");

    // Write output
    double dt=0.0;
    double norm=0;
    for(int l=0;l<N;l++) norm=max(norm,prob(l));
    for(int m=0;m<N;m++){
        file << dt << " " << prob(m)/norm << endl;
        dt=dt+delta;
    }

    cout << "Finished. Output printed in Results.txt" << endl;
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
