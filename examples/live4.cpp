/** \example live4.cpp
*   \brief <b>EXAMPLE 4</b>: Entanglement swapping protocol using a physical quantum dot as a non-ideal Bell emitter.<br>
*   Entanglement swapping protocol as presented in ref. [1] in which a quantum dot is used as a non-ideal Bell emitter.
*   We consider random noise errors due spin scattering, cross dephasing and fine structure splitting FSS in the emission of the photons. <br>
*
*   <br>
*   <b>Description</b>:<br>
*   In the entanglement swapping protocol two pairs of photons are emitted in entangled Bell states. One photon of each pair is sent to a beamsplitter
*   and the result at the output of the beamsplitter is measured. If one photon is detected in each output channel of the beamsplitter then the remaining
*   two photons are also entangled. However, if the two pairs are not ideal Bell states misdetections can happen and the resulting density matrix
*   of the two photons not traveling by the beamsplitter will not be the one of a pure state of a perfectly entangled pair of photons. In this simulation we
*   will perform the protocol in the same way than depicted in [1]. This is, the two photons not traveling trough the beamsplitter will arrive
*   at the detector through the same channel but at different times that here are labeled as 0 and 1. The calculation is repeated with different instances of the state of the QD emitter
*   to capture the effect on the output of the distribution of noise in the input. We consider random noise, cross-dephasing noise and FSS effects on the photons that may give
*   different imperfections to the input state in each emission. The result is compared with the analytical result for the current parameters. Note that
*   we are neglecting pure dephasing effects in this example (this is T<sub>2</sub><sup>*</sup> = 0) because they are computer costly but they can be
*   calculated in a longer time with small modifications to this code. We also use an ideal beamsplitter for simplicity. <br>
*   <br>
*   <br>
*   [1] F. Basso Basset et Al. <i>Entanglement swapping with photons generated on demand by a quantum dot.</i> <b>Phys. Rev. Lett., 123:160501</b> (2019)<br>
*   <br>
*
*   @param const int N Number of iterations to calculate the result.
*   @param const int prntn Number of iterations to print a progreess if the calculation message
*
*   <br>
*   Output: <br>
*   | SOQCS result |
*   | :----------- |
*   | \image html live4.png |
*   <br>
*
*   <br>
*   Analytical result ref. [1]: <br>
*   | Analytic result |
*   | :-------------- |
*   | \image html live4_analytic.png |
*   <br>
*   Note that in this output the outcomes are labeled by the polarization of the photons H/V, their time of arrival to the detector in parenthesis (in this case 0=first and 1=second) and their channels. In this
*   case both photons arrive in channel zero but delayed in time.<br>
*   <br>
*/

#include "soqcs.h"
#include <chrono>


//  Configuration constants
const int    N      =  10000;        // Number of iterations
const int    prntn  =   1000;        // Number of iterations to print.

int main()
{
    cout << "* Example 4: Entanglement swapping protocol" << endl;
    cout << endl;

    // Configure SOQCS
    cfg_soqcs(4);

    // Create simulation structures
    mati ch;
    ch.resize(3,4);
    auto example = new qocircuit(3,2,4,1);
    auto sim= new simulator();
    auto apd=new dmatrix();


    auto V=0.0;
    cout << "Start run of: " << N << endl;
    for(int i=0;i<N;i++){
        // After prntn iteration print run number
        if(i%prntn==0) cout << "Running:" << i << endl;

        // build circuit
        example->reset();
        example->def_packet(0,        0.0, 10000.0, 1.0);    // XX1
        example->def_packet(1,       16.4, 10000.0, 1.0);    // XX2
        example->def_packet(2,      16.71, 10000.0, 1.0);    //  X1(+dt)
        example->def_packet(3,       16.5, 10000.0, 1.0);    //  X2
        example->emitter('E',0);
        example->beamsplitter(1,2,45.0,0.0);
        example->detector(0);
        example->detector(1,1);
        example->detector(2,1);

        // QD emission configuration
        ch  << 0, 1, 0, 2, // Channels
               0, 2, 1, 3, // H Polarization packets
               0, 2, 1, 3; // P Polarization packets
        auto input=new state(example->num_levels());
        input->QD(ch,0.8,1.0,1.0,1.0,example);

        //Run
        auto output=sim->run(input,example);

        // Conditional detection
        apd->add_state(output,example);

        //Average visibility calculation
        V=V+example->emitted_vis(2,3);

        delete output;
        delete input;
    }
    cout << "End run" << endl;
    cout << endl;

    // Calculate and print average visibility
    cout << "Print visibility:" << endl;
    cout << "V: " << V/(double)N << endl << endl;

    // Print matrix
    cout << "Print matrix:" << endl<< endl;
    example->set_prnt_flag(2);
    apd->normalize();
    auto partial=apd->calc_measure(example);
    partial->prnt_mtx(0.01,example);
}

//*********************************************************************************************************
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
