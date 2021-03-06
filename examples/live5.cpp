/** \example live5.cpp
*   \brief <b>EXAMPLE 5</b>: Boson sampling example.<br>
*   Example of boson sampling for a for a randomly generated circuit. The results are compared with the exact calculation.<br>
**   <br>
*
*   <b>Description</b>:<br>
*   In this example a random circuit of four channels is generated and two photons are created, one in each of the first two channels. An exact calculation of
*   the output is then carried out using one of the SOQCS cores/backends an also an approximated calculation by sampling using Clifford A algorithm [1]. Both results are
*   printed for comparison. <br>
*    <b>Warning!</b> The results may not be printed in the same order for the two different calculations.<br>
*   <br>
*   <br>
*   [1] Peter Clifford and Raphael Clifford. <i>The Classical Complexity of Boson Sampling</i>, <b> arXiv:1706.01260 </b> pages 146:155. <br>
*   <br>
*/


#include "soqcs.h"
#include <chrono>

int main()
{
    cout << "* Example 5:  Boson sampling example." << endl;
    cout << endl;

    // Set the maximum number of photons above the default value.
    cfg_soqcs(2);

    // Build circuit
    auto example = new qocircuit(4);
    auto photons=new ph_bunch(example->num_levels(),1);
    photons->add_photons(1,0, example);
    photons->add_photons(1,1, example);
    photons->send2circuit('G',0,example);
    example->random_circuit();
    example->detector(0);
    example->detector(1);


    // Run exact simulation and sample
    auto sim= new simulator();
    auto apdexact=sim->run(photons,example);
    auto apdsample=sim->sample(photons,example,1000000);

    // Print exact results
    cout << "Exact result:" << endl;
    apdexact->prnt_bins(example,0.001);
    cout << endl;

    // Print sampled results
    cout << "Sampled result" << endl;
    apdsample->prnt_bins(example,0.001);
    cout << endl;

}

//*********************************************************************************************************
// Copyright ? 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
