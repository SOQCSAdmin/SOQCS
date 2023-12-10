/** \example benchlive.cpp
*   \brief <b>BENCHMARK</b>: Single permanent calculation time in SOQCS<br>
*    In this example various random circuits of N photons and 2N channels are generated. The N photons are assigned one by one to be initialized in the first N channels.
*    The amplitude to find out those photons in the same configuration at the output is calculated. The calculation of this amplitude implies the calculation of a single permanent.
*    The time required to calculate this permanent with different methods is registered and plotted as function of the number of photons and channels. The permanent is calculated using
*    the Balasubramanian/Bax/Franklin/Glynn formula "Glynn", the Ryser formula "Ryser" or the Ryser formula making use of parallelization in the way suggested in ref. [1] "Ryser 10".
*    In the last case 10 processor cores are used. <br>
*   <br>
*   [1] P.H. Lundow, K. Markstrom, <i> Efficient computation of permanents, with applications to Boson sampling and random matrices</i>, <b>Journal of Computational Physics</b>, Volume 455, 2022,110990. <br>
*   <br>
*   <b>Description</b>:<br>
*   This program measures the cost in time to execute SOQCS for an increasing number of channels and photons.
*   <br>
*   Results on an Intel I7-10750H@2.60GHz 16GB Ram:<br>
*   |Photons       |Channels      |Glynn (ms)    |Ryser (ms)    |Ryser 10 (ms) |
*   |:-------------|:-------------|:-------------|:-------------|:-------------|
*   |10            |20            |0.23          |0.30          |0.09          |
*   |11            |22            |0.23          |0.14          |0.01          |
*   |12            |24            |0.44          |0.28          |0.05          |
*   |13            |26            |1.07          |0.62          |0.09          |
*   |14            |28            |2.06          |1.26          |0.17          |
*   |15            |30            |4.29          |2.62          |0.34          |
*   |16            |32            |9.24          |9.29          |0.70          |
*   |17            |34            |18.44         |11.97         |1.45          |
*   |18            |36            |32.28         |22.61         |3.11          |
*   |19            |38            |64.16         |44.80         |6.82          |
*   |20            |40            |137.84        |95.61         |14.39         |
*   <br>
*   Note that reserving (and initializing) memory for the output state adds a constant overhead to the execution times that may affect the benchmark result.
*/


#include <soqcs.h>
#include <chrono>

int main()
{
    cout << "* SOQCS Benchmark: Full distribution calculation for a random circuit" << endl;
    cout << endl;

    //Create simulator
    auto sim= new simulator(1);

    // Print table header
    cout << left   << setw(14) << setfill(' ') << "Photons";
    cout << left   << setw(14) << setfill(' ') << "Channels";
    cout << left   << setw(14) << setfill(' ') << "Glynn (ms)";
    cout << left   << setw(14) << setfill(' ') << "Ryser (ms)";
    cout << left   << setw(14) << setfill(' ') << "Ryser 10 (ms)";
    cout << endl;

    // Main loop
    for(int nph=10;nph<=20;nph++){
        int nch=2*nph;

        //Prepare circuit
        auto circuit = new qocircuit(nch);
        circuit->random_circuit();


        // Create input state and output list with nph photons
        auto input=new state(nph,circuit->num_levels());
        auto olist=new ket_list(circuit->num_levels());
        int *occ=new int[circuit->num_levels()]();
        for(int i=0;i<nph;i++) occ[i]=1;
        input->add_term(1.0,occ);
        olist->add_ket(occ);

        // Glynn
        auto start1 = chrono::steady_clock::now();
        sim->run(input,olist,circuit,2);
        auto end1 = chrono::steady_clock::now();
        auto tGlynn=chrono::duration_cast<chrono::microseconds>(end1 - start1).count()/1000.0;

        // Ryser
        auto start2 = chrono::steady_clock::now();
        sim->run(input,olist,circuit,4);
        auto end2 = chrono::steady_clock::now();
        auto tRyser=chrono::duration_cast<chrono::microseconds>(end2 - start2).count()/1000.0;

        // Ryser 10
        auto start3 = chrono::steady_clock::now();
        sim->run(input,olist,circuit,4,10);
        auto end3 = chrono::steady_clock::now();
        auto tORyser=chrono::duration_cast<chrono::microseconds>(end3 - start3).count()/1000.0;

        // Print benchmark result
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        cout << left   << setw(14) << setfill(' ') << nph;
        cout << left   << setw(14) << setfill(' ') << nch;
        cout << left   << setw(14) << setfill(' ') << setprecision(2)<< tGlynn;
        cout << left   << setw(14) << setfill(' ') << setprecision(2)<< tRyser;
        cout << left   << setw(14) << setfill(' ') << setprecision(2)<< tORyser;
        cout << endl;
    }
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
