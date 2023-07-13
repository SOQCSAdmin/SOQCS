/** \example live9.cpp
*   \brief <b>EXAMPLE 9</b>: Dielectric as a balanced beamsplitter with losses<br>
*    Calculation of the output probabilities for a small dielectric film acting as a beamsplitter. The resulting values coincide with the ones reported in ref. [1].<br>
*   <br>
*
*   <b>Description</b>:<br>
*    We consider a circuit made of single dielectric thin film as studied in [1]. We reproduce numerically the results presented in Figs 2 and 3 of the same ref. [1]
*    obtained by means of analytical calculations to validate the loss model in SOQCS. Each of the figures correspond to the cases where two photons are injected in
*    the dielectric from the same direction or from opposite ones. Here, both situations are considered as two different input channels  therefore we plot the different
*    outcome probabilities as function of the transmission amplitude |t| for each of the  two cases | 2, 0 > and | 1, 1 >.<br>
*   <br>
*   <br>
*   [1] Stephen M. Barnett, John Jeffers, Alessandra Gatti, and Rodney Loudon. Quantum optics of lossy beam splitters. Phys. Rev. A, 57:2134–2145, 1998.<br>
*   <br>
*
*   @param const int och0   Occupation channel 0
*   @param const int och1   Occupation channel 1
*   @param const cmplx  C   Complex constant t=C*r
*   @param const cmplx maxt Maximum transmission amplitude
*   @param const int N      Maximum number of photons
*
*   <br>
*   <b> Test 1 </b> Output example t=ir and input |2,0>:
*   @code
    const int och0   =             2;
    const int och1   =             0;
    const cmplx  C   =            jm;
    const cmplx maxt = 1.0/sqrt(2.0);
    const int N      =           100;
*   @endcode
*
*   <br>
*   | Test 1       |
*   | :----------- |
*   | \image html live9_1.png |
*
*   <br>
*   <b> Test 2 </b> Output example t=r and input |1,1>:
*   @code
    const int och0   =       1;
    const int och1   =       1;
    const cmplx  C   =     1.0;
    const cmplx maxt = 1.0/2.0;
    const int N      =     100;
*   @endcode
*
*   <br>
*   | Test 2       |
*   | :----------- |
*   | \image html live9_2.png |
*   <br>
*   Note that conservation of probability imposes <b>|t &plusmn; r|<sup>2</sup><=1</b>.  It is for that reason that in <b>Test 2</b> we can not have physical values larger than
*   2|t|<sup>2</sup>=0.5 because |t + r|<sup>2</sup>>1 will break the previous conservation of probability condition.<br>
*   <br>
*/

#include "soqcs.h"
#include <fstream>


int main()
{
//------------------------------------------------------------------------
//  Configuration constants
    const int och0   = 2;                  // Occupation channel 0
    const int och1   = 0;                  // Occupation channel 1
    const cmplx  C   = jm;                 // Complex constant t=C*r
    const cmplx maxt = 1.0/sqrt(2.0);      // Maximum transmission amplitude
    const int N      = 100;                // Maximum number of photons
    const cmplx  dt   = maxt/(cmplx)(N-1); // Step increment in transmission amplitude of probability
//------------------------------------------------------------------------


    cout << "* Test 9b: Dielectric as a balanced beamsplitter with losses" << endl;
    cout << endl;

    // Create circuit.
    auto example = new qocircuit(2,1,1,1,0.0,0,0,true,'G');
    auto sim= new simulator();

    // Set-up the initial state
    auto input= new state(example->num_levels());
    hterm in_term;
    in_term.resize(1,2);
    in_term << och0, och1;
    input->add_term(1.0,in_term,example);

    // Create output file
    ofstream file;
    file.open ("Results.txt");

    // The first two loops are to sweep possible output configurations.
    cout << "Calculating output" << endl;
    for(int l=0;l<=och0+och1;l++){    // ... with l photons
        for(int k=0;k<=l;k++){
            cmplx t=0.0;
            hterm out_term;
            out_term.resize(2,2);
            out_term <<   0,     1,
                        l-k,     k;
            file << "Out: " << l-k << " " << k << endl;
            for(int i=0;i<(N-1);i++){
                // Create circuit
                example->reset();
                example->dielectric(0,1,t,C*t);
                example->detector(0);
                example->detector(1);

                // Simulate
                auto output=sim->run(input,example,0);

                // Compute output and obtain probability
                auto outcome=new p_bin(example->num_levels());
                outcome->add_state(output);
                auto measure=outcome->calc_measure(example);
                auto prob=measure->prob(out_term,example);

                // Print output
                file << 2*real(conj(t)*t) << " "  << prob << endl;

                // Increase step
                t=t+dt;

                // Free memory
                delete output;
                delete outcome;
                delete measure;
            }
            file << endl;
        }
    }
    cout << "Finished. Output printed in Results.txt" << endl;
}

//*********************************************************************************************************
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    *
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    *
// root directory of this source tree. Use of the source code in this file is only permitted under the    *
// terms of the licence in LICENCE.TXT.                                                                   *
//*********************************************************************************************************
