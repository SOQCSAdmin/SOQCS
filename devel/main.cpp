/*
*
*   SOQCS TEMPLATE FILE
*
*/

#include "soqcs.h"


int main()
{
    // Configure SOQCS
    cfg_soqcs(4);

    // Create circuit and photon bunches
    auto example = new qocircuit(2);
    auto photons= new ph_bunch(example->num_levels());
    // Build circuit

    /*
    *   BUILD YOUR CIRCUIT HERE
    */

    // Create a simulator
    simulator *sim= new simulator();
    // Simulate
    auto measured=sim->run(photons,example);

    // Print measures
    cout <<endl;
    cout << "Probability outcome: " << endl << endl;
    measured->prnt_bins();
    cout << endl;
}
