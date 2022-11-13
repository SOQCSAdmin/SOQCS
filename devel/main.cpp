/*
*
*   SOQCS TEMPLATE FILE
*
*/

#include "soqcs.h"


int main()
{
    // Create circuit
    auto example = new qodev(4, 2);
    
    // Build circuit

    /*
    *   BUILD YOUR CIRCUIT HERE
    */

    // Create a simulator
    simulator *sim= new simulator();
    // Simulate
    auto measured=sim->run(example);

    // Print measures
    cout << endl;
    cout << "Probability outcome: " << endl << endl;
    measured->prnt_bins();
    cout << endl;
}
