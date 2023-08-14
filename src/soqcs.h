/**************************************************************************
* @file soqcs.h

* @version 5.2
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
*            The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
*            root directory of this source tree. Use of the source code in this file is only permitted under the
*            terms of the licence in LICENCE.TXT.
*
* @title Stochastic Optical Quantum Circuit Simulator (SOQCS) library
*
***************************************************************************/

#include "mthread.h"


/** \mainpage Introduction
 * \section SOQCS
 * <b>SOQCS is a C++ library (with a Python port) to simulate optical circuits for quantum light modeled as Fock wavepackets.</b> Optical circuits are defined from non-ideal basic components connected by a lossy medium.
 * The library also provides support for non-ideal emitters and physical detectors considering detection efficiency, dead time, dark counts and noise effects. Detectors can also be configured to establish post-selection
 * conditions on the circuit. Circuit measurements provide detection statistics in the form of probability outcomes and density matrices.<br>
 *
 *  <b>SOQCS features</b>:
 *      - Building of circuits from <b>non-ideal physical devices</b> that may contain losses.
 *      - Connecting the circuit elements through <b>lossy media</b>
 *      - Creating sets of <b>partially distinguishable photons</b> from their physical parameters (this is shape, frequency, width, etc)
 *      - <b>Simulate the output of the circuit</b> using different core methods.
 *      - Perform detection using a model of <b>physical detectors</b> that considers effects of efficiency, dead time, dark counts and noise.
 *      - Establish <b>post-selection conditions</b> in the detector configuration.
 *      - Basic <b>boson sampling support</b>.
 *      - A <b>Python port</b>.
 *
 *
 * \image html FigIntro1.png  <br>
 * <br>
 *
 * \section Pub Related publications
 *  Please cite the most appropriate of these works if you make use of this library:<br>
 * <br>
 *  Javier Osca and Jiri Vala. <i style="color:blue;">Implementation of photon partial distinguishability in a quantum optical circuit simulation</i>. <br>
 *  <b> Comput. Phys. Commun. Volume 289, 108773 (2023). </b> <br>
 *  <br>
 *  Javier Osca and Jiri Vala. <i style="color:blue;">Implementation of a Stochastic Optical Quantum Circuit Simulator ( SOQCS ) </i>. <br>
 *  <b> arXiv:2307.06965 </b>. <br>
 * <br>
 *
 * \section StructureP Library structure
 * \image html Public.png
 * <br>
 * - <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen 3</a>: External library for matrix manipulation required by SOQCS. <b>It is not part of the SOQCS code.</b>
 * - <b>util </b>: Module of utilities.
 * - <b>qocircuit</b>: Circuit definition.
 * - <b>state</b>:  Quantum bosonic state definition.
 * - <b>qodev</b>: Quantum optical device definition. Quantum device=Input photons + Circuit definition.
 * - <b>pbins</b>: Set if probability bins. Output statistics.
 * - <b>dmat</b>: Density matrix. Output statistics.
 * - <b>simulator</b>: Simulator module. This is the core of SOQCS.
 * - <b>mthread</b>: Parallelism support for SOQCS simulations.
 * - <b>SOQCS</b>: The full C++ library
 * - <b>pySOQCS</b>: Interface for the SOQCS port into Python.
 * <br>
 * <br>
 *
 *  \section FolderP Source tree structure
 *  Inside the main SOQCS folder the next subfolders can be found in alphabetical order:
 *  - <b>devel</b>: Here can be found templates for C++ and Python simulations. <br>
 *  - <b>doc</b>: Documentation folder. Here can be found the HTML files of the documentation. <br>
 *  - <b>examples</b>: Examples of SOQCS simulations can be found here. <br>
 *  - <b>py_src</b>: Source tree of the SOQCS Python interface. <br>
 *  - <b>src</b>: Source tree of the SOQCS library. <br>
 * <br>
 *
 *   \section HistoryP Version release history
 *   Version <b>RV1.0</b>:
 *          - Framework for circuit simulation.<br>
 *          - Non-idealities in the emitter and detector statistics.<br>
 *          - Partial distinguishability of photons and photon shape model.<br>
 *          - Parallel execution support. <br>
 *          - Density matrix and fidelity measurements. <br>
 *          - Losses support. <br>
 *          - Basic Sampling. <br>
 *          - Basic python support. <br>
 *
 *   Version <b>RV1.1</b>:
 *          - Extended python support.<br>
 *          - QOL improvements.<br>
 *          - Various bugs solved.<br>
 *          - Simplified configuration.<br>
 *          - Basic MacOsX support.<br>
 *
 *   Version <b>RV1.2</b>:
 *          - Qubit encoding.<br>
 *          - More examples.<br>
 *          - QOL improvements. <br>
 *          - Various bugs solved.<br>
 *          - Automated configuration. <br>
 *          - Extended MacOsX support. <br>
 *
 *   Version <b>RV1.3</b>:
 *          - Qubit polarizarion encoding. <br>
 *          - Post selection by polarization. <br>
 *          - Use of circuits as custom gates. <br>
 *          - QOL improvements. <br>
 *          - Various bugs solved.<br>
 *
 *   Version <b>RV1.4</b>:
 *          - Added metropolis sampler. <br>
 *          - Solved bug introduced in V1.3. <br>
 *          - Various bugs solved.<br>
 * <br>
 *
 * \section lic License and copyright
 *  Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
 *  This library and its related files are subject to the licence terms detailed in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a> .<br>
 *  Use of SOQCS is only permitted under the terms of the licence in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a>.<br>
 * <br>
 *
 * \page install Compilation and installation
 * \section requisites Requirements
 *
 *  - <b>Linux/Unix or MacOsX operating system</b>. <br>
 *  - <b>C++ compiler</b>.
 *  - <b>GNU Make</b>.
 *  - <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page"> Eigen 3</a> library. (Follow instructions of installation below) <br>
 *  - <b>wget</b> or <b>curl</b> ( to automatically download Eigen with the provided scripts)
 * <br>
 *
 * \section build How to build it?
 * \subsection step1 Step 1: Download the library.
 * Download SOQCS library from the github folder.
 * \subsection step2 Step 2: Decompress.
 * Decompress the zip file with your favorite decompressor or type in the command line:
 *  @code
    unzip SOQCS-main.zip
 *  @endcode
 * and move the resulting library folder to the desired location.
 * \subsection step3 Step 3: Configure the library.
 * The configuration script <i> config.sh </i> automatically downloads and installs Eigen3 external library and creates symbolic links within the SOQCS source tree that
 * are needed to build SOQCS library. It also unzips the SOQCS HTML documentation stored in doc.zip. Note that you may need to give execution
 * permissions to the scripts present in SOQCS root folder.
 * @code
    chmod 744 *.sh
    ./config.sh
 * @endcode
 *
 * \subsection step4 Step 4: Build the library.
 * Inside the SOQCS library main folder type <i>make</i>. This will build the library and all the examples. The examples may be found in the subfolder <i>/SOQCS root folder/examples</i>.
 *
 *
 * \section start How to use it?
 * \subsection startcpp Quick start in C++
 *  To build a SOQCS simulation you can create or modify a file named <i>main.cpp</i> in the folder <i>/SOQCS root folder/devel</i>. This file has to contain:
 *  @code
        #include "soqcs.h"

        main{
            "Your simulation instructions"
        }
 *  @endcode
 *  Type <i>make</i> within that folder to build the simulation and execute it with <i>./main.x</i>. You will find a few guidelines about how to program in SOQCS below. Additionally your will
 *  find some examples in the examples folder.<br>
 *
 * \subsection startpython Quick start in Python
 *  To build a SOQCS simulation in Python you can create or modify a file named main.py in the folder <i>/SOQCS root folder/devel</i>. This file has to contain:
 *  @code
        import pysoqcs

        "Your simulation instructions"
 *  @endcode
 *  Use your preferred editor or command line client to execute it. Alternatively the library can also be used in jupyter notebooks.<br>
 *  We recommend to look at the python examples and follow the same guidelines than the C++ programs.
 *
 *  \subsection advancedstart Importing SOQCS to a project
 *  SOQCS has two build targets:
 *  - libsoqcs.a  : Is the C++ library.
 *  - libSOQCS.so : Is the C++ side of the connection of the SOQCS port into Python.
 *
 *  <br>
 *  <b>Importing SOQCS C++ library.</b> <br>
 *  In your programs include the line
 * @code
    #include "soqcs.h"
 * @endcode
 * and compile your program with
 * @code
        g++ your_program.cpp -I /"SOQCS root folder"/src -L /"SOQCS root folder"/src -lsoqcs -lpthread -O2 -std=c++17  -fopenmp -o executable_name.x
 * @endcode
 *  Alternatively you may copy libsoqcs.a and the .h files to another folder. In this case the program has to be compiled as
 * @code
        g++ your_program.cpp -I /"library folder path" -L /"library folder path" -lsoqcs -lpthread -O2 -std=c++17  -fopenmp -o executable_name.x
 * @endcode
 * <br>
 *
 * <b>Importing SOQCS port in Python </b>. <br>
 *  Copy pysoqcs.py, libSOQCS.so and the configuration file pysoqcs.conf to your working folder and start your programs with,
 *  @code
        import pysoqcs
 *  @endcode
 * these files can be found in the <i>/SOQCS root folder/py_src</i> folder. pysoqcs.conf contains the path of libSOQCS.so .
 * <br>
 *
 * \page coding Starting with SOQCS
 * There are two ways of working with SOQCS. The first one is to define a device, set up its initial photons and then measure the outcomes of that circuit. The second one implies to calculate the output state of a circuit from
 * its input state. Note that the next guidelines are given for ideal circuits. We refer the reader to the examples and the rest of the manual to see how partial distinguishability and other effects are declared.
 *
 * \section guidelines Simulating an elementary a circuit.
 * The basic steps to program in SOQCS are shown next. An example of a simple SOQCS simulation using these steps can be found in live1.cpp . Check the examples section of this documentation
 * for this and more examples. <br>
 *
 * - <b>Configure SOQCS</b>.<br>
 *   Create the circuit. A circuit is defined setting at least a maximum number of photons and the number of channels but more options are available depending the complexity of the simulation.
 *  @code
    auto circuit = new qodev(2,2);
 *  @endcode
 *  <br>
 * - <b>Create the initial set of photons (for example two photons in channel 1)</b>.<br>
    @code
    circuit->add_photons(2, 1);
 *  @endcode
 *  <br>
 *
 * - <b>Build the optical circuit</b>.<br>
 *   It is possible to build an optical circuit from the elements in the catalog. Check the documentation for more information about how to configure them. Every circuit has to end with detectors.
 *  @code
        circuit->beamsplitter(0,1,theta,phi);
        circuit->detector(0);
        circuit->detector(1);
 *  @endcode
 *  <br>
 * - <b>Create a simulator and run the simulation</b>.<br>
 *   The simulator class is the core of this library. It converts the set of input photons of a circuit  into a set of outcomes and their probabilities. <br>
 *  @code
    auto sim= new simulator();
    auto measured=sim->run(circuit);
 *  @endcode
 *  <br>
 *  Various core methods are available when creating the simulator. <br>
 *  <br>
 * - <b>Print and analyze the output </b>.<br>
 *  The most straightforward thing to do is to print an look at the output outcomes
 *  @code
        measured->prnt_bins();
    @endcode
 * <br>
 *  There are many ways to configure the notation of this output in the library. We refer the attentive reader to our documentation. <br>
 *  <br>
 * - <b>Alternative: Use density matrices </b>.<br>
 *  It is also possible to create and print density matrices with a more complete information about the outputs and to obtain the statistics of mixed states.
 *  In live8.cpp can be found an example of the use of density matrices in SOQCS.
 * <br>
 * \section guidelines2 Performing an output state calculation.
 *  Alternatively, we can simulate directly quantum bosonic states. Please, check the section about states in the manual to know more about the addition of terms to a state (see also live1.cpp).
 *  @code
        auto input= new state(circuit->num_levels());
        in_term.resize(4,2);
        in_term << 0, 0, // Channels
                   H, V, // Polarization
                   0, 1, // Packet number
                   3, 2; // Occupation
        input->add_term(1.0,in_term,circuit);
 *  @endcode
 *
 *  A simulation can be run with a similar syntax than in the example above. In this case, the output is a state. For example,
 *  @code
    auto sim= new simulator();
    auto output=sim->run(input,example);
 *  @endcode
 *  and then the output state is printed.
 *  @code
        output->prnt_state();
 *  @endcode
 *
 * \section Examples Examples of SOQCS programs
 * Nine examples in C++ can be found in the <i>/SOQCS root folder/examples/</i> subfolder. They can be compiled with the whole library following the instructions above or
 * typing <i>make</i> within the examples folder after compiling the library.
 *
 * - Live1.cpp: Elementary example program to show a basic simulation in SOQCS. <br>
 * - Live2.cpp: Example of a CNOT gate simulation in SOQCS. <br>
 * - Live3.cpp: Example of a CSign gate simulation in SOQCS. <br>
 * - Live4.cpp: An example of HOM visibility using a beamsplitter and physical detectors. <br>
 * - Live5.cpp: An example of partial distinguishability. <br>
 * - Live6.cpp: An example of the delay gate. <br>
 * - Live7.cpp: A boson sampling example. <br>
 * - Live8.cpp: A simulation of the entanglement swapping protocol. Example of use of density matrices in SOQCS. <br>
 * - Live9.cpp: An example of a dielectric film simulation in SOQCS including losses. <br>
 *
 * <b> All the examples have their corresponding version in Python</b> that can be found in the same folder as Jupyter notebooks.
 *
 * \section catalog Circuit elements catalog.
 * This section is a brief summary of the catalog of optical circuit elements available in SOQCS. We refer to the documentation of each one of them for details.
 * This list is not exhaustive therefore we also refer the user to the circuit or device classes documentation in the modules section for more information.<br>
 *
 * * <b>Basic elements:</b>
 *     - <b style="color:blue;">beamsplitter</b>(int i_ch1, int i_ch2, double theta, double phi): Adds an ideal beamsplitter to the circuit. <br>
 *     - <b style="color:blue;">phase_shifter</b>(int i_ch, double phi): Adds a phase shifter to the circuit. <br>
 *     - <b style="color:blue;">delay(ch,dt) </b>: Adds a delay in a channel. <br>
 *     - <b style="color:blue;">rewire</b>(int i_ch1, int ch2): Swaps to channels. <br>
 *     - <b style="color:blue;">dielectric</b>(int i_ch1, int i_ch2, cmplx t, cmplx r): Adds a physical dielectric beamsplitter. <br>
 *     - <b style="color:blue;">loss</b>(int i_ch, double l): Adds a lossy medium. <br>
 *     - <b style="color:blue;">MMI2</b>(int i_ch1, int ch2): Adds a 2x2 ideal MMI beamsplitter. <br>
 *     - <b style="color:blue;">NSX</b>(int i_ch1, int i_ch2, int i_ch3): Adds a built-in NSX circuit (Post-selection has to be configured with the detectors). <br>
 *     - <b style="color:blue;">random_circuit()</b>: Adds a random circuit. <br>
 * <br>
 * * <b>Polarization elements:</b>
 *     - <b style="color:blue;">rotator</b>(int i_ch, double d_theta, double d_phi): Adds a rotator to the circuit. <br>
 *     - <b style="color:blue;">pol_beamsplitter</b>(int i_ch1, int i_ch2, int P): Adds a polarizing beamsplitter. <br>
 *     - <b style="color:blue;">pol_phase_shifter</b>(int i_ch, int P, double phi): Adds a polarized phase shifter. <br>
 *     - <b style="color:blue;">pol_filter</b>(int i_ch, int P): Adds a filter that removes polarization P. <br>
 *     - <b style="color:blue;">half</b>(int i_ch, int P, double alplha, double gamma): Adds a half-waveplate. <br>
 *     - <b style="color:blue;">quarter</b>(int i_ch, int P, double alplha, double gamma): Adds a quarter-waveplate. <br>
 * <br>
 * * <b>Detection elements:</b>
 *     - <b style="color:blue;">ignore</b>(int i_ch): Ignores a channel and removes it from the output. <br>
 *     - <b style="color:blue;">detector</b>(int i_ch): Adds a detector. <br>
 *     - <b style="color:blue;">detector</b>(int i_ch, int cond): Adds a conditional detection. <br>
 *     - <b style="color:blue;">detector</b>(int i_ch, int cond, double eff, double blnk, double gamma): Adds a general physical detector. <br>
 *     - <b style="color:blue;">noise</b>(double stdev2): Adds noise to the output. <br>
 *
 * \section cores Simulator core methods
 *   The simulator can be configured to use four different core methods:
 *     - <b style="color:blue;">Direct method</b>: The calculation is performed similarly on how it is done analytically. This method is recommended when the <b>number of photons <=4</b><br>
 *     - <b style="color:blue;">Direct restricted method</b>: Same as the direct method but considering only output kets whose occupations by level are zero or one. The calculation is faster because less
 *       possible outputs have to be considered but restricts the circuits that can be simulated.<br>
 *     - <b style="color:blue;">Glynn method</b>: The calculation is performed using permanent calculations for each possible output ket. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented
 *       in gray code to calculate the permanent. This method is recommended when the <b>number of photons >=4</b><br>
 *     - <b style="color:blue;">Glynn restricted method</b>: Same as the Glynn method but considering only the output states with occupations by level of zero or one. This restricts but speeds up the output. <br>
 * <br>
 *     It is also possible to call the simulator in a <b style="color:blue;">"manual mode"</b> where only the amplitudes of a pre-determined list of kets are calculated if this list is provided to the simulator.
 *     It is possible to know more about how to configure these methods in the simulator manual page. Additionally SOQCS provide basic support for boson sampling with an implementation of the Clifford A [1] and a metropolis algorithm [2]. <br>
 *   <br>
 *   [1] Peter Clifford and Raphael Clifford. <i>The Classical Complexity of Boson Sampling</i>, Proceedings of the 2018 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA). Page 146-155. SIAM Publications Library (2018). <br>
 *   [2] Alex Neville et Al. <i>Classical boson sampling algorithms with superior performance to near-term experiments<i>, Nature Physics 13, 1153-1157 (2017)
 */
