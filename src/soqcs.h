/**************************************************************************
* @file soqcs.h
* @version 3.7
* @date 19/06/2022
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
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
 * <b>SOQCS is a C++ library (with a port in Python) to simulate optical circuits for quantum light modeled as Fock wavepackets.</b> Optical circuits are defined from non-ideal basic components connected by a lossy medium.
 * The library also provides support for non-ideal emitters and physical detectors considering detection efficiency, dead time, dark counts and noise effects. Detectors can also be configured to establish post-selection
 * conditions on the circuit. Circuit measurements provide detection statistics in the form of probability outcomes and density matrices.<br>
 *
 *  <b>SOQCS features</b>:
 *      - Building of circuits from <b>non-ideal physical devices</b> that may contain losses.
 *      - Connecting the circuit elements through <b>lossy media</b>
 *      - Creating sets of <b>partially distinguishable photons</b> from their physical parameters (this is shape, frequency, width, etc)
 *      - <b>Simulate the output of the circuit</b> using different cores/backends.
 *      - Perform detection using a model of <b>physical detectors</b> that considers effects of efficiency, dead time, dark counts and noise.
 *      - Establish <b>post-selection conditions</b> in the detector configuration.
 *      - Basic <b>boson sampling support</b>.
 *      - A <b>Python port</b> (early version, still in development).
 *
 *
 * \image html FigIntro1.png  <br> * <br>
 *
 * \section Pub Related publications
 *  Please cite the most appropriate of these works if you make use of this library:<br>
 * <br>
 *  Javier Osca and Jiri Vala. <i style="color:blue;">Implementation of photon partial distinguishability in a quantum optical circuit simulation</i>.  <b>In preparation</b>. <br>
 * <br>
 * \section StructureP Library structure
 * \image html Public.png
 * <br>
 * - <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen 3</a>: External library for matrix manipulation required by SOQCS. <b>It is not part of the SOQCS code.</b>
 * - <b>util </b>: Module of utilities.
 * - <b>qocircuit</b>: Circuit definition.
 * - <b>state</b>:  Quantum states definition.
 * - <b>photons</b>: Abstract photon definition.
 * - <b>pbins</b>: Probability bins. Output statistics.
 * - <b>dmat</b>: Density matrix.
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
 *  - <b>doc</b>: Documentation folder. The files of the HTML documentation are stored here. <br>
 *  - <b>examples</b>: Some examples of SOQCS simulations can be found here. <br>
 *  - <b>py_src</b>: Source tree of the SOQCS Python interface. <br>
 *  - <b>src</b>: Source tree of the SOQCS library. <br>
 * <br>
 *
 *   \section HistoryP Version release history
 *   Version RV1.0:
 *          - Framework for circuit simulation.<br>
 *          - Non-idealities in the emitter and detector statistics.<br>
 *          - Partial distinguishability of photons and photon shape model.<br>
 *          - Parallel execution support. <br>
 *          - Density matrix and fidelity measurements. <br>
 *          - Losses support. <br>
 *          - Basic Sampling. <br>
 *          - Basic python support. <br>
 * <br>
 *
 * \section lic License and copyright
 *  Copyright (c) 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
 *  This library and its related files are subject to the licence terms detailed in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a> .<br>
 *  Use of SOQCS is only permitted under the terms of the licence in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a>.<br>
 * <br>
 *
 *
 * \page install Compilation and installation
 * \section requisites Requirements
 *
 *  - <b>Linux/Unix operating system</b>. <br>
 *  - <b>C++ compiler</b>.
 *  - <b>GNU Make</b>.
 *  - <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page"> Eigen 3</a> library. (Follow the instruction installation below) <br>
 * <br>
 *
 * \section build How to build it?
 * \subsection step1 Step 1: Download the library.
 * Download the library from the github folder.
 * \subsection step2 Step 2: Decompress
 * Decompress the zip file with your favorite decompressor or type in the command line:
 *  @code
    unzip "soqcs zip filename".zip
 *  @endcode
 * and move the resulting library folder to the desired location.
 * \subsection step3 Step 3: Install Eigen 3.
 * SOQCS needs Eigen 3 external library to be built. Eigen 3 can be automatically downloaded and installed for SOQCS using the script in SOQCS root folder:
 * @code
    /"SOQCS root folder"/config.sh
 * @endcode
 *
 * <b>Note for advanced users</b>: The script is configured to download V3.4 of Eigen. If V3.4 becomes unavailable the script can be reconfigured for other versions changing the number in the configuration variable
 * <i>VER</i> within the script. <br>
 *
 * \subsection step3_Alt Step 3 (alternate): Install Eigen 3 manually.
 * If the automated installation of Eigen does not work Eigen 3 library can be installed manually following the next steps.
 *  - Download Eigen 3 library manually from its <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page"> homepage</a>.<br>
 *  - Decompress the library.
 *  - Copy the subfolder <i>Eigen</i> to its corresponding location in the SOQCS building source tree.
 *  @code
        cp -R /"Eigen library folder"/Eigen/ *  /"SOQCS root folder"/src/Eigen/
    @endcode
 * \subsection step4 Step 4: Build the library and all the examples.
 * Inside the SOQCS library main folder type <i>make</i>. This will build the library and all the examples. The examples may be found in the subfolder <i>/SOQCS root folder/examples</i>.
 * \subsection step4_Alt Step 4 (alternate): Build only the library.
 * Inside the library main folder go to <i>/SOQCS root folder/src</i> folder and type <i>make</i>. This will build only the library without the examples. <br>
 * \subsection Trouble Troubleshooting
 * If some of the previous steps fail it is recommended to check the configuration file <i>conf.inc</i> in SOQCS root directory. If the command of the compiler or the linker in your computer has a different name or some compilation
 * flags are troublesome in your system they can be reconfigured in this file.
 * <br>
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
 *  <b style="color:red;"> The python port is still under development</b> and it does not have yet all the same functionality than the parent C++ library.  We recommend
 *  to look at the python examples and follow the same guidelines than the C++ programs.
 *
 *  \subsection advancedstart Importing SOQCS to a project
 *  SOQCS has two build targets:
 *  - libsoqcs.a  : Is the C++ library
 *  - libsoqcs.so : Is the C++ side of the connection of the SOQCS port into Python
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
 *  Copy pysoqcs.py and libsoqcs.so to your working folder and start your programs with,
 *  @code
        import pysoqcs
 *  @endcode
 * these files can be found in the <i>/SOQCS root folder/py_src</i> folder
 * <br>
 *
 * \page coding Programming with SOQCS
 * \section guidelines Guidelines to program with SOQCS
 * The basic steps to program in SOQCS are shown next. An example of a simple SOQCS simulation using these steps can be found in live1.cpp . Check the examples section of this documentation
 * for this and more examples. <br>
 *
 * - <b>Configure SOQCS</b>.<br>
 *   Set the maximum number of photons to be used in the simulation and create the circuit.
 *  @code
        cfg_soqcs(2);
        auto circuit = new qocircuit(2);
 *   @endcode
 *  <br>
 * - <b>Create a set of photons and send them to the circuit</b>.<br>
    @code
    auto photons= new ph_bunch(circuit->num_levels());
    photons->add_photons(2, 1, circuit);
    photons->send2circuit(circuit);
 *  @endcode
 *  <br>
 *  <b>Alternatively initialize an input state and configure an emitter</b>.<br>
 *  The state can be created manually (next) or automatically with a QD model like in live4.cpp
 *  Please check the section about states in the manual to know more about the addition of terms to a state.
 *  @code
        auto input= new state(circuit->num_levels());
        in_term.resize(4,2);
        in_term << 0, 0, // Channels
                   H, V, // Polarization
                   0, 1, // Packet number
                   3, 2; // Occupation
        input->add_term(1.0,in_term,circuit);
 *  @endcode
 *  <br>
 *  Then the packets are defined manually (also like in live4.cpp) and a emitter is created in the circuit. Note that the final packet numbers may change from the user suggestion.
 *  @code
        circuit->def_packet(0,        0.0, 10000.0, 1.0);
        circuit->def_packet(1,       16.4, 10000.0, 1.0);
        circuit->emitter('G',0);
 *  @endcode
 *  <br>
 * - <b>Build the optical circuit</b>.<br>
 *   It is possible to build an optical circuit from the elements in the catalog. Check the documentation for more information about how to configure them. Every circuit has to end with the detectors.
 *  @code
        example->beamsplitter(0,1,theta,phi);
        example->detector(0);
        example->detector(1);
 *  @endcode
 *  <br>
 * - <b>Create a simulator and run the simulation</b>.<br>
 *   The simulator class is the core of this library. It converts a set of input photons into a set of outcomes and their probabilities. <br>
 *  @code
    simulator *sim= new simulator();
    auto measured=sim->run(photons,example);
 *  @endcode
 *  <br>
 *  It is also possible to run a simulation for an input state instead of a set of photons with the same syntax. In this case, the output is an output state. Detector effects may be calculated afterwards using the tools
 *  in the library.
 *  @code
    simulator *sim= new simulator();
    auto output=sim->run(input,example);
 *  @endcode
 *  Various cores/backends are available when creating the simulator. <br>
 *  <br>
 * - <b>Print and analyze the outptup </b>.<br>
 *  The most straightforward thing to do is to print an look at the output outcomes
 *  @code
        measured->prnt_bins();
    @endcode
 *  or the output state.
 *  @code
        output->prnt_state();
    @endcode
 *  There are many ways to configure the notation of this output in the library. We refer the attentive reader to our documentation. <br>
 *  <br>
 * - <b>Alternative: Use density matrices </b>.<br>
 *  It is also possible to create and print density matrices with a more complete information about the outputs and to obtain the statistics of mixed states.
 *  In live4.cpp can be found an example of the use of density matrices in SOQCS.
 * <br>
 *
 * \section Examples Examples of SOQCS programs
 * Five examples in C++ can be found in the <i>/SOQCS root folder/examples/</i> subfolder. They can be compiled with the whole library following the instructions above or
 * typing <i>make</i> within the examples folder after compiling the library.
 *
 * - Live1.cpp: Elementary example program to show a basic simulation in SOQCS. <br>
 * - Live2.cpp: An example of HOM visibility using a beamsplitter and physical detectors. <br>
 * - Live3.cpp: An example of the delay gate. <br>
 * - Live4.cpp: A simulation of the entanglement swapping protocol. Example of use of density matrices in SOQCS. <br>
 * - Live5.cpp: A boson sampling example. <br>
 *
 * Additionally, there are also versions of those same examples in Python that can be found in the same folder as Jupyter notebooks. The Python port is still in development
 * and it does not contain yet all the functionality of the parent C++ library therefore, currently not all examples in C++ have their Python counterpart.
 *
 * \section catalog Circuit elements catalog.
 * This section is a brief summary of the catalog of optical circuit elements available in SOQCS. We refer to the documentation of each one of them for details.
 * This list is not exhaustive therefore we also refer the user to the circuit class documentation in the modules section for more information.<br>
 *
 * * <b>Basic elements:</b>
 *     - <b style="color:blue;">beamsplitter</b>(int i_ch1, int i_ch2, double theta, double phi): Adds an ideal beamsplitter to the circuit. <br>
 *     - <b style="color:blue;">phase_shifter</b>(int i_ch, double phi): Adds a phase_shifter to the circuit. <br>
 *     - <b style="color:blue;">dielectric</b>(int i_ch1, int i_ch2, cmplx t, cmplx r): Adds a physical dielectric beamsplitter. <br>
 *     - <b style="color:blue;">loss</b>(int i_ch, double l): Adds a lossy medium. <br>
 *     - <b style="color:blue;">MMI2</b>(int i_ch1, int ch2): Adds a 2x2 ideal MMI beamsplitter. <br>
 *     - <b style="color:blue;">NSX</b>(int i_ch1, int i_ch2, int i_ch3): Adds a built-in NSX circuit (Post-selection has to be configured with the detectors). <br>
 *     - <b style="color:blue;">random_circuit()</b>: Adds a random circuit. <br>
 * <br>
 * * <b>Polarization elements:</b>
 *     - <b style="color:blue;">rotator</b>(int i_ch, double d_theta, double d_phi): Adds a rotator to the circuit. <br>
 *     - <b style="color:blue;">polbeamsplitter</b>(int i_ch1, int i_ch2, int P): Adds a polarizing beamsplitter. <br>
 *     - <b style="color:blue;">half</b>(int i_ch, int P, double alplha, double gamma): Adds a half-waveplate. <br>
 *     - <b style="color:blue;">quarter</b>(int i_ch, int P, double alplha, double gamma): Adds a quarter-waveplate. <br>
 * <br>
 * * <b>Detection elements:</b>
 *     - <b style="color:blue;">detector</b>(int i_ch): Adds a detector. <br>
 *     - <b style="color:blue;">detector</b>(int i_ch, int cond): Adds a conditional detection. <br>
 *     - <b style="color:blue;">detector</b>(int i_ch, int cond, double eff, double blnk, double gamma): Adds a general physical detector. <br>
 *     - <b style="color:blue;">noise</b>(double stdev2): Adds noise to the output. <br>
 * <br>
 * * <b>Emitter and distinguishability model:</b>
 *      In general is better to create the photons and send them to the circuit. But the emitter can also be configured manually.
 *     - <b style="color:blue;">def_packet</b>(int n, double t, double f, double w): Creates a new packet definition.  <br>
 *     - <b style="color:blue;">emitted_vis</b>(int i,int j): Returns the probability of two wave packets to overlap. <br>
 *     - <b style="color:blue;">emitter</b> (char ckind, int rand): Adds an emitter to the circuit using the packet definitions given by def_packet. <br>
 *
 * \section cores Simulator cores/backends
 *   The simulator can be configured to use four different cores/backends:
 *     - <b style="color:blue;">Direct method</b>: The calculation is performed similarly on how it is done analytically. This method is recommended when the <b>number of photons <=4</b><br>
 *     - <b style="color:blue;">Direct restricted method</b>: Same as the direct method but considering only output kets whose occupations by level are zero or one. The calculation is faster because less
 *       possible outputs have to be considered but restricts the circuits that can be simulated.<br>
 *     - <b style="color:blue;">Glynn method</b>: The calculation is performed using permanent calculations for each possible output ket. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented
 *       in gray code to calculate the permanent. This method is recommended when the <b>number of photons >=4</b><br>
 *     - <b style="color:blue;">Glynn restricted method</b>: Same as the Glynn method but considering only the output states with occupations by level of zero or one. This restricts but speeds up the output. <br>
 * <br>
 *     It is also possible to call the simulator in a <b style="color:blue;">"manual mode"</b> where only the amplitudes of a pre-determined list of kets are calculated if this list is provided to the simulator.
 *     It is possible to know more about how to configure these backends in the simulator manual page. Additionally SOQCS provide basic support for boson sampling with an implementation of the Clifford A algorithm [1]. <br>
 *   <br>
 *   [1] Peter Clifford and Raphael Clifford. <i>The Classical Complexity of Boson Sampling</i>, pages 146 to 155.
 *
 * @see simulator::simulator(const char* i_back,int i_mem) .
 */
