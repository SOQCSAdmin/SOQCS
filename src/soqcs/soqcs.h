/**************************************************************************
* @file soqcs.h

* @version 6.0
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
 * Stochastic Optical Quantum Circuit Simulator (SOQCS) is a C++ and Python library which offers a framework to define, simulate and study quantum linear optical circuits in presence of various imperfections typically encountered in experiments.
 * Optical circuits are defined from non-ideal basic components connected by a lossy medium. The library also provides support for non-ideal emitters and physical detectors considering detection efficiency, dead time, dark counts and noise 
 * effects. Detectors can also be configured to establish post-selection conditions on the circuit. Circuit measurements provide detection statistics in the form of probability outcomes and density matrices. <br>
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
 * \section Pub Related publications
 *  Please cite the most appropriate of these works if you make use of this library:<br>
 *  <br>
 *  Javier Osca and Jiri Vala. <i style="color:blue;">Stochastic Optical Quantum Circuit Simulator ( SOQCS ) </i>. <br>
 *  <b> SoftwareX. Volume 25, 101603 (2024)</b>. <br>
 * <br>
 *  Javier Osca and Jiri Vala. <i style="color:blue;">Implementation of photon partial distinguishability in a quantum optical circuit simulation</i>. <br>
 *  <b> Comput. Phys. Commun. Volume 289, 108773 (2023). </b> <br>
 * <br>
 *
 *  <img src="../assets/FigIntro1.png" width=40% align="left"> <br>
 *  <div style="clear: both"></div>
 *  <br>
 *
 * \section catalog Circuit elements catalog
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
 * \section cores Simulator cores
 *    The simulator can be configured to use four different cores for exact calculation:
 *     - <b style="color:blue;">Direct</b>: The calculation is performed similarly on how it is done analytically. This method is recommended when the <b>number of photons <=4</b>.<br>
 *     - <b style="color:blue;">Glynn</b>: We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanent [1]. This method is recommended  when the  <b>number of photons >=4 and <=15 </b>.<br>
 *     - <b style="color:blue;">Ryser</b>: We use the Ryser formula [2] to calculate the permanent using parallelization in the ways described in ref [3]. This method is recommended  when the <b>number of photons >=7</b>.<br>
 *     - <b style="color:blue;">Fast Ryser</b>:  Same as the Ryser method but automatically restricting the output distribution to the non-zero contributions after post-selection.
 * <br>
 *    It is also possible to call the simulator in a <b style="color:blue;"> manual mode</b> where only the amplitudes of a pre-determined list of kets are calculated if this list is provided to the simulator. <br>
 *
 *    The simulator can also be used for sampling with two possible methods:
 *     - <b style="color:blue;">Clifford A</b>: [4]
 *     - <b style="color:blue;">Metropolis*</b>:[5]
 *   <br>
 *   <br>
 *   [1] D. G. Glynn, *The permanent of a square matrix*, European Journal of Combinatorics 31 (7) (2010) 1887?1891 <br>
 *   [2] H. J. Ryser, *Combinatorial mathematics*. Volume14. American Mathematical Society. (1963). <br>
 *   [3] P. Lundow, K. Markström, *Efficient computation of permanents, with applications to boson sampling and random matrices*, Journal of Computational Physics 455 (2022) 110990. <br>
 *   [4] Peter Clifford and Raphael Clifford. *The Classical Complexity of Boson Sampling*, Proceedings of the 2018 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA). Page 146-155. SIAM Publications Library (2018). <br>
 *   [5] Alex Neville et Al. *Classical boson sampling algorithms with superior performance to near-term experiments*, Nature Physics 13, 1153-1157 (2017) <br>
 *
 * \page install Compilation and installation
 * \section requisites Requirements
 *
 *  - <b>Linux or MacOsX operating system</b>
 *  - <b>C++ Compiler</b>
 *  - <b>GNU Make</b>
 *  - <b>ar tool</b>
 *  - <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page"> Eigen 3</a> library.
 * <br>
 * <br>
 * These are the C++ requirements. The python port requires some extra of its own.
 *
 * \section build Installation
 * \subsection step1 Step 1: Install compilation tools.
 * SOQCS is a C++ library with a Python port. Sources will be compiled automatically as part of the installation.
 * The system needs a C++ compiler, the make and ar tools and the Eigen 3 library for linear algebra. To install Eigen 3  
 * type in your command line,
 * @code
   sudo apt install libeigen3-dev
 * @endcode
 * or
 * @code
   sudo zypper install eigen3-devel
 * @endcode
 * depending on your distribution. The names of the libraries may also depend on the distribution. The rest of the tools
 * are standard compilation tools that can be found in your favourite package repository. 
 *
 * \subsection step2 Step 2: Install the library.
 * The C++ sources will be downloaded and compiled along with the Python port.
 * @code
   pip install git+https://github.com/SOQCSAdmin/SOQCS
 * @endcode
 *
 * \section start Using SOQCS C++ library. 
 *  In your programs include the line
 * @code
    #include <soqcs.h>
 * @endcode
 * and compile your program with
 * @code
        g++ your_program.cpp -lsoqcs -lpthread -O2 -std=c++17  -fopenmp -fPIC -o executable_name.x
 * @endcode
 *  Alternatively you may copy libsoqcs.a and the .h files to another folder. In this case the program has to be compiled as
 * @code
        g++ your_program.cpp -I /"library folder path" -L /"library folder path" -lsoqcs -lpthread -O2 -std=c++17  -fopenmp -fPIC -o executable_name.x
 * @endcode 
 *
 * For examples and the library API check the rest of the documentation.
 *
 * \page details_pub Software details
 * \section StructureP Library structure
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
 *  <img src="../assets/Public.png" width=40% align="left"> <br>
 *  <div style="clear: both"></div>
 *  <br>
 * 
 *  \section FolderP Source tree structure
 *  Inside the main SOQCS folder the next subfolders can be found in alphabetical order:
 *  - <b>doc</b>: Documentation folder. Here can be found the HTML files of the documentation. <br>
 *  - <b>examples</b>: Examples of SOQCS simulations can be found here. <br>
 *  - <b>src/soqcs</b>: Source tree of the SOQCS library. <br>
 * <br>
 *
 *   \section HistoryP Version release history
 *   Version <b>RV1.5</b>:
 *          - New benchmark example.
 *          - Added Ryser method to calculate permanents
 *          - Added Ryser method with parallelism as shown in: P.H. Lundow, K. Markström, Journal of Computational Physics, Volume 455, 2022,110990.
 *          - Added optimization in the amount of permanents to be calculated when post-selection is present.
 *          - Samplers now return a list of samples instead of just their probability distribution.
 *          - C++ Muti-thread server functionality now available also on Python.
 *          - C++ restriction on the output selection functionality now available also on Python.
 *          - New packaging. pip install should be enough to install the library.
 *          - Improved documentation.
 *
 *   Version <b>RV1.4.1</b>:
 *          - Various bugs solved 
 *
 *   Version <b>RV1.4</b>:
 *          - Added metropolis sampler. 
 *          - Solved bug introduced in V1.3. 
 *          - Various bugs solved. 
 *
 *   Version <b>RV1.3</b>:
 *          - Qubit polarization encoding. 
 *          - Post selection by polarization. 
 *          - QOL improvements. 
 *          - Various bugs solved.
 *
 *   Version <b>RV1.2</b>:
 *          - Qubit encoding. 
 *          - More examples.
 *          - QOL improvements.
 *          - Various bugs solved.
 *          - Automated configuration.
 *          - Extended MacOsX support.
 *
 *   Version <b>RV1.1</b>:
 *          - Extended Python support.
 *          - QOL improvements.
 *          - Various bugs solved.
 *          - Simplified configuration.
 *          - Basic MacOsX support.
 *
 *   Version <b>RV1.0</b>:
 *          - Framework for circuit simulation.
 *          - Non-idealities in the emitter and detector statistics.
 *          - Partial distinguishability of photons and photon shape model.
 *          - Parallel execution support.
 *          - Density matrix and fidelity measurements.
 *          - Losses support.
 *          - Basic Sampling.
 *          - Basic Python support. 
 *
 * <br>
 * 
 * \page lic License and copyright
 *  Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
 *  This library and its related files are subject to the licence terms detailed in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a> .<br>
 *  Use of SOQCS is only permitted under the terms of the licence in <a href="../assets/LICENCE.TXT">LICENCE.TXT</a>.<br>
 * <br>
 *
 */
