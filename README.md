# Stochastic Optical Quantum Circuit Simulator (SOQCS) #

<p align="justify"> SOQCS is a C++ library (with a port in Python) to simulate optical circuits for quantum light modeled as Fock wavepackets. Optical circuits are defined from non-ideal basic components connected by a lossy medium. The library also provides support for non-ideal emitters and physical detectors considering detection efficiency, dead time, dark counts and noise effects. Detectors can also be configured to establish post-selection conditions on the circuit. Circuit measurements provide detection statistics in the form of probability outcomes and density matrices. </p>

**SOQCS features**:

* Building of circuits from **non-ideal physical devices** that may contain losses.
* Connecting the circuit elements through **lossy media**.
* Creating sets of **partially distinguishable photons** from their physical parameters (this is shape, frequency, width, etc).
* **Simulate the output of an optical circuit** using different cores/backends.
* Perform detection using a model of **physical detectors** that considers effects of efficiency, dead time, dark counts and noise.
* Establish **post-selection** conditions in the detector configuration.
* Basic **boson sampling support**.
* A **Python port**.


# Related Publications #
Please cite the most appropriate of these works if you make use of this library:

*  Javier Osca and Jiri Vala.  <span style="color:blue"> <i>Implementation of photon partial distinguishability in a quantum optical circuit simulation</i></span>. <br>
   **Comput. Phys. Commun. Volume 289, 108773 (2023).** <br>
*  Javier Osca and Jiri Vala.  <span style="color:blue"> <i>Implementation of a Stochastic Optical Quantum Circuit Simulator ( SOQCS ) </i></span>. <br>
   **In preparation.**

 
# 1. Requirements #

* Linux or MaxOsX operating system
* C++ Compiler
* GNU Make
* [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* **wget** or **curl** ( to automatically download Eigen with the provided scripts)


# 2. Quick Start #
# 2.1 How to build it? #
**Step 1**: Download the library.<br>
Download SOQCS library from the github folder.

**Step 2**: Decompress. <br>
Decompress the zip file with your favorite decompressor.

**Step 3**: Configure the library. <br>
The configuration script *config.sh* automatically downloads and installs Eigen3 external library and creates symbolic links within the SOQCS source tree that are 
needed to build SOQCS library. It also unzips the SOQCS HTML documentation stored in doc.zip. Note that you may need to give execution permissions to the scripts present in SOQCS root folder.

```bash
chmod 744 *.sh
./config.sh
```  

**Step 4**: Build the library. <br>
Inside the library main SOQCS folder type <i>make</i>. This will build the library and all the examples. 


# 2.2 How to use it? #
<p align="justify"> Nine examples in C++ can be found in the <i>examples</i> subfolder. They can be compiled with the whole library following the instructions above or typing <i>make</i> within the examples folder after compiling the library.
Additionally, there are also versions of those same examples in Python that can be found in the same folder as Jupyter notebooks. 
</p>

* **Example 1  [C++] [Python]**: Elementary example program to show a basic simulation in SOQCS.
* **Example 2  [C++] [Python]**: Example of a CNOT gate simulation in SOQCS.
* **Example 3  [C++] [Python]**: Example of a CSign gate simulation in SOQCS.
* **Example 4  [C++] [Python]**: An example of HOM visibility using a beamsplitter and physical detectors.
* **Example 5  [C++] [Python]**: An example of partial distinguishability.
* **Example 6  [C++] [Python]**: An example of the delay gate.
* **Example 7  [C++] [Python]**: A boson sampling example.
* **Example 8  [C++] [Python]**: A simulation of the entanglement swapping protocol. Example of use of density matrices in SOQCS.
* **Example 9  [C++] [Python]**: An example of a dielectric film simulation in SOQCS including losses.

For extended information about how to use SOQCS library in your own projects check the documentation.
# 3. Documentation #
<p align="justify"> For more details about how to program with SOQCS library, add it to your project or learn the details of the available methods and classes consult the available documentation that can be found in the root folder of SOQCS
after executing the configuration script. Otherwise, it will be found compressed in the file doc.zip. The documentation consists in two manuals, one for the C++ library that can be accessed by clicking on <b>SOQCS_cpp.html</b> and one for the Python port that can be accesed by clicking on <b>SOQCS_phy.html</b>.</p>

# 4. Authorship #
<b>Javier Osca</b> <br>
javier.oscacotarelo@mu.ie

<b>Jiri Vala</b> <br>
jiri.vala@mu.ie

# 5. License and Copyright #
Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved. This library and its related files are subject to the licence terms detailed in LICENCE.TXT .
Use of SOQCS is only permitted under the terms of the licence in [LICENCE.TXT](./LICENCE.TXT). 

# Version release history #


* Version RV1.4:

    * Added metropolis sampler
    * Solved bug introduced in RV1.3
    * Various bugs solved
    
* Version RV1.3:

    * Post-selection by polarization
    * Qubit polarization encoding
    * Use of circuits as custom gates
    * QOL improvements
    * Various bugs solved
    
* Version RV1.2:

    * Qubit codification.
    * More examples.
    * QOL improvements.
    * Various bugs solved.
    * Automated configuration.
    * Extended MacOsX support.
 
* Version RV1.1:

    * Extended python support.
    * QOL improvements
    * Various bugs solved.
    * Simplified configuration. 
    * Basic MacOsX support.
    
* Version RV1.0:

    * Framework for circuit simulation.
    * Non-idealities in the emitter and detector statistics.
    * Partial distinguishability of photons and photon shape model.
    * Parallel execution support.
    * Density matrix and fidelity measurements.
    * Losses support.
    * Basic Sampling.
    * Basic python support.
