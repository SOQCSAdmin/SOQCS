# Stochastic Optical Quantum Circuit Simulator (SOQCS) #

 <p align="justify"> Stochastic Optical Quantum Circuit Simulator (SOQCS) is a C++ and Python library which offers a framework to define, simulate and study quantum linear optical circuits in presence of various imperfections typically encountered in experiments. Optical circuits are defined from non-ideal basic components connected by a lossy medium. The library also provides support for non-ideal emitters and physical detectors considering detection efficiency, dead time, dark counts and noise effects. Detectors can also be configured to establish post-selection conditions on the circuit. Circuit measurements provide detection statistics in the form of probability outcomes and density matrices. </p>

**SOQCS features**:

* Building of circuits from **non-ideal physical devices** that may contain losses.
* Connecting the circuit elements through **lossy media**.
* Creating sets of **partially distinguishable photons** from their physical parameters (this is shape, frequency, width, etc).
* **Simulate the output of an optical circuit** using different cores.
* Perform detection using a model of **physical detectors** that considers effects of efficiency, dead time, dark counts and noise.
* Establish **post-selection** conditions in the detector configuration.
* **Boson sampling support**.
* A **Python port**.


# Related Publications #
Please cite the most appropriate of these works if you make use of this library:

*  Javier Osca and Jiri Vala.  <span style="color:blue"> <i>Stochastic Optical Quantum Circuit Simulator ( SOQCS ) </i></span>. <br>
   **SoftwareX. Volume 25, 101603 (2024)**. 
   
*  Javier Osca and Jiri Vala.  <span style="color:blue"> <i>Implementation of photon partial distinguishability in a quantum optical circuit simulation</i></span>.<br> 
   **Comput. Phys. Commun. Volume 289, 108773 (2023).**
 
# 1. Requirements #

* Linux or MacOsX operating system
* C++ Compiler
* GNU Make
* ar tool
* [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* Python 3
* matplotlib (installed automatically)
* numpy (installed automatically)


# 2. Installation #

**Step 1**: Install compilation tools 

SOQCS is a C++ library with a Python port. Sources will be compiled automatically as part of the installation.
The system needs a c++ compiler, the make and ar tools and the Eigen3 library for linear algebra. To install Eigen3 
type in your command line,

sudo apt install libeigen3-dev

or 

sudo zypper install eigen3-devel

depending on your distribution. The names of the libraries may also depend on the distribution. The rest of the tools
are standard compilation tools that can be found in your favourite package repository. 


**Step 2**: Install the library.

pip install git+https://github.com/SOQCSAdmin/SOQCS

# 3. Documentation #
For more details about how to use SOQCS library in Python consult the available documentation that can be found in https://SOQCSADmin.github.io/SOQCS/index.html

The C++ API documentation can be found in https://soqcsadmin.github.io/SOQCS/indexcpp.html

# 4. Authorship #
<b>Javier Osca</b> <br>
soqcslib@gmail.com

<b>Jiri Vala</b> <br>
jiri.vala@mu.ie

# 5. License and Copyright #
Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved. This library and its related files are subject to the licence terms detailed in LICENCE.TXT .
Use of SOQCS is only permitted under the terms of the licence in [LICENCE.TXT](./LICENCE.TXT). 

# Version release history #

* Version RV1.5:
    * New benchmark example.
    * Added Ryser method to calculate permanents
    * Added Ryser method with parallelism as shown in: J. Comp. Phys., 455 (2022), 110990
    * Added optimization in the amount of permanents to be calculated when post-selection is present.
    * Samplers now return a list of samples instead of just their probability distribution.
    * C++ Muti-thread server functionality now available also on Python.
    * C++ restriction on the output selection functionality now available also on Python.
    * New packaging. pip install should be enough to install the library.
    * Improved documentation.

* Version RV1.4.1:
    * Various bugs solved
    
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
    * Basic Python support.
