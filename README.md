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
* A **Python port** (early version, still in development).


# Related Publications #
Please cite the most appropriate of these works if you make use of this library:

* **<span style="color:blue"> Implementation of photon partial distinguishability in a quantum optical circuit simulation. </span>**. <i>Javier Osca and Jiri Vala</i>.  In preparation. 

# 1. Requirements #

* Linux operating system
* C++ Compiler
* GNU Make
* [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page)


# 2. Quick Start #
# 2.1 How to build it? #
**Step 1**: Download the library.<br>
Download the library from the github folder.

**Step 2**: Decompress. <br>
Decompress the zip file with your favorite decompressor.

**Step 3**: Install Eigen 3, unzip the documentation and create symbolic links. <br>
SOQCS needs Eigen 3 external library to be built. Eigen 3 can be automatically downloaded and installed using the configuration script <i>config.sh</i> in SOQCS root folder.
This configuration script also unzips the HTML documentation of the library and creates a few symbolic link within the SOQCS source tree needed to build the library. It also creates
a symbolic link to the HTML documentation in the library root folder named <b>Documentation.html</b><br>
<br>
Note that you may need to give execution permission to the configuration script:

```bash
chmod 744 config.sh
```  

**Step 4**: Build the library. <br>
Inside the library main SOQCS folder type <i>make</i>. This will build the library and all the examples. 

# 2.2 How to use it? #
<p align="justify"> Five examples in C++ can be found in the <i>examples</i> subfolder. They can be compiled with the whole library following the instructions above or typing <i>make</i> within the examples folder after compiling the library.
Additionally, there are also versions of those same examples in Python that can be found in the same folder as Jupyter notebooks. The Python port is still in development and it does not contain yet all the functionality of the parent C++ library
therefore, currently not all examples in C++ have their Python counterpart.
</p>

* **Example 1  [C++] [Python]**: Elementary example program to show a basic simulation in SOQCS.
* **Example 2  [C++] [Python]**: An example of HOM visibility using a beamsplitter and physical detectors.
* **Example 3  [C++] [Python]**: An example of the delay gate.
* **Example 4  [C++]**: A simulation of the entanglement swapping protocol. Example of use of density matrices in SOQCS.
* **Example 5  [C++] [Python]**: A boson sampling example.

For extended information about how to use SOQCS library in your own projects check the documentation.
# 3. Documentation #
For more details about how to program with the SOQCS library, add it to your project or learn the details of the available methods and classes consult the available documentation that can be found in the root folder of SOQCS. It can be accessed by clicking on **Documentation.html** after downloading the library and executing the configuration script. Otherwise, it will be found compressed in the file doc.zip .

# 4. Authorship #
<b>Javier Osca</b> <br>
javier.oscacotarelo@mu.ie

<b>Jiri Vala</b> <br>
jiri.vala@mu.ie

# 5. License and Copyright #
Copyright (c) 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved. This library and its related files are subject to the licence terms detailed in LICENCE.TXT .
Use of SOQCS is only permitted under the terms of the licence in [LICENCE.TXT](./LICENCE.TXT). 
