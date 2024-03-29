##########################################################################################################
#                                                                                                        #
#  SOQCS interface with Python. Python side.                                                             #
#                                                                                                        #
#                                                                                                        #
# Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.    #
# The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    #
# root directory of this source tree. Use of the source code in this file is only permitted under the    #
# terms of the licence in LICENCE.TXT.                                                                   #
##########################################################################################################

import os
import math 
import cmath 
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from ctypes import cdll,c_char,c_int,c_long,c_double,POINTER

#------------------------------------------------------------------------------#      
# CPP library configuration
#------------------------------------------------------------------------------#      
cpplib='libpysoqcs.so'

#------------------------------------------------------------------------------#      
# Open SOQCS C++ interface                                                   #
#------------------------------------------------------------------------------#      
soqcs = cdll.LoadLibrary(os.path.dirname(__file__) + '/' + cpplib)

#------------------------------------------------------------------------------#     
# Memory management auxiliary methods.                                         #
# Free the memory of an object created in C++                                  #
# Not valid for classes                                                        #
#------------------------------------------------------------------------------#      
def free_ptr(array):
    """

    Free memory assigned to pointers defined in C++
    
    :array(POINTER): Vector pointer
    
    """    
    array.restype=POINTER(c_char)
    soqcs.free_ptr(array)

#------------------------------------------------------------------------------#      
# Converter MTX->int*                                                          #
#------------------------------------------------------------------------------#      
def to_int_ptr(mtx):
    """

    Converts a mtx into a integer pointer.
    
    :mtx(list[][]): Matrix
    
    """  
    n=len(mtx)
    m=len(mtx[0])
    length=n*m


    send=(c_int*length)()
    for i in range(0,n):
        for j in range(0,m):            
            send[i*m+j]=mtx[i][j]

    return send,n,m

#------------------------------------------------------------------------------#      
# Converter VEC->int*                                                          #
#------------------------------------------------------------------------------#      
def to_int_vec(vec):
    """

    Converts a list into a integer pointer.
    
    :vec(list[]): List
    
    """  

    n=len(vec)

    send=(c_int*n)()
    for i in range(0,n):
        send[i]=vec[i]

    return send,n

#------------------------------------------------------------------------------#
# Wrapper for C++ SOQCS configuration method                                   #
# Configuration of the maximum number of photons of the simulation             #
#------------------------------------------------------------------------------#     
def cfg_soqcs(nph):
    """

    It defines the maximum number of photons to be used by SOQCS.
    
    :nph (int): Maximum number of photons used by SOQCS.
    
    """  
    soqcs.all_cfg_soqcs(nph);   

    
#------------------------------------------------------------------------------#      
# Wrapper for C++ SOQCS class qocircit                                         #
# Definition of a quantum optical circuit                                      #
#------------------------------------------------------------------------------#     
class qocircuit(object):
    """

    A quantum circuit consist in a set of interconnected quantum optical elements.

    :nch (int): Number of channels
    :nm (optional[int]): Number of modes
    :ns (optional[int]): Number of packets
    :np (optional[int]): Number of periods
    :dtp (optional[double]): Length of the periods
    :clock (optional[int]):  Detector behavior configuration: |br|
                            0: Counter. The detectors behave as counters. |br|
                            1: Time. The detectors are able distinguish arrival times but they are blind to frequency. (This mode may change packet numeration).|br|
                            2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. |br|
                            3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. |br|
                            4: Period.  Detectors can only distinguish the period in which photons arrive. |br|
    :R (optional[int]): Number of iterations in the calculation of detector dead time and dark counts effects.
    :loss (optional[bool]): Explicit computation of losses.
    :ckind (optional[char]): Packet shape: |br|
                             'G': Gaussian |br|
                             'E': Exonetial |br|
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """

    #---------------------------------------------------------------------------      
    # Create circuit
    #---------------------------------------------------------------------------     
    def __init__(self, nch, nm=1, ns=1, np=1, dtp=-1.0, clock=0, R=0, loss=False, ckind='G', dummy=False):
        if(ckind=='G'): 
            ikind=0
        else:
            ikind=1
            
        func=soqcs.qoc_new_qocircuit
        func.restype=c_long    
        if(dummy==False):
            self.obj = func(nch,nm,ns,np,c_double(dtp),clock,R,loss,ikind)

    #---------------------------------------------------------------------------              
    # Delete circuit        
    #---------------------------------------------------------------------------      
    def __del__(self):       
        soqcs.qoc_destroy_qocircuit(c_long(self.obj))

    #---------------------------------------------------------------------------              
    # Return nunmber of levels          
    #---------------------------------------------------------------------------      
    def num_levels(self):
        """

        Returns the total number of levels of the circuit..

        :return(int): Number of levels of the circuit.
        """
        return soqcs.qoc_num_levels(c_long(self.obj))

    #---------------------------------------------------------------------------          
    # Adds to the circuit a random circuit
    #---------------------------------------------------------------------------      
    def random_circuit(self):
        """
        
        Creates a circuit defined by a random unitary matrix.
        
        
        """
        soqcs.qoc_random_circuit(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a NSX subcircuit
    #---------------------------------------------------------------------------      
    def NSX(self, ch1, ch2, ch3):
        """
            
        Adds a NSX circuit element. Post-selection still has to be carried out to obtain the proper functionality.
        
        :ch1 (int): NSX input channel 1.
        :ch2 (int): NSX input channel 2.
        :ch3 (int): NSX input channel 3.
        
        
        """
        soqcs.qoc_NSX(c_long(self.obj),ch1,ch2,ch3)

    #---------------------------------------------------------------------------      
    # Adds to the circuit an ideal beamsplitter
    #---------------------------------------------------------------------------      
    def beamsplitter(self, ch1, ch2, theta, phi):
        """
            
        Adds a beamsplitter to the quantum circuit attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        
        """
        soqcs.qoc_beamsplitter(c_long(self.obj),ch1,ch2,c_double(theta),c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a dielectric beamsplitter
    #---------------------------------------------------------------------------      
    def dielectric(self, ch1, ch2, t, r):
        """
            
        Adds a dieletric film attached to channels ch1 and ch2.
        It may also work as a dieletric beamsplitter.
        
        :ch1 (int): Dielectric input channel 1.
        :ch2 (int): Dielectric input channel 2.
        
        """
        soqcs.qoc_dielectric(c_long(self.obj),ch1,ch2,c_double(t.real),c_double(t.imag),c_double(r.real),c_double(r.imag))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a MMI 2x2 beamsplitter
    #---------------------------------------------------------------------------      
    def MMI2(self, ch1, ch2):
        """
            
        Adds a 2x2 MMI to the circut attached to channels ch1 and ch2.
        
        :ch1 (int): MMI2 input channel 1.
        :ch2 (int): MMI2 input channel 2.
        
        """
        soqcs.qoc_MMI2(c_long(self.obj),ch1,ch2)

    #---------------------------------------------------------------------------      
    # Adds to the circuit a rewire gate
    #---------------------------------------------------------------------------      
    def rewire(self, ch1, ch2):
        """
            
        Adds a swap gate between two channels
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        
        """
        soqcs.qoc_rewire(c_long(self.obj),ch1,ch2)

    #---------------------------------------------------------------------------          
    # Adds to the circuit a general homogeneous medium
    #---------------------------------------------------------------------------      
    def general_medium(self, ch, t):
        """
            
        Adds a general homogeneous medium to the circuit
        
        :ch (int): Input channel.
        :t (complex): Transmitance
        
        """
        soqcs.qoc_phase_shifter(c_long(self.obj),ch,c_double(t.real),c_double(t.imag))
        
    #---------------------------------------------------------------------------              
    # Adds to the circuit a phase shifter
    #---------------------------------------------------------------------------      
    def phase_shifter(self, ch, phi):
        """
            
        Adds a phase shifter to the circuit in channel ch.
        
        :ch (int): Phase shifter input channel.
        :phi (float): Angle phi in degrees.
        
        """
        phir=math.pi*phi/180
        t=cmath.exp(1j*phir)
        self.general_medium(ch,t)

    #---------------------------------------------------------------------------          
    # Adds to the circuit a lossy medium
    #---------------------------------------------------------------------------      
    def loss(self, ch, l):
        """
            
        Adds a lossy medium with loss probability l to the circuit in channel ch.
        
        :ch (int): Lossy medium input channel.
        :l  (float): Loss probability.
        
        """
        self.general_medium(ch, math.sqrt(1.0-l))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a rotator
    #---------------------------------------------------------------------------      
    def rotator(self, ch, theta, phi):
        """
        
        Adds a polarization rotation element to the circuit attached to channel ch.
        
        :ch (int):  Rotator input channel.
        :theta (double):  Angle theta in degrees.
        :phi   (double):  Angle phi in degrees.
                
        """
        soqcs.qoc_rotator(c_long(self.obj),ch,c_double(theta),c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a polarized beamsplitter
    #---------------------------------------------------------------------------      
    def pol_beamsplitter(self, ch1, ch2, P, theta):
        """
            
        Adds a polarized beamsplitter attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        :P   (int): Polarization to which the beamsplitter is sensitive.
        :theta (float): Effectiveness of the beamsplitter. 90º => 50/50 beamsplitter. 0º => No sensitivity.
        
        """
        soqcs.qoc_pol_beamsplitter(c_long(self.obj),ch1, ch2, P,c_double(theta))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a polarized phase shifter
    #---------------------------------------------------------------------------      
    def pol_phase_shifter(self, ch, P, phi):
        """
            
        Adds a polarized phase shifter attached to channel ch.
        
        :ch  (int): Phase shifter input channel.
        :P   (int): Polarization to which the phase shifter is sensitive.
        :phi (double):  Phase in degrees.
        
        """
        soqcs.qoc_pol_phase_shifter(c_long(self.obj),ch, P,c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a polarization filter
    #---------------------------------------------------------------------------  
    def pol_filter(self, ch, P):
        """
            
        Adds a polarization filter attached to channel ch.
        
        :ch  (int): Polarization filter input channel.
        :P   (int): Polarization to be filtered.
        
        """
        soqcs.qoc_pol_filter(c_long(self.obj),ch, P)
        
    #---------------------------------------------------------------------------      
    # Adds to the circuit a waveplate
    #---------------------------------------------------------------------------      
    def waveplate(self, ch, alpha, gamma):
        """
            
        Adds a general waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double): Rotation angle in degrees.
        :gamma (double): Wavenumber delay of one component with respect the other in degrees.
        
        """
        soqcs.qoc_waveplate(c_long(self.obj),ch, c_double(alpha), c_double(gamma))

    #---------------------------------------------------------------------------      
    # Adds to the circuit a half waveplate
    #---------------------------------------------------------------------------      
    def half(self, ch, alpha):  
        """
            
        Adds a half waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double): Rotation angle in degrees.
        
        """
        self.waveplate(ch, alpha, 90.0)

    #---------------------------------------------------------------------------      
    # Adds to the circuit a quarter waveplate
    #---------------------------------------------------------------------------      
    def quarter(self, ch, alpha):   
        """
            
        Adds a quarter waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double):  Rotation angle in degrees.
        
        """
        self.waveplate(ch, alpha, 45.0)

    #---------------------------------------------------------------------------      
    # Adds to the circuit a detector
    #---------------------------------------------------------------------------      
    def detector(self, ch, cond=-1, pol=-1, mpi=-1, mpo=-1, eff=1.0, blnk=0.0, gamma=0.0):
        """
            
        Adds a detector to a channel of the circuit.
        
        :ch   (int): Detector channel.
        :cond (optional[int]): Detection condition. |br|
         cond>=0:  Readings in the remaining channels are considered only by calculations in probability bins and density matrices if the number of photons in this channel are equal to cond. |br|
         cond=-1:  There is no condition and works as a normal detector. |br|
         cond=-2:  The channel is ignored by outcome calculations in probability bins and density matrices. |br|
        :pol (optional[int]): Polarization condition. If cond>=0, pol determines the polarization of the photons to fulfill the condition. Note that if pol=-1 no assumption about the polarization of those photons is made.
        :mpi (optional[int]): Initial period of the detection window (if -1 takes the first one as default).
        :mpo (optional[int]): Final period of the detection window (if -1 takes the first one as default).
        :eff   (optional[float]): Efficiency of the detector.
        :blnk  (optional[float]): Ratio of time in which the detector is inactive due other detections.
        :gamma (optional[float]): Average rate of dark counts in this channel.
        
        """
        soqcs.qoc_detector(c_long(self.obj),c_int(ch),c_int(cond),c_int(pol),c_int(mpi),c_int(mpo),c_double(eff),c_double(blnk),c_double(gamma))

    #---------------------------------------------------------------------------      
    # Ignore detections on a channel
    #---------------------------------------------------------------------------      
    def ignore(self, ch):
        """
            
        Flags the channel to be ignored in the outcome calculations. No detector is placed in this channel.
        
        :ch (int): Ignored channel.
        
        """
        soqcs.qoc_detector(c_long(self.obj),c_int(ch),c_int(-2),c_int(-1),c_int(-1),c_int(-1),c_double(1.0),c_double(0.0),c_double(0.0))

    #---------------------------------------------------------------------------              
    # Adds to the output a random Gaussian whote noise of stdev^2 dispersion
    #---------------------------------------------------------------------------      
    def noise(self, stdev2):
        """
            
        Adds Gaussian white noise to the output.
        
        :stdev2 (float): Dispersion of the Gaussian noise.
        
        """
        soqcs.qoc_noise(c_long(self.obj),c_double(stdev2))

    #---------------------------------------------------------------------------              
    # Adds to the circuit a packet definition
    #---------------------------------------------------------------------------      
    def def_packet(self, n,t,f,w): 
        """
            
        Adds to the circuit a new definition of a wave packet. |br|
        **Warning!**  Note that when the detectors are not configured as counters the packet numeration may change on emission.

        
        :n(int):   Suggested wavepacket number.
        :t(float): Wavepacket characteristic emission time.
        :f(float): Wavepacket characteristic frequency.
        :w(float): Width or decay length depending if the packet shape model is Gaussian or Exponential.                
        
        """
        func=soqcs.qoc_def_packet
        func.restype=c_int
        return func(c_long(self.obj),n,c_double(t),c_double(f),c_double(w))
      

    #---------------------------------------------------------------------------      
    # Calculates the emitted packet visibility
    #---------------------------------------------------------------------------      
    def emitted_vis(self, i,j):
        """
        
        Overlap probability between two wave packets.
        
        :i (int):  Packet i. 
        :j (int):  Packet j. 
        :return(float): Probability of overlap.
        
        """
        func=soqcs.qoc_emitted_vis
        func.restype=c_double
        return func(c_long(self.obj),i,j)                

    #---------------------------------------------------------------------------      
    # Adds to the circuit an emitter
    #---------------------------------------------------------------------------      
    def emitter(self):
        """
        
        Adds an emitter to the circuit using the packet definition given by def_packet.
        
        
        """
        soqcs.qoc_emitter(c_long(self.obj))     

    #---------------------------------------------------------------------------          
    #  Adds to the circuit a delay
    #---------------------------------------------------------------------------      
    def delay(self, ch):
        """
        
        Increases the optical path of a channel by a quantity equal to a period.
        
        :ch (int):  Channel where the delay is introduced.
                
        """
        soqcs.qoc_delay(c_long(self.obj),ch)     

    #---------------------------------------------------------------------------      
    # Print circuit
    #---------------------------------------------------------------------------      
    def prnt(self, fmt=0):  
        """
        
        Prints the circuit matrix.
        
       :fmt(int): Flag that controls the print style. |br|
        0 = Prints numerically. |br|
        N = Prints the polarization using the alphabet (H/V). |br|
                
        """
        soqcs.qoc_prnt(c_long(self.obj),fmt)
        sys.stdout.flush()
        

#------------------------------------------------------------------------------#      
# Wrapper for C++ SOQCS class state                                            #
# Definition of a quantum photonic state                                       #
#------------------------------------------------------------------------------#      
class state(object):
    """

    Quantum bosonic state definition. A quantum state is a sum of kets multiplied by amplitudes.

    :level (int): Number of levels to describe the state.
    :nph(optional(int)): Maximum number of photons.
    :maxket(optional(int)): Maximum number of different terms in the summation. (Internal memory).
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """   
    #---------------------------------------------------------------------------      
    # Create a state
    #---------------------------------------------------------------------------      
    def __init__(self, level, nph=-1, maxket=50,dummy=False):
        func=soqcs.st_new_state
        func.restype=c_long
        if(dummy==False):
            self.obj = func(nph, level,maxket)

    #---------------------------------------------------------------------------              
    # Delete a state
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.st_destroy_state(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Calculate braket coefficient between two states
    #---------------------------------------------------------------------------      
    def braket(self,state):   
        """

        Performs the bra-ket operation <bra|state> using the complex conjugate of this state as bra.
        
        :state(state): State in the right hand side of the braket operation.         
        :return (complex): The complex number result of the projection.
            
        """  
        func=soqcs.st_braket
        func.restype=POINTER(c_double)
        array_ptr=2*c_double
        array_ptr=func(c_long(self.obj),c_long(state.obj))
        value=complex(c_double(array_ptr[0]).value,c_double(array_ptr[1]).value)        
        free_ptr(array_ptr)
        return value

    #---------------------------------------------------------------------------      
    # Normalizes a state
    #---------------------------------------------------------------------------      
    def normalize(self):
        """

        Normalizes the state.         
        
        """  
        soqcs.st_normalize(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Rephase
    #---------------------------------------------------------------------------      
    def rephase(self,mtx,qoc):
        """

        Changes the global phase of the state to make real the coefficient of the reference ket defined by mtx.
         
        :mtx(list[][]):  Matrix that defines the ket which coefficient is used as reference. Each column defines the configuration of one level. |br|
                         There are four different ways to create these configurations depending on the number of rows
                         of the term matrix: |br|
                         |br|
                         4-Row: Channels, polarization, wavepacket and occupation in this order. |br|
                         3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases. |br|
                         2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases. |br|
                         1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering. |br|
                         |br|
                         Except in the 1-Row case the order is irrelevant. Furthermore, levels not configured are assumed to be initialized by zero.
                 
        :qoc(qocircuit): Circuit to which the term is related.
        
        """  
        param=to_int_ptr(mtx)
        soqcs.st_rephase(c_long(self.obj), param[0],param[1],param[2],c_long(qoc.obj)) 

    #---------------------------------------------------------------------------              
    # Adds a new term (ket + amplitude) to a state using a state description
    #---------------------------------------------------------------------------      
    def add_term(self, ampl, mtx,qoc):
        """

        Adds a new term to the state using a definition of the term.
         
        :ampl(complex):  Amplitude of the new term.
        :mtx(list[][]):  Matrix that defines the new term. Each column defines the configuration of one level. |br|
                         There are four different ways to create these configurations depending on the number of rows
                         of the term matrix: |br|
                         |br|
                         4-Row: Channels, polarization, wavepacket and occupation in this order. |br|
                         3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases. |br|
                         2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases. |br|
                         1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering. |br|
                         |br|
                         Except in the 1-Row case the order is irrelevant. Furthermore, levels not configured are assumed to be initialized by zero.
                 
        :qoc(qocircuit): Circuit to which the term is related.
        
        """    
        param=to_int_ptr(mtx)
        soqcs.st_add_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), param[0],param[1],param[2],c_long(qoc.obj)) 

    #---------------------------------------------------------------------------      
    #  Adds a new term (ket + amplitude) to a state directly from level definition
    #---------------------------------------------------------------------------      
    def add_ket(self, ampl, vec):   
        """

        Adds a new term to the state using a vector with the level occupation.
         
        :ampl(complex):  Amplitude of the new term.
        :vel(list[int]): List with the occupation of each level in the new term.
        :qoc(qocircuit): Circuit to which the term is related.
        
        """  
        length=len(vec)
        send=(c_int*length)()
        for i in range(0,length):
            send[i]=vec[i]
            
        soqcs.st_add_raw_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), send) 

    #---------------------------------------------------------------------------      
    # Post-selection by a projector
    #---------------------------------------------------------------------------      
    def post_selection(self, prj):   
        """

        Post-selection over states by the condition defined in the the "projector".
         
        :prj(projector): Projector with the description of the levels (and their occupations) to be post-selected.
        
        """            
        func=soqcs.st_post_selection
        func.restype=c_long
        obj=func(c_long(self.obj), c_long(prj.obj))
        newstate=state(1,True) 
        newstate.obj=obj
        return newstate

    #---------------------------------------------------------------------------      
    # Encode a state into qubit encoding. (Path encodig version)
    # Those kets that can not be encoded are ignored and the result is normalized
    # WARNING: It is responsability of the library user to make sure that not 
    # repetitions arise after encoding.
    #---------------------------------------------------------------------------      
    def encode(self, qmap,  qoc):
        """
        Encode a state into a qubit encoding (Path encoding version). Those kets that can not be encoded are dismissed and the result is normalized. |br|
        **Only for ideal circuits (nm=1 and ns=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.


        :qmap(list[][]): 2xn matrix with the qubit definitions. Each column has two entries
                         specifying the channels that define the qubit.
        :qoc(qocircuit): Circuit which the state is related.
        :return output(state): An encoded state.
        
        """  
        param=to_int_ptr(qmap)
        func=soqcs.st_encode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[2],c_long(qoc.obj))  
        aux=state(1,True)
        aux.obj=obj
        return aux 

   #----------------------------------------------------------------------------      
   # Decode a state from a qubit encoding into photon modes. (Path encoding version)
   #----------------------------------------------------------------------------
    def decode(self, qmap, ancilla, qoc):
        """
        Decode a state from a qubit encoding into photon encoding. (Path encoding version). |br|
        **Only for ideal circuits (nm=1 and ns=1)**. |br|        
        Note that we need to define the values of the extra ancilla channels that are not included into the qubit encoding
        but they are needed to the circuit to work in the intended way


        :qmap(list[][]): 2xn matrix with the qubit definitions. Each column has two entries
                         specifying the channels that define the qubit.
        :ancilla(list[]): List with the values of the ancilla channels (from smaller to large channel number).
        :qoc(qocircuit): Circuit which the state is related.
        :return output(state): A decoded state.
        
        """  
        
        nq=len(qmap[0])
        nanz=len(ancilla)
        nch=2*nq+nanz   

      
        # Compute ancillas
        aux=[1]*nch
        for i in range(0,nq):
            aux[qmap[0][i]]=0
            aux[qmap[1][i]]=0
            
        k=0
        occ=[]
        for i in range(0,nch):
            if aux[i]==0:
                occ.append(0)
            if aux[i]==1:
                occ.append(ancilla[k])
                k=k+1
                
        anzstate= state(nch,maxket=1)
        occxt=[occ]
        anzstate.add_term(1.0, occxt,qoc)
        
        
        param=to_int_ptr(qmap)
        func=soqcs.st_decode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[2],c_long(anzstate.obj),c_long(qoc.obj))  
        aux=state(1,True)
        aux.obj=obj
        return aux 

    #---------------------------------------------------------------------------      
    # Encode a state into qubit encoding. (Polarization encodig version)
    # Those kets that can not be encoded are ignored and the result is normalized
    # WARNING: It is responsability of the library user to make sure that not 
    # repetitions arise after encoding.
    #---------------------------------------------------------------------------      
    def pol_encode(self, qvec,  qoc):
        """
        Encode a state into a qubit encoding (Polarization encoding version). Those kets that can not be encoded are dismissed and the result is normalized. |br|
        **Only for ideal circuits (ns=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.


        :qvec(list[]): List with the qubit definitions. 
        :qoc(qocircuit): Circuit which the state is related.
        :return output(state): An encoded state.
        
        """  
        param=to_int_vec(qvec)
        func=soqcs.st_pol_encode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(qoc.obj))  
        aux=state(1,True)
        aux.obj=obj
        return aux 

   #----------------------------------------------------------------------------      
   # Decode a state from a qubit encoding into photon modes. (Polarization encoding version)
   #----------------------------------------------------------------------------
    def pol_decode(self, qvec, ancilla, qoc):
        """
        Decode a state from a qubit encoding into photon encoding. (Polarization encoding version).|br|
        **Only for ideal circuits (ns=1)**. |br|        
        Note that we need to define the values of the extra ancilla channels that are not included into the qubit encoding
        but they are needed to the circuit to work in the intended way


        :qvec(list[]): List with the qubit definitions. 
        :ancilla(list[]): List with the values of the ancilla channels (from smaller to large channel number). 0="Horizontal"/1="Vertical". 
        :qoc(qocircuit): Circuit which the state is related.
        :return output(state): A decoded state.
        
        """  
        
        nq=len(qvec)
        nanz=len(ancilla)
        nch=2*nq+nanz   

      
        # Compute ancillas
        aux=[1]*nch
        for i in range(0,nq):
            aux[qvec[i]]=0

            
        k=0
        l=0
        occ=np.zeros((3,nch*2))
        for i in range(0,nch):
            for j in range(0,2):            
                if aux[i]==0:
                    occ[0][l]=i
                    occ[1][l]=j
                    occ[1][l]=0
                    
                if aux[i]==1:
                    occ[0][l]=i
                    occ[1][l]=j
                    occ[1][l]=ancilla[j][k]
                    k=k+1
            l=l+1
            
        anzstate= state(nch,maxket=1)
        anzstate.add_term(1.0, occ,qoc)
        
        
        param=to_int_vec(qvec)
        func=soqcs.st_pol_decode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(anzstate.obj),c_long(qoc.obj))  
        aux=state(1,True)
        aux.obj=obj
        return aux 
    
    #---------------------------------------------------------------------------      
    # Print state
    #---------------------------------------------------------------------------      
    def prnt_state(self,format=0,column=0,loss=False,qoc=None):
        """
            Prints a state in one of the specified the formats. If loss=True loss channels are printed in blue.
            
            :format (optional[int]): Format of the output.
            :column (optional[int]): Print in row or column. 0=Rows/1=Columns
            :loss (optional[bool]): Print loss channels in different color (True=Yes, False=No)
            :qoc (optional[qocircuit]): Quantum circuit to whom the state is referred.
        """  
        if(qoc==None):
            aux=0
        else:
            qoc.restype=qocircuit
            aux=qoc.obj
                                    
        soqcs.st_prnt_state(c_long(self.obj), format, column, loss, c_long(aux))
        sys.stdout.flush()
                                    

#---------------------------------------------------------------------------                  
# Wrapper for C++ SOQCS class projector
# Definition of a projector to perform post-selection
#---------------------------------------------------------------------------      
class projector(object):
    """

    Condition that have to be met to accept a state in post-selection.

    :level (int): Number of levels
    :nph(optional(int)): Maximum number of photons.
    :maxket(optional(int)): Maximum number of terms to define the projector. (Internal memory)
    
    """  
    #---------------------------------------------------------------------------      
    # Create a projector
    #---------------------------------------------------------------------------      
    def __init__(self,nlevel,nph=-1,maxket=10):
        func=soqcs.prj_new_projector
        func.restype=c_long    
        self.obj = func(nph,nlevel,maxket)

    #---------------------------------------------------------------------------          
    # Delete a projector
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.prj_destroy_projector(c_long(self.obj))

    #---------------------------------------------------------------------------                
    # Add a new term (ket + amplitude) to a projector    
    #---------------------------------------------------------------------------      
    def add_term(self, ampl, mtx,qoc):
        """

        Adds a new term to the projector.
         
        :ampl(complex):  Amplitude of the new term.
        :mtx(list[][]):  Matrix that defines the new tern. Each column defines the configuration of one level.|br|
                         There are four different ways to create these configurations depending on the number of rows
                         of the term matrix: |br|
                         |br|
                         4-Row: Channels, polarization, wavepacket and occupation in this order. |br|
                         3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases. |br|
                         2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases. |br|
                         1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering. |br|
                         |br|
                         Except in the 1-Row case the order is irrelevant.
                 
        :qoc(qocircuit): Circuit to which the term is related.
        
        """  
        param=to_int_ptr(mtx)
        soqcs.prj_add_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), param[0],param[1],param[2],c_long(qoc.obj)) 


#------------------------------------------------------------------------------#
# Wrapper for C++ SOQCS class p_bin                                            #
# Definition of a set of probability bins                                      #
#------------------------------------------------------------------------------#      
class p_bin(object):
    """

    Outcomes. Set of bins where in each bin it is stored the probability of one particular outcome.
    
    :level (int): Number of levels to describe the bins.
    :nph(optional(int)): Maximum number of photons.    
    :maxket(optional(int)): Maximum number of bins possible. (Internal memory)
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """  
    #---------------------------------------------------------------------------      
    # Create a p_bin
    #---------------------------------------------------------------------------      
    def __init__(self, level, nph=-1, maxket=50,dummy=False):
        func=soqcs.pb_new_pbin
        func.restype=c_long
        if(dummy==False):
            self.obj = func(nph,level,maxket)

    #---------------------------------------------------------------------------                      
    # Delete a p_bin
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.pb_destroy_pbin(c_long(self.obj))

    #---------------------------------------------------------------------------                
    # Add the statistics of a state to a pbin
    #---------------------------------------------------------------------------      
    def add_state(self, state):
        """

        Adds the statistics of a quantum state.
    
        :state (state): Input state.
        
        """
        soqcs.pb_add_state(c_long(self.obj),c_long(state.obj)) 

    #---------------------------------------------------------------------------      
    # Total probability of a set of bins. Maybe inferior to one because of the
    # post-selection
    #---------------------------------------------------------------------------              
    def trace(self):
        """
        
        Obtains the total probability of the set of bins. It may be inferior to one due post-selection or encoding procedures.
        
        :return(double): Normalization value of the set of probability bins.
    
        """
        func=soqcs.pb_trace
        func.restype=c_double
        return func(c_long(self.obj))

    #---------------------------------------------------------------------------              
    # Normalize the p_pin. All the probabilities sum to one. 
    # This may not be true if post-selection has been applied.
    #---------------------------------------------------------------------------      s
    def normalize(self):
        """
        
        Normalizes the total probability of the set of bins to one.
    
        """
        soqcs.pb_normalize(c_long(self.obj)) 

    #---------------------------------------------------------------------------      
    # Calculate the measurement
    #---------------------------------------------------------------------------      
    def calc_measure(self, qoc):
        """
        
        Calculates the effect of the detectors defined in a circuit over the outcome contained in this set of probability bins.
        These effect are post-selection-conditions, detection window ,dark counts, detector dead time, losses and circuit noise.
    
        :qoc (qocircuit):  Circuit where the detectors are defined.
        :return (p_bin): Outcomes after detector effects calculations.
            
        """
        func=soqcs.pb_calc_measure
        func.restype=c_long    
        obj=func(c_long(self.obj),c_long(qoc.obj)) 
        aux=p_bin(1,True)
        aux.obj=obj
        return aux

    #---------------------------------------------------------------------------          
    # Return the number of bins
    #---------------------------------------------------------------------------      
    def nbins(self):
        """
           Returns the number of bins stored.        
    
            :return (int): Number of bins.
        """
        return soqcs.pb_nbins(c_long(self.obj))

    #---------------------------------------------------------------------------          
    # Return nunmber of levels          
    #---------------------------------------------------------------------------      
    def num_levels(self):
        """
           Returns the number of levels that define each bin.
    
            :return (int): Number of levels.
        """
        return soqcs.pb_num_levels(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Return tag of the bin ("bin name")
    #---------------------------------------------------------------------------      
    def tag(self, index):
        """
        
        Returns the occupation of the bin referred by the index in a string format
    
        :index (int):  Index of a bin.
        :return (int): Occupation in string format.
        
        """
        nlevel=self.num_levels()
        func=soqcs.pb_tag
        func.restype=POINTER(c_char)
        array_ptr=c_char*nlevel
        array_ptr=func(c_long(self.obj),index)
        array=[0]*nlevel
        for i in range(0, nlevel):
            array[i]=array_ptr[i].decode('UTF-8')
        free_ptr(array_ptr)
        return "".join(array)

    #---------------------------------------------------------------------------      
    # Return probability of a bin
    #---------------------------------------------------------------------------      
    def prob_idx(self, index):
        """
        
        Returns the probability of the outcome stored in the bin referred by the index.
    
        :index (int):  Index of a bin
        :return (int): Returns the probability of the referred bin.
        
        """
        func=soqcs.pb_prob
        func.restype=c_double
        return func(c_long(self.obj),index)

    #---------------------------------------------------------------------------          
    # Return probability of a bin from its definition (using a circuit)
    #---------------------------------------------------------------------------      
    def prob_qoc(self, mtx,qoc):
        """
        Returns the probability of the outcome stored in the bin described by the provided definition (Circuit version). |br|

        :mtx (lilst[][]): Matrix that defines the bin. Each column defines the configuration of one level. |br|
         There are four different ways to create a definition depending on the number of rows
         of the matrix: |br|
         4-Row: Channels, polarization, wavepacket and occupation in this order. |br|
         3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases. |br|
         2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases. |br|
         1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering. |br|
        :qoc (qocircuit): Circuit to which the outcomes are referred.
        :return (int): Returns the probability of the referred bin.
           
        """
        func=soqcs.pb_prob_def_qoc
        func.restype=c_double
        param=to_int_ptr(mtx)
        return func(c_long(self.obj), param[0],param[1],param[2],c_long(qoc.obj))    

    #---------------------------------------------------------------------------      
    # Return probability of a bin from its definition 
    # (using a generalized metacircuit)
    #---------------------------------------------------------------------------      
    def prob(self, mtx,dev):
        """
        Returns the probability of the outcome stored in the bin described by the provided definition (Device version). |br|

        :mtx (lilst[][]): Matrix that defines the bin. Each column defines the configuration of one level. |br|
         There are four different ways to create a definition depending on the number of rows
         of the matrix: |br|
         4-Row: Channels, polarization, wavepacket and occupation in this order. |br|
         3-Row: Channels, polarization and occupation in this order. Wavepackets are assumed to be the same in all cases. |br|
         2-Row: Channels and occupation in this order. Polarization and wavepackets are assumed to be the same in all cases. |br|
         1-Row: Occupation. In this case the user must provide just a list of occupations for each channel in level ordering. |br|
        :dev (qodevice): Device to which the outcomes are referred.
        :return (int): Returns the probability of the referred bin.
           
        """
        func=soqcs.pb_prob_def
        func.restype=c_double
        param=to_int_ptr(mtx)
        return func(c_long(self.obj), param[0],param[1],param[2],c_long(dev.circ.obj))   

    #---------------------------------------------------------------------------      
    # Translate the labels of a probability bins into qubit encoding. (Path encoding version).
    # Those that can not be encoded are ignored and the result is normalized.
    #---------------------------------------------------------------------------
    def translate(self, qmap,   dev):
        """
        Translates the labels of a set of probability bins into a qubit encoding. (Path encoding version).
        Those that can not be encoded are ignored and the result is normalized. |br|
        **Only for ideal circuits (nm=1 and nd=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.

        :qmap(list[][]): 2xn matrix with the qubit definitions. Each column has two entries
                         specifying the channels that define the qubit.
        :dev(qodev): Device to which the set of probability bins is related.
        :return out(p_bin): A translated set of probability bins.
        
        """          
        param=to_int_ptr(qmap)
        func=soqcs.pb_translate
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[2],c_long(dev.circ.obj))  
        aux=p_bin(1,True)
        aux.obj=obj
        return aux

    #---------------------------------------------------------------------------      
    # Translate the labels of a probability bins into qubit encoding.  (Polarization encoding version).
    # Those that can not be encoded are ignored and the result is normalized.
    #---------------------------------------------------------------------------    
    def pol_translate(self, qvec,   dev):
        """
        Translates the labels of a set of probability bins into a qubit encoding. (Polarization encoding version).
        Those that can not be encoded are ignored and the result is normalized. |br|
        **Only for ideal circuits (nm=1 and nd=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.

        :qvec(list[]): List with the qubit definitions. 
        :dev(qodev): Device to which the set of probability bins is related.
        :return out(p_bin): A translated set of probability bins.
        
        """          
        param=to_int_vec(qvec)
        func=soqcs.pb_pol_translate
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(dev.circ.obj))  
        aux=p_bin(1,True)
        aux.obj=obj
        return aux
    
    #---------------------------------------------------------------------------          
    # Print a set of probability bins (using a circuit)
    #---------------------------------------------------------------------------      
    def prnt_bins_qoc(self,format=0,thresh=0.0,loss=False,qoc=None):
        """
            Prints a set of probability bins (Circuit version).
            
            :format (optional[int]): Format of the output.
            :thresh (optional[float]): Minimum value of a bin to be printed.
            :loss (optional[bool]): Print loss channels in different color (True=Yes, False=No)
            :qoc (optional[qocircuit]): Quantum circuit to whom the bins are referred.
        """     
        if(qoc==None):
            aux=0
        else:
            qoc.restype=qocircuit
            aux=qoc.obj
                        
        soqcs.pb_prnt_bins_qoc(c_long(self.obj), format, c_double(thresh), loss, c_long(aux))
        sys.stdout.flush()            

    #---------------------------------------------------------------------------      
    # Print a set of probability bins (using a generalized metacircuit)
    #---------------------------------------------------------------------------      
    def prnt_bins(self,format=0,thresh=0.0,loss=False,dev=None):
        """
            Prints a set of probability bins (Device version).
            
            :format (optional[int]): Format of the output.
            :thresh (optional[float]): Minimum value of a bin to be printed.
            :loss (optional[bool]): Print loss channels in different color (True=Yes, False=No)
            :dev (optional[qodev]): Quantum circuit to whom the bins are referred.
        """  
        if(dev==None):
            aux=0
        else:
            dev.restype=qodev
            aux=dev.circ.obj
                                    
        soqcs.pb_prnt_bins(c_long(self.obj), format, c_double(thresh), loss, c_long(aux))            
        sys.stdout.flush()

    #---------------------------------------------------------------------------      
    # Plots a set of probability bins graphically
    #---------------------------------------------------------------------------      
    def show(self, pmax=-1.0,sizex=8, sizey=5, dpi=100, angle=70, font=14):
        """
        
            Plots the outcomes in a bar diagram.
        
            :sizex(optional[float]): Size of the plot in the x direction in inches.
            :sizey(optional[float]): Size of the plot in the y direction in inches.
            :dpi(optional[float]): Density of points  by inch.
            :angle(optional[float]): Angle of the horizontal axis labels.
            :font(optional[fint]): Size of the font of the horizontal axis labels.
            
        """  
        nbins=self.nbins()
        taglist=list()
        problist=list()
        for i in range(0,nbins):
            taglist.append(self.tag(i))
            problist.append(self.prob_idx(i))
  
        zipped_lists = zip(taglist, problist)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        staglist, sproblist = [ list(tuple) for tuple in  tuples]
        
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(sizex,sizey)
        fig.set_dpi(dpi)
        plt.xticks(rotation = angle)
        if pmax>0.0:
            plt.ylim(0,pmax)
        
        y_pos = np.arange(len(staglist))
        pbar=ax.bar(y_pos, sproblist, align='center', alpha=0.5)
        ax.bar_label(pbar,fmt='%.3f')
        ax.tick_params(axis='x', which='major', labelsize=font)
        ax.set_axisbelow(True)
        plt.grid(visible=True, which='major', axis='y', color='0.8', linestyle='--', linewidth=2)
        plt.xticks(y_pos, staglist)
        plt.ylabel('Frequency',fontsize=font)
        plt.title('Output')
       
        plt.draw()
        plt.show()
           

#------------------------------------------------------------------------------#          
# Wrapper for C++ SOQCS class dmatrix                                          #
# Definition of a density matrix                                               #
#------------------------------------------------------------------------------#     
class dmatrix(object):
    """

    Outcomes. Density matrix of states.
    
    :meme(optional(int)): Maximum number of rows(or/and columns) possible. (Internal memory)
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """ 
    #---------------------------------------------------------------------------          
    # Create a density matrix
    #---------------------------------------------------------------------------      
    def __init__(self, mem=1000,dummy=False):
        func=soqcs.dm_new_dmat
        func.restype=c_long
        if(dummy==False):
            self.obj = func(mem) ## Problem here

    #---------------------------------------------------------------------------              
    # Delete a density amtriz
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.dm_destroy_dmat(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Return trace
    #---------------------------------------------------------------------------      
    def trace(self):
        """

        Calculates the trace of the density matrix.
    
        :return(double): Trace of the density matrix,
    
        """ 
        func=soqcs.dm_trace
        func.restype=c_double
        return func(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Normalize the density matrix (diagonal=1)
    #---------------------------------------------------------------------------      
    def normalize(self):
        """

        Normalizes the density matrix to trace = 1.
    
        """ 
        soqcs.dm_normalize(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Calculates the fidelity of the matrix with a referene state
    #---------------------------------------------------------------------------      
    def fidelity(self, state):
        """

        Calculates the fidelity of a density matrix with respect a reference state.
        
        :state(state): Reference state to compare with the density matrix.
        :return(double): Fidelity value.
    
        """ 
        func=soqcs.dm_fidelity
        func.restype=c_double
        return func(c_long(self.obj),c_long(state.obj))

    #---------------------------------------------------------------------------                
    # Add the statistics of a state to a density matrix (using a circuit)
    #---------------------------------------------------------------------------      
    def add_state_qoc(self, state, qoc):
        """

        Adds a new state to the density matrix using the post-selection condition defined by 
        the circuit detectors (Circuit version).
        
        :state(state): State to be added to the set of states used to calculate the density matrix entries.
        :qoc(qocircuit): Circuit where the detectors are defined.        
    
        """ 
        soqcs.dm_add_state_qoc(c_long(self.obj),c_long(state.obj),c_long(qoc.obj)) 

    #---------------------------------------------------------------------------      
    # Calculate the measurement (using a circuit)
    #---------------------------------------------------------------------------      
    def calc_measure_qoc(self, qoc):
        """

        Calculate measure (Circuit version). It means to remove degrees of freedom depending on the 
        circuit configuration of the detectors.
        
        :qoc(qocircuit): Circuit to which the density matrix is related.     
        :return(dmat): Reduced matrix.
    
        """ 
        func=soqcs.dm_calc_measure_qoc
        func.restype=c_long    
        obj=func(c_long(self.obj),c_long(qoc.obj)) 
        aux=dmatrix(1,True)
        aux.obj=obj
        return aux

    #---------------------------------------------------------------------------      
    # Print state (using a circuit)
    #---------------------------------------------------------------------------      
    def prnt_mtx_qoc(self,format=0,thresh=0.0,qoc=None):
        """

        Prints a relevant density submatrix (Circuit version). 
        Only the rows and columns with sum values greater than "thresh" are printed.
        
        :format (optional[int]): Format of the output.
        :thresh (optional[double]): Threshold value to print.
        :qoc (optional[qocircuit]): Circuit to which the density matrix is related.
    
        """ 
        if(qoc==None):
            aux=0
        else:
            qoc.restype=qocircuit
            aux=qoc.obj
                                    
        soqcs.dm_prnt_mtx_qoc(c_long(self.obj), format, c_double(thresh), c_long(aux))
        sys.stdout.flush()

    #---------------------------------------------------------------------------      
    # Add the statistics of a state to a density matrix (using a generalized metacircuit)
    #---------------------------------------------------------------------------      
    def add_state(self, state, dev):
        """

        Adds a new state to the density matrix using the post-selection condition defined by 
        the circuit detectors (Device version).
        
        :state(state): State to be added to the set of states used to calculate the density matrix entries.
        :dev(qodev): Device where the detectors are defined     
    
        """ 
        soqcs.dm_add_state(c_long(self.obj),c_long(state.obj),c_long(dev.circ.obj)) 

    #---------------------------------------------------------------------------      
    # Calculate the measurement (using a generalized metacircuit)
    #---------------------------------------------------------------------------      
    def calc_measure(self, dev):
        """

        Calculate measure (Device version). It means to remove degrees of freedom depending on the 
        circuit configuration of the detectors.
        
        :dev(qodev): Device to which the density matrix is related.
        :return(dmat): Reduced matrix.
    
        """ 
        func=soqcs.dm_calc_measure
        func.restype=c_long    
        obj=func(c_long(self.obj),c_long(dev.circ.obj)) 
        aux=dmatrix(1,True)
        aux.obj=obj
        return aux

    #---------------------------------------------------------------------------      
    # Print state (using a generalized metacircuit)
    #---------------------------------------------------------------------------      
    def prnt_mtx(self,format=0,thresh=0.0,dev=None):
        """

        Prints a relevant density submatrix (Device version). 
        Only the rows and columns with sum values greater than "thresh" are printed.
        
        :format (optional[int]): Format of the output.
        :thresh (optional[double]): Threshold value to print.
        :dev (optional[qodev]): Device to which the density matrix is related.
    
        """ 
        if(dev==None):
            aux=0
        else:
            dev.circ.restype=qocircuit
            aux=dev.circ.obj
                                    
        soqcs.dm_prnt_mtx(c_long(self.obj), format, c_double(thresh), c_long(aux))
        sys.stdout.flush()

    #---------------------------------------------------------------------------      
    # Translate the labels of a probability bins into qubit encoding. (Path encoding version).
    # Those that can not be encoded are ignored and the result is normalized.
    #---------------------------------------------------------------------------
    def translate(self, qmap, dev):
        """
        Translates the labels of a set of probability bins into a qubit encoding. (Path encoding version).
        Those that can not be encoded are ignored and the result is normalized. |br|
        **Only for ideal circuits (nm=1 and nd=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.

        :qmap(list[][]): 2xn matrix with the qubit definitions. Each column has two entries
                         specifying the channels that define the qubit.
        :dev(qodev): Device to which the set of probability bins is related.
        :return out(p_bin): A translated set of probability bins.
        
        """          
        param=to_int_ptr(qmap)
        func=soqcs.dm_translate
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[2],c_long(dev.circ.obj))  
        aux=dmatrix(1,True)
        aux.obj=obj
        return aux

    #---------------------------------------------------------------------------      
    # Translate the labels of a probability bins into qubit encoding. (Polarization encoding version).
    # Those that can not be encoded are ignored and the result is normalized.
    #---------------------------------------------------------------------------    
    def pol_translate(self, qvec, dev):
        """
        Translates the labels of a set of probability bins into a qubit encoding. (Polarization encoding version).
        Those that can not be encoded are ignored and the result is normalized. |br|
        **Only for ideal circuits (nm=1 and nd=1)**. |br|        
        **Warning!** post-selection over ignored channels may not result in a pure state. 
        Encoding is not possible under those circunstances and  a warning will be printed.

        :qvec(list[]): List with the qubit definitions. 
        :dev(qodev): Device to which the set of probability bins is related.
        :return out(p_bin): A translated set of probability bins.
        
        """          
        param=to_int_vec(qvec)
        func=soqcs.dm_pol_translate
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(dev.circ.obj))  
        aux=dmatrix(1,True)
        aux.obj=obj
        return aux
    
            

#------------------------------------------------------------------------------#                      
# Wrapper for C++ SOQCS class qodev                                            #
# Definition of a generalized device                                           #
#------------------------------------------------------------------------------#      
class auxqodev(object):
    """

    Auxiliary class implementing the functionality of a Quantum Optical Device.
    It is initialized with the same parameters than qodev.
    Intended for the internal use of the library.
    
    """
    #---------------------------------------------------------------------------      
    # Create a qodev
    #---------------------------------------------------------------------------      
    def __init__(self, nph, nch, nm=1, ns=1, np=1, dtp=-1.0, clock=0, R=0, loss=False, ckind='G', maxket=1):
        if(ckind=='G'): 
            ikind=0
        else:
            ikind=1            
        func=soqcs.dev_new_qodev
        func.restype=c_long
        self.obj = func(nph,nch,nm,ns,np,c_double(dtp),clock,R,loss,ikind,maxket)

    #---------------------------------------------------------------------------                     
    # Delete a qodev
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.dev_destroy_qodev(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Concatenates to quantum devices (with some limitations)
    #---------------------------------------------------------------------------              
    def concatenate(self,dev):
        """

        Appends the contect of a quantum device with some limitations. All input 
        photons have to be defined in the first circuit while all detectors have
        to be defined in the last. Devices with delays may also habe restrictions
        to be concatenated.
   
        :dev (qodev): Device to be concatenated.   
    
        """
        soqcs.dev_concatenate(c_long(self.obj),c_long(dev.obj))

    #---------------------------------------------------------------------------      
    # Adds a gate using other device as the gate definition
    #---------------------------------------------------------------------------    
    def add_gate(self, chlist, dev):
        """

        Adds a gate using other device as the gate definition.
   
        :chlist (list[]): List of channels to which the new gate is attached
        :dev (qodev): Device defining the gate
    
        """

        param=to_int_vec(chlist)
        soqcs.dev_add_gate(c_long(self.obj),param[0],param[1],c_long(dev.obj))
        

    #---------------------------------------------------------------------------      
    # Add photons to the device
    #---------------------------------------------------------------------------      
    def add_photons(self, N, ch, P=0, t=0.0 ,f=1.0, w=1.0):
        """
        
        Adds N photons to a device.
   
        :N (int): Number of photons to be added.
        :ch (int): Channel where the photons are added.
        :P (optional[int]): Polarization of the photons.
        :t (optional[float]): Central time of the photons wave-packet
        :f (optional[float]): Central frequency of the photons wave-packet
        :w (optional[float]): Width (Gaussians) or decay time (Exponentials) of the photon packet.
        :return  [int]: Packet number if success, a negative value if failure.
 
        """
        return soqcs.dev_add_photons(c_long(self.obj),c_int(N),c_int(ch),c_int(P),c_double(t),c_double(f),c_double(w)) 

    #---------------------------------------------------------------------------      
    # Add a QD emitter to the device
    #---------------------------------------------------------------------------      
    def add_QD(self, ch1, ch2, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0, S=0.0, k=1.0, tss=1000000.0, thv=1000000.0, cascade=0):
        """
        
        Adds a bi-excition XX and excition X photons cascade with states compatible 
        with a physical quantum dot emitter.
   
        :ch1 (int): Channel where the bi-exciton XX photon is emitted.
        :ch2 (int): Channel where the exciton X photon is emitted.
        :t1  (optional[float]): Central time of the XX photon wave-packet
        :f1  (optional[float]): Central frequency of the XX photon wave-packet
        :w1  (optional[float]): Decay time of the XX photon wave-packet
        :t2  (optional[float]): Central time of the X photon wave-packet
        :f2  (optional[float]): Central frequency of the X photon wave-packet
        :w2  (optional[float]): Decay time of the X photon wave-packet
        :S   (optional[float]): Fine Structure Splitting (FSS)
        :k   (optional[float]): Signal to noise ratio.
        :tss (optional[float]): Spin scattering time.
        :thv (optional[float]): Cross dephasing time.
        :cascade (optional[int]):  The second photon is considered to be emitted randomly after t2 to simulate a casacase 0=No/1=Yes
                                                                 
        """

        soqcs.dev_add_QD(c_long(self.obj),c_int(ch1),c_int(ch2),c_double(t1),c_double(f1),c_double(w1),c_double(t2),c_double(f2),c_double(w2),c_double(S),c_double(k),c_double(tss),c_double(thv),cascade) 

    #---------------------------------------------------------------------------      
    # Initializes the device with a path encoded Bell state        
    #---------------------------------------------------------------------------      
    def add_Bell(self, ch1, ch2, ckind, phi=0.0, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0):
        """
        
        Adds a path encoded Bell state in two channels. |br|
        
        ||Phi>=||ch1,ch2>+e^(iphi)||ch1,ch2>
   
        :ch1 (int): Channel 1 where the Bell state is defined.
        :ch2 (int): Channel 2 where the Bell state is defined.
        :ckind(optional[char]): Kind og Bell state to be created. |br|
                                '+'=||00> + ||11>  |br|
                                '-'=||00> - ||11>  |br|
                                'p'=||01> + ||10>  |br|
                                'm'=||01> - ||10>  |br|
        :phi  (optional[float]): Relative phase between the first and second ket in the definition of the Bell state.    
        :t1  (optional[float]): Central emission time of the photon in channel 1.
        :f1  (optional[float]): Central emission frequency of the photon in channel 1.
        :w1  (optional[float]): Width or decay time of the the photon in channel 1.
        :t2  (optional[float]): Central emission time of the photon in channel 2.
        :f2  (optional[float]): Central emission frequency of the photon in channel 2.
        :w2  (optional[float]): Width or decay time of the the photon in channel 2.
                                                                 
        """
        kind=0
        if(ckind=='+'): 
            kind=0
        if(ckind=='-'): 
            kind=1
        if(ckind=='p'): 
            kind=2
        if(ckind=='m'): 
            kind=3

        soqcs.dev_add_Bell(c_long(self.obj),c_int(ch1),c_int(ch2),kind,c_double(phi),c_double(t1),c_double(f1),c_double(w1),c_double(t2),c_double(f2),c_double(w2))

    #---------------------------------------------------------------------------      
    # Initialized the device with a polarization encoded Bell state
    #---------------------------------------------------------------------------      
    def add_BellP(self, ch1, ch2, ckind, phi=0.0, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0):
        """
        
        Adds a polarization encoded Bell state in two channels. |br|
        
        ||Phi>=||ch1,ch2>+e^(iphi)||ch1,ch2>
   
        :ch1 (int): Channel 1 where the Bell state is defined.
        :ch2 (int): Channel 2 where the Bell state is defined.
        :ckind(optional[char]): Kind og Bell state to be created. |br|
                                '+'=||HH> + ||VV>  |br|
                                '-'=||HH> - ||VV>  |br|
                                'p'=||HV> + ||VH>  |br|
                                'm'=||HV> - ||VH>  |br|
        :phi  (optional[float]): Relative phase between the first and second ket in the definition of the Bell state.    
        :t1  (optional[float]): Central emission time of the photon in channel 1.
        :f1  (optional[float]): Central emission frequency of the photon in channel 1.
        :w1  (optional[float]): Width or decay time of the the photon in channel 1.
        :t2  (optional[float]): Central emission time of the photon in channel 2.
        :f2  (optional[float]): Central emission frequency of the photon in channel 2.
        :w2  (optional[float]): Width or decay time of the the photon in channel 2.
                                                                 
        """
        kind=0
        if(ckind=='+'): 
            kind=0
        if(ckind=='-'): 
            kind=1
        if(ckind=='p'): 
            kind=2
        if(ckind=='m'): 
            kind=3

        soqcs.dev_add_BellP(c_long(self.obj),c_int(ch1),c_int(ch2),kind,c_double(phi),c_double(t1),c_double(f1),c_double(w1),c_double(t2),c_double(f2),c_double(w2))
    
    #---------------------------------------------------------------------------      
    # Returns the initial state
    #---------------------------------------------------------------------------      
    def input(self):  
        """
        
        Returs the photonic initial state of the quantum device.
        
        :return (state): A copy of the initial state.
                                                                 
        """
        func=soqcs.dev_input
        func.restype=c_long
        obj=func(c_long(self.obj))
        newstate=state(1,dummy=True)
        newstate.obj=obj;
        return newstate        

    #---------------------------------------------------------------------------      
    # Return the internal circuit definition
    #---------------------------------------------------------------------------      
    def circuit(self):  
        """
        
        Returs the the circuits stored in the quantum device.
        
        :return (state): A copy of the defined circuit.
                                                                 
        """
        func=soqcs.dev_circuit
        func.restype=c_long
        obj=func(c_long(self.obj))
        qoc=qocircuit(1,dummy=True)
        qoc.obj=obj
        return qoc

    #---------------------------------------------------------------------------      
    # Changes Gram-Schmidt packet order
    #---------------------------------------------------------------------------      
    def repack(self,vec):
        """
        
        Changes the Gram-Schmidt orthonormalization order of the photon packets
        
        :vec (list):  Vector with the preferred packet order.
        
        """
        param=to_int_vec(vec)
        soqcs.st_add_term(c_long(self.obj),param[0],param[1]) 

    #---------------------------------------------------------------------------      
    # Calculates the emitted packet visibility
    #---------------------------------------------------------------------------      
    def emitted_vis(self, i, j):
        """
        
        Overlapping probability of two wave packets.
        
        :i (int):  Packet i. 
        :j (int):  Packet j. 
        
        """
        func=soqcs.dev_emitted_vis
        func.restype=c_double
        return func(c_long(self.obj), i ,j)  

    #---------------------------------------------------------------------------      
    # Print packet configuration
    #---------------------------------------------------------------------------      
    def prnt_packets(self):
        """
        
        Prints the packet configuration of the quantum device.
        
        
        """
        soqcs.dev_prnt_packets(c_long(self.obj))  

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a random circuit
    #---------------------------------------------------------------------------      
    def random_circuit(self):
        """
        
        Creates a circuit defined by a random unitary matrix.
        
        
        """
        soqcs.dev_random_circuit(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a NSX subcircuit
    #---------------------------------------------------------------------------      
    def NSX(self, ch1, ch2, ch3):
        """
            
        Adds a NSX circuit element. This is the built-in version of the NSX circuit.
        Post-selection still has to be carried out to obtain the proper functionality.
        
        :ch1 (int): NSX input channel 1.
        :ch2 (int): NSX input channel 2.
        :ch3 (int): NSX input channel 3.
        
        
        """
        soqcs.dev_NSX(c_long(self.obj),c_int(ch1),c_int(ch2),c_int(ch3))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit an ideal beamsplitter
    #---------------------------------------------------------------------------      
    def beamsplitter(self, ch1, ch2, theta, phi):
        """
            
        Adds a beamsplitter to the quantum device attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        
        """
        soqcs.dev_beamsplitter(c_long(self.obj),c_int(ch1),c_int(ch2),c_double(theta),c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a dielectric beamsplitter
    #---------------------------------------------------------------------------      
    def dielectric(self, ch1, ch2, t, r):
        """
            
        Adds a dieletric film attached to channels i_ch1 and i_ch2
        It may also work as a dieletric beamsplitter.
        
        :ch1 (int): Dielectric input channel 1.
        :ch2 (int): Dielectric input channel 2.
        
        """
        soqcs.dev_dielectric(c_long(self.obj),c_int(ch1),c_int(ch2),c_double(t.real),c_double(t.imag),c_double(r.real),c_double(r.imag))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a MMI 2x2 beamsplitter
    #---------------------------------------------------------------------------      
    def MMI2(self, ch1, ch2):
        """
            
        Adds a 2x2 MMI device to the circuit attached to channels ch1 and ch2.
        
        :ch1 (int): MMI2 input channel 1.
        :ch2 (int): MMI2 input channel 2.
        
        """
        soqcs.dev_MMI2(c_long(self.obj),c_int(ch1),c_int(ch2))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a rewire gate
    #---------------------------------------------------------------------------      
    def rewire(self, ch1, ch2):
        """
            
        Adds a swaps gate between two channels
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        
        """
        soqcs.dev_rewire(c_long(self.obj),c_int(ch1),c_int(ch2))

    #---------------------------------------------------------------------------          
    # Adds to the metacircuit a phase shifter
    #---------------------------------------------------------------------------      
    def phase_shifter(self, ch, phi):
        """
            
        Adds a phase shifter to the circuit in channel ch.
        
        :ch (int): Phase shifter input channel.
        :phi (float): Angle phi in degrees.
        
        """
        soqcs.dev_phase_shifter(c_long(self.obj),c_int(ch),c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a hmogeneous lossy medium
    #---------------------------------------------------------------------------      
    def loss(self, ch, l):
        """
            
        Adds a lossy medium with loss probability l to the circuit in channel ch.
        
        :ch (int): Lossy medium input channel.
        :l  (float): Loss probability.
        
        """
        soqcs.dev_loss(c_long(self.obj),c_int(ch),c_double(l))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a delay
    #---------------------------------------------------------------------------      
    def delay(self, ch):
        """
        
        Adds a delay of one period to a channel
        
        :ch (int):  Channel where the delay is introduced.
                
        """
        soqcs.dev_delay(c_long(self.obj),c_int(ch))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a rotator
    #---------------------------------------------------------------------------      
    def rotator(self, ch, theta, phi):
        """
        
        Adds a polarization rotation device to the circuit attached to channel ch.
        
        :ch (int):  Rotator input channel.
        :theta (double):  Angle theta in degrees.
        :phi   (double):  Angle phi in degrees.
                
        """
        soqcs.dev_rotator(c_long(self.obj),c_int(ch),c_double(theta),c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a polarized beamsplitter
    #---------------------------------------------------------------------------      
    def pol_beamsplitter(self, ch1, ch2, P, theta):
        """
            
        Adds a polarized beamsplitter attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        :P   (int): Polarization to which the beamsplitter is sensitive.
        :theta (float): Effectiveness of the beamsplitter. 90º => 50/50 beamsplitter. 0º => No sensitivity.
        
        """
        soqcs.dev_pol_beamsplitter(c_long(self.obj),c_int(ch1), c_int(ch2), c_int(P),c_double(theta))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a polarized phase shifter
    #---------------------------------------------------------------------------      
    def pol_phase_shifter(self, ch, P, phi):
        """
            
        Adds a polarized phase shifter attached to channels ch.
        
        :ch  (int): Phase shifter input channel.
        :P   (int): Polarization to which the phase shifter is sensitive.
        :phi (double): Phase in degrees
        
        """
        soqcs.dev_pol_phase_shifter(c_long(self.obj),c_int(ch), c_int(P), c_double(phi))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a polarization filter
    #--------------------------------------------------------------------------- 
    def pol_filter(self, ch, P):
        """
            
        Adds a polarization filter attached to channel ch.
        
        :ch  (int): Polarization filter input channel.
        :P   (int): Polarization to be filtered.
        
        """
        soqcs.dev_pol_filter(c_long(self.obj),c_int(ch), c_int(P))
        
    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a half waveplate
    #---------------------------------------------------------------------------      
    def half(self, ch, alpha):
        """
            
        Adds a half waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double):  Rotation angle in degrees.
        
        """
        soqcs.dev_half(c_long(self.obj), c_int(ch), c_double(alpha))

    #---------------------------------------------------------------------------      
    # Adds to the metacircuit a quarter waveplate
    #---------------------------------------------------------------------------
    def quarter(self, ch, alpha):
        """
            
        Adds a half waveplate attached to channels ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double):  Rotation angle in degrees.
        
        """
        soqcs.dev_quarter(c_long(self.obj), c_int(ch), c_double(alpha))

    #---------------------------------------------------------------------------              
    # Adds to the metacircuit a detector
    #---------------------------------------------------------------------------      
    def detector(self, ch, cond=-1, pol=-1, mpi=-1, mpo=-1, eff=1.0, blnk=0.0, gamma=0.0):
        """
            
        Adds a detector to a channel of the circuit.
        :ch   (int):    Detector channel.
        :cond (int):    Detection condition.|br|
        cond>=0: The readings in the rest of the channels are accepted only in the number of photons in this channel is equal to cond. |br|
        cond=-1:  There is no condition and works as a normal detector. |br|
        cond=-2:  The channel is ignored. |br|
        :pol (optional[int]): Polarization condition. If cond>=0, pol determines the polarization of the photons to fulfill the condition. Note that if pol=-1 no assumption about the polarization of those photons is made.
        :mpi (optional[int]): Initial period of the detection window (if -1 takes the first one as default).
        :mpo (optional[int]): Final period of the detection window (if -1 takes the first one as default).
        :eff   (float):  Efficiency of the detector.
        :blnk  (float):  Ratio of time in which the detector is inactive due other detections.
        :gamma (float):  Average rate of dark counts in this channel.
        
        """
        soqcs.dev_detector(c_long(self.obj),c_int(ch),c_int(cond),c_int(pol),c_int(mpi),c_int(mpo),c_double(eff),c_double(blnk),c_double(gamma))

    #---------------------------------------------------------------------------      
    # Ignores a channel
    #---------------------------------------------------------------------------      
    def ignore(self, ch):
        """
            
        Flags the channel to be ignored in the outcome calculations. No detector is placed in this channel.
        
        :ch (int): Ignored channel.
        
        """
        soqcs.dev_detector(c_long(self.obj),c_int(ch),c_int(-2),c_int(-1),c_int(-1),c_int(-1),c_double(1.0),c_double(0.0),c_double(0.0))

    #---------------------------------------------------------------------------              
    # Adds to the output a random Gaussian whote noise of stdev^2 dispersion
    #---------------------------------------------------------------------------      
    def noise(self, stdev2):
        """
            
        Adds Gaussian white noise to the output.
        
        :stdev2 (float): Dispersion of the Gaussian noise.
        
        """
        soqcs.dev_noise(c_long(self.obj),c_double(stdev2))
    
    #---------------------------------------------------------------------------      
    #  Apply a *single ket* post-selection condition defined by the detectors 
    # ( only valid for ideal circuits).
    #---------------------------------------------------------------------------      
    def apply_condition(self, inputst, ignore):
        """
            
        Apply the post-selection condition defined by the detectors (to ideal circuits). |br|
        **Warning!**. This can be only applied to ideal circuits ns=1 where detector naturally 
        define a single ket projector. Otherwise the conditions lead to a density matrix output and
        there are other tools available in SOQCS for that purpose. Note that ignored channels are 
        not transcribed therefore post-selection may lead to collision between otherwise different kets.
    
        :return (state): Post-selected state
    
        """
        func=soqcs.dev_apply_condition
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(inputst.obj),ignore)
        newstate=state(1,dummy=True)
        newstate.obj=obj;
        return newstate


#------------------------------------------------------------------------------#
# Wrapper for C++ SOQCS class simulator                                        #
# Definition of the quantum optical circuit simulator                          #
#------------------------------------------------------------------------------#     
class simulator(object):
    """

    Simulator that can be used to calculate the outcomes of a quantum device or
    the output state of a circuit given an input state.

    :mem(optional([int]): Reserved memory for the output expressed as a maximum number of terms. (Internal memory)
    
    """
    #---------------------------------------------------------------------------          
    # Create a simulator
    #---------------------------------------------------------------------------      
    def __init__(self, mem=1000):
        func=soqcs.sim_new_simulator
        func.restype=c_long    
        self.obj = func(mem)

    #---------------------------------------------------------------------------      
    # Delete a simulator
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.sim_destroy_simulator(c_long(self.obj))

    #---------------------------------------------------------------------------               
    # Run a simulator for a metacircuit
    #---------------------------------------------------------------------------      
    def run(self, dev, method=0, nthreads=-1):
        """

        Calculates an output outcome from a device using the selected backend/core and the physical detectors definitions established in that device description.


        :dev(qodev): Input quantum device.
        :method (optional[int]): Core method selected |br|
                            0 = Direct method: The calculation is performed similarly on how it is done analytically. |br|
                            1 = Direct restricted: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            2 = Glynn method:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents. |br|
                            3 = Glynn restricted: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            4 = Ryser method:  The calculation is performed using permanents. We use the Ryser formula implemented in gray code to calculate the permanents. |br|
                            5 = Ryser restricted: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            6 = Fast Ryser method: Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output. |br|
                            7 = Fast Ryser restricted: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. |br|                            
        :nthreads (optional[int]): Number of threads to be used by Ryser methods.
        :return(p_bin): Device outcome.
    
        """
        func=soqcs.sim_run
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(dev.circ.obj),method,nthreads) 
        newoutcome=p_bin(1,True)
        newoutcome.obj=obj
        return newoutcome
   
    #---------------------------------------------------------------------------      
    # Run a simulator for an input state and a circuit
    #---------------------------------------------------------------------------      
    def run_st(self, istate,qoc,method=0,nthreads=-1,st_list=-1):
        """

        Calculates an output state as a function of an input initial state using the selected backend/core according to the rules established by a quantum circuit.

        :istate(state): Initial state.
        :qoc(qocircuit): Input quantum circuit.
        :method (optional[int]): Core method selected |br|
                            0 = Direct method: The calculation is performed similarly on how it is done analytically. |br|
                            1 = Direct restricted: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            2 = Glynn method:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents. |br|
                            3 = Glynn restricted: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            4 = Ryser method:  The calculation is performed using permanents. We use the Ryser formula implemented in gray code to calculate the permanents. |br|
                            5 = Ryser restricted: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                            6 = Fast Ryser method: Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output. |br|
                            7 = Fast Ryser restricted: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. |br|                            
        :nthreads (optional[int]): Number of threads to be used by Ryser methods.
        :st_list (optional[state]): State that contains a list of ket (of any amplitude) to be calculated by run_st. If no list is provided the full output state is obtained.
        :return(state): Output state.
    
        """

        if st_list==-1:
            func=soqcs.sim_run_state
            func.restype=c_long
            obj=func(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj),method,nthreads) 
        else:
            func=soqcs.sim_run_list
            func.restype=c_long
            obj=func(c_long(self.obj),c_long(istate.obj),c_long(st_list.obj),c_long(qoc.obj),method,nthreads) 
        newstate=state(1,True)
        newstate.obj=obj
        return newstate

    #---------------------------------------------------------------------------      
    # Run a Clifford A sampling for a device
    #---------------------------------------------------------------------------      
    def sample(self, dev, N):
        """

        Sampling of a device using Clifford A algorithm. |br|
        Proceedings of the 2018 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA). Page 146-155. SIAM Publications Library (2018). |br|
        **Warning!** Clifford A is defined to be used with a single input term. Therefor neither Bell or QD initializations are recommended.

        :dev(qodev): Input quantum device.
        :N(int): Number of samples.
        :return(p_bin): Device outcome.

        """
        func=soqcs.sim_sample
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(dev.circ.obj),N) 
        newoutcome=p_bin(1,True)
        newoutcome.obj=obj
        return newoutcome



    #---------------------------------------------------------------------------      
    # Run a Metropolis sampling for a device
    #---------------------------------------------------------------------------      
    def metropolis(self, dev, mode, N, Nburn=0, Nthin=1):
        """

        Sampling of a device using a metropolis algorithm. |br|
        Nature Physics 13, 1153-1157 (2017). |br|
        **Warning!** Metropolis is defined to be used with a single input term. Therefor neither Bell or QD initializations are recommended.

        :dev(qodev): Input quantum device.
        :mode(int): Sampling mode. |br|
                      0: f Classical. |br|
                      1: g Uniform. |br|
                      2: g Classical. |br|
                      3: f Classical ( Restricted ). |br|
                      4: g Uniform  ( Restricted ). |br|
                      5: g Classical  ( Restricted ). |br|
        :N(int): Number of samples.
        :optional(Nburn(int)): Number of initial samples to be skipped.
        :optional(Nthin(int)): Number of thinning samples.
        :return(p_bin): Device outcome.

        """
        func=soqcs.sim_metropolis
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(dev.circ.obj),mode,N,Nburn,Nthin) 
        newoutcome=p_bin(1,True)
        newoutcome.obj=obj
        return newoutcome

    #---------------------------------------------------------------------------      
    # Get a sample using one of the sampling methods available
    #---------------------------------------------------------------------------      
    def get_sample(self, dev,N,method, Nburn=0,Nthin=1,mode=0):
        """

        Get a list of samples of a device . |br|
        Nature Physics 13, 1153-1157 (2017). |br|
        **Warning!** Metropolis is defined to be used with a single input term. Therefor neither Bell or QD initializations are recommended.

        :dev(qodev): Input quantum device.
        :N(int): Number of samples.
        :method(int): Sampling method. |br|
                      0: Clifford A. |br|
                      1: Metropolis. |br|
                      5: g Classical  ( Restricted ). |br|

        :Nburn optional([int]): Number of initial samples to be skipped.
        :Nthin optional([int]): Number of thinning samples.
        :mode optional([int]): Sampling mode of the Metropolis sampler. |br|
                              0: f Classical. |br|
                              1: g Uniform. |br|
                              2: g Classical. |br|
                              3: f Classical ( Restricted ). |br|
                              4: g Uniform  ( Restricted ). |br|
                              5: g Classical  ( Restricted ). |br|        
                              
        :return(list[]): List of N samples

        """
        sample=[]                
        nlevel=dev.circuit().num_levels()
        
        if(method==0):
            func=soqcs.sim_get_sample
        if(method==1):
            func=soqcs.sim_get_metro        
        func.restype=POINTER(c_char)
        array_ptr=c_char*nlevel
        
        i=0
        istored=0
        while istored<N:
            array_ptr=func(c_long(self.obj),c_long(dev.circ.obj),mode)
            i=i+1                
            if (i>=Nburn) and (i%Nthin)==0:
                istored=istored+1
                array=[0]*nlevel
                for j in range(0,nlevel):
                    array[j]=array_ptr[j].decode('UTF-8')                
                sample.append("".join(array))

            free_ptr(array_ptr)
        
        
        
        return sample
    

#------------------------------------------------------------------------------#
# Wrapper for C++ SOQCS class mthread                                          #
# Definition of a served for paralel jobs                                      #
#------------------------------------------------------------------------------#     
class thread_server(object):
    """

    Thread server used to launch jobs in parallel.

    :mem(optional([int]): Reserved memory for the output expressed as a maximum number of terms. (Internal memory)
    
    """
    #---------------------------------------------------------------------------          
    # Create a simulator
    #---------------------------------------------------------------------------      
    def __init__(self, mem=1000):
        func=soqcs.mt_new_mthread
        func.restype=c_long    
        self.obj = func(mem)

    #---------------------------------------------------------------------------      
    # Delete a simulator
    #---------------------------------------------------------------------------      
    def __del__(self):
        soqcs.mt_destroy_mthread(c_long(self.obj))
    
    def send_st(self,istate,qoc, method=0):
        """

        Send job to the server ( Circuit version )

        :istate(int): Initial state
        :qoc(int): Quantum optical circuit to which the circuit is referred.
        :method (int): Core method selected |br|
                    0 = Direct method: The calculation is performed similarly on how it is done analytically. |br|
                    1 = Direct restricted: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    2 = Glynn method:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents. |br|
                    3 = Glynn restricted: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    4 = Ryser method:  The calculation is performed using permanents. We use the Ryser formula implemented in gray code to calculate the permanents. |br|
                    5 = Ryser restricted: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    6 = Fast Ryser method: Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output. |br|
                    7 = Fast Ryser restricted: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. |br|                            
                    
        """            
        soqcs.mt_send_work(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj),method) 


    def send(self, dev, method=0):
        """

        Send job to the server ( Device version )

        :istate(int): Initial state
        :qoc(int): Quantum optical circuit to which the circuit is referred.
        :method (int): Core method selected |br|
                    0 = Direct method: The calculation is performed similarly on how it is done analytically. |br|
                    1 = Direct restricted: Same as the direct method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    2 = Glynn method:  The calculation is performed using permanents. We use the Balasubramanian/Bax/Franklin/Glynn formula implemented in gray code to calculate the permanents. |br|
                    3 = Glynn restricted: Same as the Glynn method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    4 = Ryser method:  The calculation is performed using permanents. We use the Ryser formula implemented in gray code to calculate the permanents. |br|
                    5 = Ryser restricted: Same as the Ryser method but considering only output states of occupations by level zero or one. This restricts but speeds up the output. |br|
                    6 = Fast Ryser method: Same as the Ryser method but only those kets with non-zero contribution to post-selection are considered. This restriction speeds up the output. |br|
                    7 = Fast Ryser restricted: Same as the Fast Ryser method but considering only output states of occupations by level zero or one. This further restricts but speeds up the output. |br|                            
                    
        """            
        istate=dev.input()
        qoc=dev.circuit()
        soqcs.mt_send_work(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj),method) 
        

    def receive(self):
        """

        Receive a job from the server 

        :return(state): Output state.                
        
        """          
        func=soqcs.mt_receive_work
        func.restype=c_long
        obj=soqcs.mt_receive_work(c_long(self.obj)) 
        newstate=state(1,True)
        newstate.obj=obj
        return newstate
        


#------------------------------------------------------------------------------#      
# Class gcircuit: Graphical circuit                                            #
# Graphical representation of a quantum optical circuit                        #
#------------------------------------------------------------------------------#      
class gcircuit(object):
    """

    Auxiliary class implementing the graphical representation of a Quantum Optical Device.
    Graphical elements definitions are stored in a list that is read when the device is shown.
    Intended for internal use of the library.

    
    """
    #---------------------------------------------------------------------------      
    # Creates a gcircuit object
    #---------------------------------------------------------------------------      
    def __init__(self):              
        self.list=[]    # List
        self.nl=0       # Numbe of elements in the list
        self.newrow=0   # Number of calls to newline. Affects the plot dimensions.

    #---------------------------------------------------------------------------      
    # Concatenates two devices
    #---------------------------------------------------------------------------      
    def concatenate(self, dev):
        """

        Appends the plot of a quantum device with the current one.
   
        :dev (qodev): Plot to be concatenated.   
    
        """
        filtered=list(filter(lambda element: element[0] != 'empty', dev.list))
        # print(dev.list)
        # print(filtered)
        auxlist =  self.list + filtered
        #auxnl= self.nl+ dev.nl 
        auxnl= self.nl+ len(filtered)
        self.list=auxlist
        self.nl=auxnl
       
    #---------------------------------------------------------------------------      
    # Adds a gate to the plot specifying input channels and detectors already included 
    # in the gate definition. The gate is plotted with a custom text.
    #---------------------------------------------------------------------------            
    def add_gate(self, chlist, inlist, outlist, text):
        """

        Adds a gate to the plot pecifying input channels and detectors already included in the gate definition. 
        The gate is plotted with a custom text.
   
        :chlist (list[]): List of channels.
        :inlist (list[]): List of channels to be initialized by the gate
        :outlist (list[]): List of channels with a detector included in the gate
        :text (string): Text to be plotted inside the gate.
    
        """
        nch=len(chlist)
        nin=len(inlist)
        nout=len(outlist)

        # Drawn inputs
        for i in range(0,nin):
            self.empty(inlist[i],1)
        
        # Rewire
        ch1=chlist[0]
        for i in range(0,nch):
            ch2=chlist[i]
            if(ch2!=ch1+i):
                self.rewire(ch2,ch1+i)    
            
        # Plot gate
        self.element(ch1,nch,text, False)
        
        # Rewire
        for i in range(0,nch):
            ch2=chlist[i]
            if(ch2!=ch1+i):
                self.rewire(ch1+i,ch2)  
                
        # Plot detectors
        for i in range(0,nout):
             self.final_dec(outlist[i],-3)
        

        
    #---------------------------------------------------------------------------      
    # Defines a wire to a channel with length depth.
    #---------------------------------------------------------------------------      
    def add_wire(self,ch,depth): 
        """

        Defines an empty wire of depth cells long.
   
        :ch (int): Channel where the wire is plotted.
        :depth (int): Length in cells of the wire.
    
        """    
        self.list.append(['add_wire',ch,depth])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Initial number of photons
    #---------------------------------------------------------------------------      
    def init_ph(self,ch,n):    
        """

        Defines input number of photons indicator in a channel
   
        :ch (int): Channel where the indicator is plotted.
        :n (int):  Number of photons
    
        """          
        self.list.append(['init_ph',ch,n])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Empty channel indicator
    #---------------------------------------------------------------------------            
    def empty(self,ch,n):    
        """

        Defines a graphical representaton of an input channel indicator of an empty channel. 
        There are three possible indicators. |br|
        - Open channel: A chanel that has been left open for concatenation. |br|
        - Emtpy channel: A zero photon channel. |br|
        - Gate defined chanel: A channel that has been defined as part of the gate.|br|
        
   
        :ch (int): Channel where the indicator is plotted.
        :n (int):  -1: "Open channel" | 0="Empty channel" | >0='Gate defined channel' 
    
        """          
        self.list.append(['empty',ch,n])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Initial Bell state
    #---------------------------------------------------------------------------      
    def bell(self,ch1,ch2,kind):    
        """

        Defines a Bell state input configuration indicator to a couple of channels.
        Automatic rewire is done if the two channels are not consecutive.
   
        :ch1 (int): Bell state channel 1.
        :ch2 (int): Bell state channel 2.
        :kind (char):  Kind of Bell stte
    
        """  
        ich1=min(ch1,ch2)
        ich2=max(ch1,ch2)
        
        if(ich2!=ich1+1):
            self.rewire(ich2,ich1+1,initial=True)
        
        self.bell_element(ich1,ich2,kind)
    
        if(ich2!=ich1+1):
            self.rewire(ich1+1,ich2)

    #---------------------------------------------------------------------------                  
    # Bell state element. (Two consecutive channels, no rewires)
    #---------------------------------------------------------------------------      
    def bell_element(self,ch,auxch,kind): 
        """

        Defines a Bell state input configuration indicator to a couple of consecutive channels.
   
        :ch (int): Bell state channel 1.
        :auxch (int): Label of Bell state channel 2. (In case we habe rewire two non-consecutive ones)
        :kind (char):  Kind of Bell stte
    
        """  
        self.list.append(['bell_element',ch,auxch,kind])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------                  
    # Defines a final detector
    #---------------------------------------------------------------------------      s
    def final_dec(self,ch,post,pol=-1):
        """

        Defines a graphical representation of a detector in a channel.
   
        :ch (int): Detector channel.
        :post (int): Post-selection condition.
    
        """  
        self.list.append(['final_dec',ch,post,pol])
        self.nl=self.nl+1
        
    #---------------------------------------------------------------------------      
    # Defines a circuit element
    #---------------------------------------------------------------------------      
    def element(self,ch,nch,text, colored=False): 
        """

        Defines a graphical representation of circuit element. 
   
        :ch (int): First channel of the element.
        :nch(int): Number of channels used by the element.
        :text (string): Text describing the element.
        :colored (bool): Prints the element in a different color than the default one.
    
        """  
        self.list.append(['element',ch,nch,text, colored])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a rewire between two channels
    #---------------------------------------------------------------------------      
    def rewire(self,ch1,ch2,initial=False): 
        """

        Defines a graphical representation of a swap operation between two channels.
        There is also used in an automatic manner when optical elements referring to
        non-consecutive channels are plotted.
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        """  
        self.list.append(['rewire',ch1,ch2,initial])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a single channel circuit element
    #---------------------------------------------------------------------------      
    def one_gate(self,ch,text, colored=False): 
        """

        Defines a graphical representation of an optical element of one input channel.
        
        :ch (int): Channel.
        :text (string): Text to be plotted describing the element.
        :colored (bool): Prints the gate in a different color than the default one.
        """
        self.element(ch,1,text, colored)        

    #---------------------------------------------------------------------------              
    # Defines a two channel circuit elements
    # Rewires for non consecutive channels are calculated automatically
    #---------------------------------------------------------------------------      
    def two_gate(self,ch1,ch2,text, colored=False): 
        """

        Defines a graphical representation of an optical element of two input channels.
        If those input channels are not consecutive automatic rewires are created to
        bring the channels together grapgically and left them in their original position 
        after the graphical element is plotted.
        
        :ch1 (int): Channel1.
        :ch2 (int): Channel2.
        :text (string): Text to be plotted describing the element.
        :colored (bool): Prints the gate in a different color than the default one.
        
        """
        if(ch2!=ch1+1):
            self.rewire(ch2,ch1+1)
        
        self.element(ch1,2,text, colored)
    
        if(ch2!=ch1+1):
            self.rewire(ch1+1,ch2)

    #---------------------------------------------------------------------------      
    # Defines a three channel circuit element
    # Rewires for non consecutive channels are calculated automatically
    #---------------------------------------------------------------------------      
    def three_gate(self,ch1,ch2,ch3,text, colored=False): 
        """

        Defines  graphical representation of an optical element of three input channels.
        If those input channels are not consecutive automatic rewires are created to
        bring the channels together grapgically and left them in their original position 
        after the graphical element is plotted.
        
        :ch1 (int): Channel1.
        :ch2 (int): Channel2.
        :ch3 (int): Channel3.
        :text (string): Text to be plotted describing the element.
        :colored (bool): Prints the gate in a different color than the default one.
        
        """
        if(ch2!=ch1+1):
            self.rewire(ch2,ch1+1)

        if(ch3!=ch1+2):
            self.rewire(ch3,ch1+2)
        
        self.element(ch1,3,text, colored)
    
        if(ch2!=ch1+1):
            self.rewire(ch1+1,ch2)

        if(ch3!=ch1+2):
            self.rewire(ch1+2,ch3)

    #---------------------------------------------------------------------------      
    # Defines a delay circuit element. (It has an specific special symbol)
    #---------------------------------------------------------------------------      
    def delay(self,ch):
        """

        Defines a graphical representation of a delay. It is similat to plot a
        one gate but with a different symbol.
        
        :ch (int): Channel.
        
        """
        self.list.append(['delay',ch])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a separator
    #---------------------------------------------------------------------------      
    def separator(self): 
        """

        Defines a graphical separator between different sections of the device definition, 
        
        """
        self.list.append(['separator'])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a newline
    #---------------------------------------------------------------------------      
    def newline(self): 
        """

        Defines a change in the printing row where the elements are being drawn.
        
        """
        self.list.append(['newline'])
        self.nl=self.nl+1
        self.newrow=self.newrow+1
        
    #---------------------------------------------------------------------------              
    # Draws an empty wire in a channel for a given depth
    #---------------------------------------------------------------------------      
    def g_add_wire(self,ch,depth):
        """

        Draws a wire of depth cells long.
   
        :ch (int): Channel where the wire is plotted.
        :depth (int): Length in cells of the wire.
    
        """
        if(self.ended[ch]==0):
            init=self.display[ch];
            end=depth
        
            if(depth>init):
                height=0.5+(1+self.pad)*ch+self.line*self.width_line            
                if init>0:                               
                    plt.plot([2*init, 2*end],[height, height], color='black', linewidth=4)
                else:
                    plt.plot([0.25,2*end],[height, height], color='black', linewidth=4)
                    text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
                    plt.text(0.0,height,str(ch),**text_kwargs)
                
                self.display[ch]=depth

    #---------------------------------------------------------------------------      
    # Draws a representation of the initial photon configuration of a channel
    #---------------------------------------------------------------------------  
    def g_init_ph(self,ch,n):
        """

        Drsws an input number of photons indicator in a channel
   
        :ch (int): Channel where the wire is plotted.
        :n (int):  Number of photons
    
        """     
        self.ended[ch]=0
        self.nph[ch]=self.nph[ch]+n
        height=0.5+(1+self.pad)*ch   
        
        init=self.display[ch]
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch]
            end=init+1
        
        self.display[ch]=end 
        posx=2*init
         
                
        # Triangle
        coord = np.array([[posx+0.5,height+0.5], [posx+0.5,height-0.5], [posx+1.5 , height]])
        poly=plt.Polygon(coord,color='orchid')
        self.ax.add_patch(poly)
        plt.plot([posx+0.5,posx+0.5],[height+0.5, height-0.5], color='black', linewidth=4)
        plt.plot([posx+0.5,posx+1.5],[height-0.5, height    ], color='black', linewidth=4)
        plt.plot([posx+1.5,posx+0.5],[height    , height+0.5], color='black', linewidth=4)
        
        plt.plot([posx+0.25,posx+0.5],[height, height], color='black', linewidth=4)
        plt.plot([posx+1.5,posx+2],[height, height], color='black', linewidth=4)
        text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
        t=plt.text(posx+0.75,height,str(int(self.nph[ch])),backgroundcolor='orchid',**text_kwargs)   
        t.set_bbox(dict(alpha=0.0))
        plt.text(posx+0.0,height,str(ch),**text_kwargs)   
        
        if(self.display[ch]<=0):
            self.display[ch]=1
            
    #---------------------------------------------------------------------------      
    # Draws an empty channel indicator
    #---------------------------------------------------------------------------              
    def g_empty(self,ch,n):
        """

        Draws and empty channel indicator. There are three possible indicators. |br|
        - Open channel: A chanel that has been left open for concatenation. |br|
        - Emtpy channel: A zero photon channel. |br|
        - Gate defined chanel: A channel that has been defined as part of the gate.|br|
        
   
        :ch (int): Channel where the indicator is drawn.
        :n (int):  -1: "Open channel" | 0="Empty channel" | >0='Gate defined channel' 
    
        """    
        
        self.ended[ch]=0
        self.nph[ch]=self.nph[ch]+n
        height=0.5+(1+self.pad)*ch   
        
        init=self.display[ch]
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch]
            end=init+1
        
        self.display[ch]=end 
        posx=2*init
         
                
        # Triangle
        coord = np.array([[posx+0.5,height+0.5], [posx+0.5,height-0.5], [posx+1.5 , height]])
        if(n<0):
            poly=plt.Polygon(coord,color='darkorchid')
        if(n==0):
            poly=plt.Polygon(coord,color='lightblue')
        if(n>0):
            poly=plt.Polygon(coord,color='darkcyan')

        self.ax.add_patch(poly)
        plt.plot([posx+0.5,posx+0.5],[height+0.5, height-0.5], color='black', linewidth=4)
        plt.plot([posx+0.5,posx+1.5],[height-0.5, height    ], color='black', linewidth=4)
        plt.plot([posx+1.5,posx+0.5],[height    , height+0.5], color='black', linewidth=4)
        
        plt.plot([posx+0.25,posx+0.5],[height, height], color='black', linewidth=4)
        plt.plot([posx+1.5,posx+2],[height, height], color='black', linewidth=4)
        text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
        if(n<0):
            t=plt.text(posx+0.75,height,'x',backgroundcolor='orchid',**text_kwargs)   
        if(n==0):            
            t=plt.text(posx+0.75,height,'0',backgroundcolor='white',**text_kwargs)               
        if(n>0):            
            t=plt.text(posx+0.75,height,'G',backgroundcolor='white',**text_kwargs)               

        t.set_bbox(dict(alpha=0.0))
        plt.text(posx+0.0,height,str(ch),**text_kwargs)   
        
        if(self.display[ch]<=0):
            self.display[ch]=1
            
    #---------------------------------------------------------------------------      
    # Draws a Bell state element. 
    #---------------------------------------------------------------------------      
    def g_bell_element(self,ch,auxch,kind): 
        """

        Draws a Bell state input configuration indicator to a couple of consecutive channels.
   
        :ch (int): Bell state channel 1.
        :auxch (int): Label of Bell state channel 2. (In case we have rewired two non-consecutive ones)
        :kind (char):  Kind of Bell stte
    
        """  
        self.ended[ch]=0
        self.ended[ch+1]=0
        nch=2
        maxdepth=0
        for i in range(ch,ch+nch):
            maxdepth=max(maxdepth,self.display[i])
          
        init=maxdepth
        end=init+1;       

        for i in range(ch,ch+nch):
            self.g_add_wire(i,init)
            self.display[i]=end                                
      
        posx=2*init
        height=0.5+(1+self.pad)*ch            
        height2=0.5+(1+self.pad)*(ch+1) 
        add=(nch-1)      
              
        rect = patches.Rectangle(((2*init+0.2), (height-0.5+0.0)),
                                  (2*end-0.2-(2*init+0.2)), 
                                  ((height+add+0.6-0.0)-(height-0.5+0.0)),
                                  linewidth=3, edgecolor='black', linestyle='--',facecolor='thistle')
        self.ax.add_patch(rect)
      
        # Triangle #1
        coord = np.array([[2*init+0.5,height+0.5], [2*init+0.5,height-0.5], [2*init+1.5 , height]])
        poly=plt.Polygon(coord,color='orchid')
        self.ax.add_patch(poly)
        plt.plot([2*init+0.5,2*init+0.5],[height+0.5, height-0.5], color='black', linewidth=4)
        plt.plot([2*init+0.5,2*init+1.5],[height-0.5, height    ], color='black', linewidth=4)
        plt.plot([2*init+1.5,2*init+0.5],[height    , height+0.5], color='black', linewidth=4)
        
      
        # Triangle #2
        coord = np.array([[2*init+0.5,height2+0.5], [2*init+0.5,height2-0.5], [2*init+1.5 , height2]])
        poly=plt.Polygon(coord,color='orchid')
        self.ax.add_patch(poly)
        plt.plot([2*init+0.5,2*init+0.5],[height2+0.5, height2-0.5], color='black', linewidth=4)
        plt.plot([2*init+0.5,2*init+1.5],[height2-0.5, height2    ], color='black', linewidth=4)
        plt.plot([2*init+1.5,2*init+0.5],[height2    , height2+0.5], color='black', linewidth=4)
        
        #Lines 1
        plt.plot([2*init+0.25,2*init+0.5],[height, height], color='black', linewidth=4)
        plt.plot([2*init+1.5,2*init+2],[height, height], color='black', linewidth=4)
      
        #Lines 2
        plt.plot([2*init+0.25,2*init+0.5],[height2, height2], color='black', linewidth=4)
        plt.plot([2*init+1.5,2*init+2],[height2, height2], color='black', linewidth=4)
      
        if(kind=='O'):
            text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
            t=plt.text(2*init+0.9,height,"0/1",backgroundcolor='orchid',**text_kwargs)   
            t.set_bbox(dict(alpha=0.0))
            plt.text(0.0,height,str(ch),**text_kwargs)   
 
            text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
            t=plt.text(2*init+0.9,height2,"0/1",backgroundcolor='orchid',**text_kwargs)   
            t.set_bbox(dict(alpha=0.0))
            plt.text(posx+0.0,height2,str(auxch),**text_kwargs)   

        if(kind=='P'):
            text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
            t=plt.text(2*init+0.9,height,"H/V",backgroundcolor='orchid',**text_kwargs)   
            t.set_bbox(dict(alpha=0.0))
            plt.text(0.0,height,str(ch),**text_kwargs)   
 
            text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
            t=plt.text(2*init+0.9,height2,"H/V",backgroundcolor='orchid',**text_kwargs)   
            t.set_bbox(dict(alpha=0.0))
            plt.text(posx+0.0,height2,str(auxch),**text_kwargs) 
                                 
    #---------------------------------------------------------------------------              
    # Draws a representation of the detector configuration in a channel
    #---------------------------------------------------------------------------      
    def g_final_dec(self,ch,post,pol):
        """

        Draws a graphical representation of a detector in a channel.
   
        :ch (int): Detector channel.
        :post (int): Post-selection condition.
    
        """ 
        self.ended[ch]=1
        init=self.display[ch]
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch]
            end=init+1

        self.display[ch]=end
        height=0.5+(1+self.pad)*ch+self.line*self.width_line                        
        
        if(post>=-1): 
            color="steelblue"
        
        if(post==-2): 
            color="red"     
        
        if(post==-3): 
            color="darkcyan"
            
            
        circle = plt.Circle((2*init+0.5,height),0.5, fc=color,ec="black",linewidth=4)
        self.ax.add_patch(circle)


        rect = patches.Rectangle((2*init, height-0.5),
                                  ((2*init+0.5)-(2*init)), 
                                  ((height+0.5)-(height-0.5)),
                                  linewidth=4, edgecolor='black', facecolor=color)
        self.ax.add_patch(rect)

        rect = patches.Rectangle((2*init+0.45, height-0.47),
                                  ((2*init+0.55)-(2*init+0.45)), 
                                  ((height+0.47)-(height-0.47)),
                                  linewidth=0, edgecolor='black', facecolor=color)
        self.ax.add_patch(rect)
        text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
        plt.text(2*init+1.5,height,str(ch),**text_kwargs)   
           
        if(pol==0): 
            strpol='H'
        if(pol==1): 
            strpol='V'
        if(pol>1): 
            strpol='P'+str(pol)
            
        if(post>=0) and (pol<0):
            plt.text(2*init+0.4,height,str(post),backgroundcolor=color,**text_kwargs)   
            
        if(post>=0) and (pol>=0):
            plt.text(2*init+0.4,height,str(post)+strpol,backgroundcolor=color,**text_kwargs)   
            
        if(post==-3):
            plt.text(2*init+0.4,height,'G',backgroundcolor=color,**text_kwargs)   

    #---------------------------------------------------------------------------      
    # Draws a representation of a general circuit element
    #---------------------------------------------------------------------------      
    def g_element(self,ch,nch,text, colored=False): 
        """

        Draws a graphical representation of circuit element. 
        Essentially a box with some text describing its functionality. 
        All the channels have to appear consecutively.
   
        :ch (int): First channel of the element.
        :nch(int): Number of channels used by the element.
        :text (string): Text describing the element.
        :colored (bool): Prints the graphical element in a different color than the default one.
    
        """  
        maxdepth=0
        for i in range(ch,ch+nch):
            maxdepth=max(maxdepth,self.display[i])
            self.ended[i]=0
            
        init=maxdepth
        end=init+1;    
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch];
            end=init+1;
        
        for i in range(ch,ch+nch):
            self.g_add_wire(i,init)
            self.display[i]=end                                
            
        height=0.5+(1+self.pad)*ch+self.line*self.width_line                        
        add=(nch-1)      

        color='mediumpurple'
        if colored==True:
            color='mediumorchid'            
        rect = patches.Rectangle(((2*init+0.2), (height-0.5+0.2)),
                                  (2*end-0.2-(2*init+0.2)), 
                                  ((height+add+0.5-0.2)-(height-0.5+0.2)),
                                 linewidth=4, edgecolor='black', facecolor=color)
        self.ax.add_patch(rect)

        
        for i in range(0,nch):
            heightch=0.5+(1+self.pad)*(ch+i)+self.line*self.width_line            
            
            plt.plot([2*init,(2*init+0.2)],[heightch, heightch], color='black', linewidth=4)
            plt.plot([(2*end-0.2),2*end],[heightch, heightch], color='black', linewidth=4)
                        
        
        text_kwargs = dict(ha='center', va='center',multialignment='left', fontsize=self.size_font-2, color='black')
        plt.text((2*init+1),(height+(nch-1)/2),str(text),**text_kwargs)

    #---------------------------------------------------------------------------      
    # Draws a rewire between two channels
    #---------------------------------------------------------------------------      
    def g_rewire(self,ch1,ch2,initial): 
        """

        Draws a graphical representation of a swap operation between two channels.
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        """  
        maxdepth=0
        chi=min(ch1,ch2)
        cho=max(ch1,ch2)
        for chaux in range(chi,cho+1):
            maxdepth=max(maxdepth,self.display[chaux])
            
        if maxdepth>0:
            init=maxdepth
            end=init+1;
            if end>(self.depth-2):
                self.g_newline()
                init=self.display[ch1];
                end=init+1;

            self.g_add_wire(chi,init)   
            self.display[chi]=end        
            self.g_add_wire(cho,init)   
            self.display[cho]=end        
            for chaux in range(chi+1,cho):
                self.g_add_wire(chaux,end)   
                self.display[chaux]=end
                                
            height1=0.5+(1+self.pad)*chi+self.line*self.width_line                       
            height2=0.5+(1+self.pad)*cho+self.line*self.width_line                       
        
            if self.ended[chi]==0:
                plt.plot([2*init,2*end],[height1, height2], color='black', linewidth=4)
            if self.ended[cho]==0:
                plt.plot([2*init,2*end],[height2, height1], color='black', linewidth=4)   
        
        aux=self.ended[chi]
        self.ended[chi]=self.ended[cho]
        self.ended[cho]=aux
        
    #---------------------------------------------------------------------------                  
    # Draws a delay circuit element. (It has an specific special symbol)
    #---------------------------------------------------------------------------      
    def g_delay(self,ch):
        """

        Draws a graphical representation of a delay. It is similat to plot a
        one gate but with a different symbol.
        
        :ch (int): Channel.
        
        """
        init=self.display[ch];
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch];
            end=init+1;
            
        self.display[ch]=end;
        height=0.5+(1+self.pad)*ch+self.line*self.width_line                            
        
        plt.plot([2*init    , 2*init+0.5],[height+0.3, height+0.3], color='black', linewidth=4)
        plt.plot([2*init+0.5, 2*init+1  ],[height    , height    ], color='black', linewidth=4)
        plt.plot([2*init+1  , 2*init+1.5],[height+0.3, height+0.3], color='black', linewidth=4)
        plt.plot([2*init+1.5, 2*init+2  ],[height    , height    ], color='black', linewidth=4)
        
        plt.plot([2*init    , 2*init    ],[height+0.3, height    ], color='black', linewidth=4)
        plt.plot([2*init+0.5, 2*init+0.5],[height+0.3, height    ], color='black', linewidth=4)
        plt.plot([2*init+1  , 2*init+1.0],[height+0.3, height    ], color='black', linewidth=4)
        plt.plot([2*init+1.5, 2*init+1.5],[height+0.3, height    ], color='black', linewidth=4)
    
    
    #---------------------------------------------------------------------------              
    # Draws a separator to clarify the plot
    #---------------------------------------------------------------------------      
    def g_separator(self): 
        """

        Draws a graphical separator between different sections of the device definition, 
        
        """
        maxdepth=0
        for i in range(0,self.nch):
            maxdepth=max(maxdepth,self.display[i])

        init=maxdepth
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[0];
            end=init+1;
               
        for i in range(0,self.nch):
            self.g_add_wire(i,end)
            self.display[i]=end
            
        plt.plot([2*init+1, 2*init+1],
                 [0+self.line*self.width_line , self.nch+self.line*self.width_line], 
                 color='gray',linestyle='dashed', linewidth=6)

    #---------------------------------------------------------------------------              
    # Draws a "hard" separator to identify changes of line
    #---------------------------------------------------------------------------      
    def g_hard_separator(self, width): 
        """

        Draws a different graphical separator to indicate a change of printing row.
        
        """
        maxdepth=0
        for i in range(0,self.nch):
            maxdepth=max(maxdepth,self.display[i])

        init=maxdepth
        end=init+1;
               
        for i in range(0,self.nch):
            self.g_add_wire(i,end)
            self.display[i]=end
            
        plt.plot([2*init+1, 2*init+1],
                 [0+self.line*self.width_line , self.nch+self.line*self.width_line], 
                 color='black',linestyle='dashdot', linewidth=width)

    #---------------------------------------------------------------------------              
    # Advances the printing row.
    #---------------------------------------------------------------------------           
    def g_newline(self): 
        """

        Changes the row where the elements are being drawn to a new line and plots
        two "hard" separators to indicate this change.
        
        """
        self.g_hard_separator(3)
        self.display=np.zeros(self.nch)
        self.line=self.line+1
        self.g_hard_separator(3)
        for i in range(0,self.nch):
            self.g_add_wire(i,1)
                    
    #---------------------------------------------------------------------------          
    # Plots and launches the widget with the canvas of the whole circuit
    #---------------------------------------------------------------------------      
    def show(self, nch, depth, sizexy, font, slin=0, rows=1):
        """
        Shows the quantum device circuit. |br|
        It reads the list of definition and draws the graphical elements. Each channels 
        works like a stack where graphical elements can be added the top. Each graphical
        elements occupies a cell. The depth of the plot is determined by the maximum number 
        of cells occupied for all channels.
        
        :nch (int): Number of channels.
        :depth (int): Depth (in cells) of the circuit to be drawn.
        :sizexy (int): Number of pixels by cell.
        :font(int): Font size of the labels.
        :slin(optional(int)): Starting lines 0='No'/1='Yes'. If starting lines is 0 no lines are shown until the channel is initialized.
        :rows(optional(int)): Numbers of printing lines to be used by the plotter. |br| 
                   If rows=1, the plot is printed in one line unless hard separators are used to manually 
                   point out the change of printing line. |br| 
                   If rows>1, the plot is drawn in a fixed quantity of lines. In this case the change of line is done automatically. |br|
                   
        """
  
        dpi=sizexy;
        
        self.nch=nch
        self.size_font=int(font)        
        self.display=np.zeros(nch)
        self.nph=np.zeros(nch)
        self.pad=0.1
        self.line=0
        self.width_line=((sizexy+sizexy/10)*(nch-1)+sizexy)/dpi+0.5*sizexy/dpi
        self.depth=depth
        if rows>1:
            self.newrow=rows-1
        self.ended=[1]*nch
        if slin==0 :
            self.ended=[1]*nch
        else:
            self.ended=[0]*nch
        
        self.fig, self.ax= plt.subplots()        
        self.fig.set_dpi(dpi)
        self.fig.set_size_inches( ((depth-1)*2*sizexy+sizexy)/dpi, (self.newrow+1)*((sizexy+sizexy/10)*nch+sizexy)/dpi, forward=True)
        self.ax.set_facecolor("white")
        plt.xlim([0, ((depth-1)*2*sizexy-0.8*sizexy)/dpi])
        plt.ylim([-0.5*sizexy/dpi, (self.newrow+1)*self.width_line])
        plt.gca().invert_yaxis()
        plt.axis('off')
        
        for i in range(0,self.nl):
            gate=self.list[i]
            
            if(gate[0]=='add_wire'):             
                self.g_add_wire(gate[1],gate[2])

            if(gate[0]=='init_ph'): 
                self.g_init_ph(gate[1],gate[2])
                
            if(gate[0]=='empty'): 
                self.g_empty(gate[1],gate[2])
                
            if(gate[0]=='bell_element'): 
                self.g_bell_element(gate[1],gate[2],gate[3])

            if(gate[0]=='final_dec'): 
                self.g_final_dec(gate[1],gate[2],gate[3])

            if(gate[0]=='element'): 
                self.g_element(gate[1],gate[2],gate[3], gate[4])

            if(gate[0]=='rewire'): 
                self.g_rewire(gate[1],gate[2],gate[3])
          
            if(gate[0]=='delay'): 
                self.g_delay(gate[1])

            if(gate[0]=='separator'): 
                self.g_separator()
                
            if(gate[0]=='newline'): 
                self.g_newline()

        plt.draw()
        plt.show()

        
#------------------------------------------------------------------------------#                 
# Class qodev: Full metacircuit                                                #
# It merges the classes auxqodev and qcircuit to create a full device          #
# representation# in python able to create both the logical and the graphical  #
# representations of the device                                                #
#------------------------------------------------------------------------------#      
class qodev(object):
    """

    A quantum optical device consists in a quantum optical circuit and a photon abstraction of its initial state.

    :nph (int): Maximum number of photons
    :nch (int): Number of channels
    :nm (optional[int]): Number of modes
    :ns (optional[int]): Number of packets
    :np (optional[int]): Number of periods
    :dtp (optional[double]): Length of the periods
    :clock (optional[int]):  Detector behavior configuration: |br|
                            0: Counter. The detectors behave as counters. |br|
                            1: Time. The detectors are able distinguish arrival times but they are blind to frequency. (This mode may change packet numeration). |br|
                            2: Full.  The detectors are able to distinguish arrival times and frequency of the outgoing photons. |br|
                            3: Manual mode.  Like 1. But the user is fully responsible of the packet definitions in order to obtain the correct calculations. |br|
                            4: Period.  Detectors can only distinguish the period in which photons arrive. |br|
    :R (optional[int]): Number of iterations in the calculation of detector dead time and dark counts effects.
    :loss (optional[bool]): Explicit computation of losses.
    :ckind (optional[char]): Packet shape: |br|
                             'G': Gaussian |br|
                             'E': Exponential |br|
    :mem(optional(int)): Number of reserved entries to define the initial state. (Internal memory)
    
    """
    #---------------------------------------------------------------------------          
    # Create a device
    #---------------------------------------------------------------------------      
    def __init__(self, nph, nch, nm=1, ns=1, np=1, dtp=-1.0, clock=0, R=0,   loss=False, ckind='G', mem=4):
        self.circ=auxqodev(nph,nch,nm,ns,np,dtp,clock,R,loss,ckind,mem)
        self.disp=gcircuit()
        self.first=1
        self.nph=nph
        self.nch=nch
        self.nm=nm     # Needed for the compiler.
        self.inlist=[]
        self.outlist=[]
        

    #---------------------------------------------------------------------------      
    # Concatenates two devices ( with some limitations)
    #---------------------------------------------------------------------------      
    def concatenate(self,dev):
        """

        Appends the content of two quantum devices with some limitations. 
        All the input photons have to be defined in the first circuit and 
        all detectors in the last. Circuits with delays cannot be concatenated either.   

        :dev (qodev): Device to be concatenated.   
    
        """
        self.disp.concatenate(dev.disp)
        self.circ.concatenate(dev.circ)
    #---------------------------------------------------------------------------      
    # Alias to concatenate
    #---------------------------------------------------------------------------       
    def __lshift__(self, dev):
        self.concatenate(dev)
    
    #---------------------------------------------------------------------------      
    # Adds a gate using other device as the gate definition
    #---------------------------------------------------------------------------   
    def add_gate(self, chlist, dev,text=''):
        """

        Adds a gate using other device as the gate definition.
   
        :chlist (list[]): List of channels to which the new gate is attached
        :dev (qodev): Device defining the gate
        :text (string): Text with the name given to the gate ( to be plotted )
    
        """
        
        inlist=[]
        for i in range(0,len(dev.inlist)):
            self.inlist.append(chlist[dev.inlist[i]])
            inlist.append(chlist[dev.inlist[i]])

        outlist=[]
        for i in range(0,len(dev.inlist)):
            self.outlist.append(chlist[dev.outlist[i]])    
            outlist.append(chlist[dev.outlist[i]])


        self.disp.add_gate(chlist,inlist,outlist,text)
        self.circ.add_gate(chlist,dev.circ)
            
    #---------------------------------------------------------------------------              
    # Adds photons to the device
    #---------------------------------------------------------------------------      
    def add_photons(self,  N, ch, P=0, t=0.0 ,f=1.0, w=1.0):
        """
        
        Adds N photons to the quantum optical device.
   
        :N (int): Number of photons to be added.
        :ch (int): Channel where the photons are added.
        :P (optional[int]): Polarization of the photons.
        :t (optional[float]): Central time of the photons wave-packet
        :f (optional[float]): Central frequency of the photons wave-packet
        :w (optional[float]): Width (Gaussians) or decay time (Exponentials) of the photon packet.
        :return  [int]: Packet number if success, a negative value if failure.
 
        """
        aux=self.circ.add_photons(N, ch, P ,t, f, w)
        self.disp.init_ph(ch,N)
        self.inlist.append(ch)
        return aux
    
    #---------------------------------------------------------------------------              
    # Adds an open channel to the device
    #---------------------------------------------------------------------------    
    def open_channel(self,  ch):
        """
        
        Adds 0 photons to the device, no packet is defined and the indicator of channel
        to be concatenated is drawn if the circuit is plotted.
   
        :ch (int): Channel that is left open.
 
        """
        self.disp.empty(ch,-1)

    #---------------------------------------------------------------------------              
    # Adds an empty channel to the device
    #---------------------------------------------------------------------------            
    def empty_channel(self,  ch):
        """
        
        Adds 0 photons to the device, no packet is defined and the indicator 
        of empty channel is drawn if the circuit is plotted. Note that this is different
        than adding 0 photons with add_photons because in this case no packet is defined.
   
        :ch (int): Channel that is left empty.
 
        """
        self.disp.empty(ch,0)
 
    #---------------------------------------------------------------------------              
    # Initialize qubits (Path encoding)
    #---------------------------------------------------------------------------            
    def qubits(self, qlist, ancilla, qmap):
        """
        
        Initialize the device with photons that represent a list of path encoded qubits.
   
        :qlist (list[int]): Qubit values
        :ancilla (list[int]): List with the values to initalize the ancillas. Values are assigned in order from small to large channel number.
        :qmap (list[][]): 2xn matrix with the qubit definitions. Each column has two entries
                          specifying the channels that define the qubit.
                          
        """
        
        channels=[0]*self.nch
        for i in range(0,len(qlist)):
            channels[qmap[0][i]]=1
            channels[qmap[1][i]]=1
            if qlist[i]==0:
                self.circ.add_photons(0, qmap[0][i])
                self.circ.add_photons(1, qmap[1][i])                
                self.disp.init_ph(qmap[0][i],0)
                self.disp.init_ph(qmap[1][i],1)
            else:
                self.circ.add_photons(1, qmap[0][i])
                self.circ.add_photons(0, qmap[1][i])                
                self.disp.init_ph(qmap[0][i],1)
                self.disp.init_ph(qmap[1][i],0)
                           
        j=0
        for i in range(0,self.nch):
            if channels[i]==0:
                self.circ.add_photons(ancilla[j], i)
                self.disp.init_ph(i,ancilla[j])
                j=j+1
        
    #---------------------------------------------------------------------------              
    # Initialize qubits (Polarization encoding)
    #---------------------------------------------------------------------------            
    def pol_qubits(self, qlist, ancilla, qmap): ## Not tested
        """
        
        Initialize the device with photons that represent a list of polarization encoded qubits.
   
        :qlist (list[int]): Qubit values
        :ancilla (list[int]): List with the values to initalize the ancillas. Values are assigned in order from small to large channel number.
        :qmap(list[]): List with the qubit definitions. 
                          
        """
        
        channels=[0]*self.nch
        for i in range(0,len(qlist)):
            channels[qmap[i]]=1
            self.circ.add_photons(1, qmap[i],qlist[i])                
            self.disp.init_ph(qmap[i],1)
   
                           
        j=0
        for i in range(0,self.nch):
            if channels[i]==0:
                self.circ.add_photons(1,i,ancilla[j])
                self.disp.init_ph(i,1)
                j=j+1
                
    #---------------------------------------------------------------------------               
    # Adds Quantum Dot emitter to the device
    #---------------------------------------------------------------------------      
    def add_QD(self, ch1, ch2, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0, S=0.0, k=1.0, tss=1000000.0, thv=1000000.0, cascade=0):
        """
        
        Adds a pair of photons to the input of a circuit as if they where emitted by a quantum dot in a XX-X cascade.
   
        :ch1 (int): Channel where the bi-exciton XX photon is emitted.
        :ch2 (int): Channel where the exciton X photon is emitted.
        :t1  (optional[float]): Central time of the XX photon wave-packet
        :f1  (optional[float]): Central frequency of the XX photon wave-packet
        :w1  (optional[float]): Decay time of the XX photon wave-packet
        :t2  (optional[float]): Central time of the X photon wave-packet
        :f2  (optional[float]): Central frequency of the X photon wave-packet
        :w2  (optional[float]): Decay time of the X photon wave-packet
        :S   (optional[float]): Fine Structure Splitting (FSS)
        :k   (optional[float]): Signal to noise ratio.
        :tss (optional[float]): Spin scattering time.
        :thv (optional[float]): Cross dephasing time.
        :cascade (optional[int]):  The second photon is considered to be emitted randomly after t2 to simulate a casacase 0=No/1=Yes                                                                 
        
        """
        self.circ.add_QD(ch1, ch2 ,t1, f1, w1, t2, f2, w2, S, k, tss, thv, cascade)
        self.disp.bell(ch1,ch2,'P')
        self.inlist.append(ch1)
        self.inlist.append(ch2)

    #---------------------------------------------------------------------------      
    # Initializes the device with a path encoded Bell state
    #---------------------------------------------------------------------------              
    def add_Bell(self, ch1, ch2, ckind, phi=0.0, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0):
        """
        
        Adds a path encoded Bell state as a device input. |br|
        
        ||Phi>=||ch1,ch2>+e^(iphi)||ch1,ch2>
   
        :ch1 (int): Channel 1 where the Bell state is defined.
        :ch2 (int): Channel 2 where the Bell state is defined.
        :ckind(optional[char]): Kind og Bell state to be created. |br|
                                '+'=||00> + ||11>  |br|
                                '-'=||00> - ||11>  |br|
                                'p'=||01> + ||10>  |br|
                                'm'=||01> - ||10>  |br|
        :phi  (optional[float]): Relative phase between the first and second ket in the definition of the Bell state.    
        :t1  (optional[float]):  Central emission time of the photon in channel 1.
        :f1  (optional[float]):  Central emission frequency of the photon in channel 1.
        :w1  (optional[float]):  Width or decay time of the the photon in channel 1.
        :t2  (optional[float]):  Central emission time of the photon in channel 2.
        :f2  (optional[float]):  Central emission frequency of the photon in channel 2.
        :w2  (optional[float]):  Width or decay time of the the photon in channel 2.
                                                                 
        """
        self.circ.add_Bell(ch1, ch2, ckind, phi ,t1,f1, w1, t2, f2, w2)
        self.disp.bell(ch1,ch2,'O')
        self.inlist.append(ch1)
        self.inlist.append(ch2)

    #---------------------------------------------------------------------------      
    # Initializes the device with a polarization encoded Bell state
    #---------------------------------------------------------------------------              
    def add_BellP(self, ch1, ch2, ckind, phi=0.0, t1=0.0, f1=1.0, w1=1.0, t2=0.0, f2=1.0, w2=1.0):
        """
        
        Adds a polarization encoded Bell state as a device input. |br|
        
        ||Phi>=||ch1,ch2>+e^(iphi)||ch1,ch2>
   
        :ch1 (int): Channel 1 where the Bell state is defined.
        :ch2 (int): Channel 2 where the Bell state is defined.
        :ckind(optional[char]): Kind of Bell state to be created. |br|
                                '+'=||HH> + ||VV>  |br|
                                '-'=||HH> - ||VV>  |br|
                                'p'=||HV> + ||VH>  |br|
                                'm'=||HV> - ||VH>  |br|
        :phi  (optional[float]): Relative phase between the first and second ket in the definition of the Bell state.    
        :t1  (optional[float]):  Central emission time of the photon in channel 1.
        :f1  (optional[float]):  Central emission frequency of the photon in channel 1.
        :w1  (optional[float]):  Width or decay time of the the photon in channel 1.
        :t2  (optional[float]):  Central emission time of the photon in channel 2.
        :f2  (optional[float]):  Central emission frequency of the photon in channel 2.
        :w2  (optional[float]):  Width or decay time of the the photon in channel 2.
                                                                 
        """
        self.circ.add_BellP(ch1, ch2, ckind, phi ,t1,f1, w1, t2, f2, w2)
        self.disp.bell(ch1,ch2,'P')
        self.inlist.append(ch1)
        self.inlist.append(ch2)

    #---------------------------------------------------------------------------              
    # Returns the input state of the device
    #---------------------------------------------------------------------------      
    def input(self):  
        """
        
        Returns a copy of the photon initial state.
        
        :return (state): A copy of the initial state.
                                                                 
        """
        return self.circ.input()

    #---------------------------------------------------------------------------      
    # Returns the internal logical circuit representation
    #---------------------------------------------------------------------------      
    def circuit(self):  
        """
        
        Returns a quantum optical circuit equivalent to the quantum device except for the initial conditions that are not defined in a circuit.
        
        :return (qocircuit): A copy of the defined circuit.
                                                                 
        """
        return self.circ.circuit()

    #---------------------------------------------------------------------------      
    # Adds a separator to the circuit to see a more clear plot
    #---------------------------------------------------------------------------      
    def separator(self): 
        """
        
        Adds a separator to the circuit to see a more clear plot.
                                                                 
        """
        self.disp.separator()    
        
    #---------------------------------------------------------------------------      
    # Adds a change of line to see a more clear plot
    #---------------------------------------------------------------------------      
    def newline(self): 
        """
        
        Adds a change in the printing row where the elements are being drawn to
        see a more clear plot.
                                                                 
        """
        self.disp.newline()    

    #---------------------------------------------------------------------------      
    # Changes Gram-Schmidt packet order
    #---------------------------------------------------------------------------      
    def repack(self,vec):
        """
        
        Changes the Gram-Schmidt orthonormalization order of the photon packets
        
        :vec (list):  List with the preferred packet order.
        
        """
        return self.circ.repack(vec)

    #---------------------------------------------------------------------------      
    # Calculates the emitted packet visibility
    #---------------------------------------------------------------------------      
    def emitted_vis(self, i, j):
        """
        
        Overlap probability of two wavepackets.
        
        :i (int):  Packet i. 
        :j (int):  Packet j. 
        :return (double): Overlap probability
        
        """
        return self.circ.emitted_vis(i, j)

    #---------------------------------------------------------------------------      
    # Print packet configuration
    #---------------------------------------------------------------------------      
    def prnt_packets(self):
        """
        
        Prints the packet configuration of the quantum device.
        
        
        """
        self.circ.prnt_packets()

    #---------------------------------------------------------------------------      
    # Adds to the device a random circuit
    #---------------------------------------------------------------------------      
    def random_circuit(self, ):
        """
        
        Creates a circuit defined by a random unitary matrix.
        
        
        """
        strf="     RND"
        self.circ.random_circuit()
        self.disp.element(0, self.nch, strf)

    #---------------------------------------------------------------------------                      
    # Adds to the device a NSX subcircuit    
    #---------------------------------------------------------------------------      
    def NSX(self, ch1, ch2, ch3):
        """
            
        Adds a NSX circuit element. Post-selection still has to be carried out to obtain the proper functionality.
        
        :ch1 (int): NSX input channel 1.
        :ch2 (int): NSX input channel 2.
        :ch3 (int): NSX input channel 3.
        
        
        """
        strf="     NSX"
        self.circ.NSX(ch1,ch2,ch3)
        self.disp.three_gate(ch1, ch2, ch3,strf)

    #---------------------------------------------------------------------------      
    # Adds to the device an ideal beamsplitter
    #---------------------------------------------------------------------------      
    def beamsplitter(self, ch1, ch2, theta, phi, colored= False ):
        """
            
        Adds a beamsplitter to the quantum device attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        :theta (double):  Angle theta in degrees.
        :phi   (double):  Angle phi in degrees.
        :colored (bool): Prints the beamsplitter in a different color than the default one.
        
        """
        str1="     BS\n \u03B8="
        str2=str(round(theta,2))
        str3="º\n \u03C6="
        str4=str(round(phi,2))
        str5="º"
        strf=str1+str2+str3+str4+str5
        self.circ.beamsplitter(ch1,ch2,theta,phi)
        self.disp.two_gate(ch1, ch2,strf, colored)

    #---------------------------------------------------------------------------      
    # Adds to the device a dielectric beamsplitter
    #---------------------------------------------------------------------------      
    def dielectric(self, ch1, ch2, t, r):
        """
            
        Adds a dieletric film attached to channels ch1 and ch2.
        It may also work as a dieletric beamsplitter.
        
        :ch1 (int): Dielectric input channel 1.
        :ch2 (int): Dielectric input channel 2.
        
        """
        str1="     Di\n y="
        str2=str(t)
        str3="\n r="
        str4=str(r)
        strf=str1+str2+str3+str4
        self.circ.dielectric(ch1,ch2,t,r)
        self.disp.two_gate(ch1, ch2,strf)

    #---------------------------------------------------------------------------      
    # Adds to the device a MMI 2x2 beamsplitter
    #---------------------------------------------------------------------------      
    def MMI2(self, ch1, ch2):
        """
            
        Adds a 2x2 MMI to the device attached to channels ch1 and ch2.
        
        :ch1 (int): MMI2 input channel 1.
        :ch2 (int): MMI2 input channel 2.
        
        """
        strf="     MMI"
        self.circ.MMI2(ch1,ch2)
        self.disp.two_gate(ch1, ch2,strf)

    #---------------------------------------------------------------------------      
    # Adds to the device a rewire gate
    #---------------------------------------------------------------------------      
    def rewire(self, ch1, ch2):
        """
            
        Adds a swap gate between two channels
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        
        """
        self.circ.rewire(ch1,ch2)
        self.disp.rewire(ch1, ch2)       

    #---------------------------------------------------------------------------                 
    # Adds to the device a phase shifter
    #---------------------------------------------------------------------------      
    def phase_shifter(self, ch, phi, colored=False):
        """
            
        Adds a phase shifter to the device in channel ch.
        
        :ch (int): Phase shifter input channel.
        :phi (float): Angle phi in degrees.
        :colored (bool): Prints the phase shifter in a different color than the default one.
        
        """
        str1="PS \u03C6="
        str2=str(round(float(phi),2))
        str3="º"
        strf=str1+str2+str3
        self.circ.phase_shifter(ch,phi)
        self.disp.one_gate(ch,strf,colored)

    #---------------------------------------------------------------------------          
    # Adds to the device an homogeneous lossy medium
    #---------------------------------------------------------------------------      
    def loss(self, ch, l):
        """
            
        Adds a lossy medium with loss probability l to the circuit in channel ch.
        
        :ch (int): Lossy medium input channel.
        :l  (float): Loss probability.
        
        """
        str1="LSS L="
        str2=str(round(l,2))
        strf=str1+str2
        self.circ.loss(ch,l)
        self.disp.one_gate(ch,strf)

    #---------------------------------------------------------------------------      
    #  Adds to the device a delay
    #---------------------------------------------------------------------------      
    def delay(self, ch):
        """
        
        Increases the optical path of a channel by a quantity equal to a period
        
        :ch (int):  Channel where the delay is introduced.
                
        """
        self.circ.delay(ch)
        self.disp.delay(ch)

    #---------------------------------------------------------------------------      
    # Adds to the device a rotator
    #---------------------------------------------------------------------------      
    def rotator(self, ch, theta, phi):
        """
        
        Adds a polarization rotation element to the device attached to channel ch.
        
        :ch (int):  Rotator input channel.
        :theta (double):  Angle theta in degrees.
        :phi   (double):  Angle phi in degrees.
                
        """
        str1="ROT \n \u03B8="
        str2=str(round(theta,2))
        str3="º\n \u03C6="
        str4=str(round(phi,2))
        str5="º"
        strf=str1+str2+str3+str4+str5        
        self.circ.rotator(ch,theta,phi)
        self.disp.one_gate(ch,strf)

    #---------------------------------------------------------------------------      
    # Adds to the device a polarized beamsplitter
    #---------------------------------------------------------------------------      
    def pol_beamsplitter(self, ch1, ch2, P, theta):
        """
            
        Adds a polarized beamsplitter attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        :P   (int): Polarization to which the beamsplitter is sensitive.
        :theta (float): Effectiveness of the beamsplitter. 90º => 50/50 beamsplitter. 0º => No sensitivity.
        
        """
        str1="  POLBS\n="
        if P==0:
            str2='H'
        else:
            str2='V'
        

        str1="POLBS("
        if P==0:
            str2='H'
        else:
            str2='V' 
               
        str3=")\n \u03B8="
        str4=str(round(theta,2))
        str5="º"
        strf=str1+str2+str3+str4+str5

        self.circ.pol_beamsplitter(ch1,ch2,P,theta)
        self.disp.two_gate(ch1, ch2,strf)

    #---------------------------------------------------------------------------                 
    # Adds to the device a polarized phase shifter
    #---------------------------------------------------------------------------      
    def pol_phase_shifter(self, ch, P, phi, colored=False):
        """
            
        Adds a polarized phase shifter to the device in channel ch.
        
        :ch (int): Phase shifter input channel.
        :P   (int): Polarization to which the phase shifter is sensitive.
        :phi (float): Angle phi in degrees.
        :colored (bool): Prints the phase shifter in a different color than the default one.
        
        """
                  
        str1="PS"
        if P==0:
            str2='(H)'
        else:
            str2='(V)'
        str3=" \u03C6="
        str4=str(round(float(phi),2))
        str5="º"
        strf=str1+str2+str3+str4+str5
        self.circ.pol_phase_shifter(ch,P,phi)
        self.disp.one_gate(ch,strf,colored)

    #---------------------------------------------------------------------------      
    # Adds to the device a polarization filter
    #---------------------------------------------------------------------------
    def pol_filter(self, ch, P):
        """
            
        Adds a polarization filter attached to channel ch.
        
        :ch  (int): Polarization filter input channel.
        :P   (int): Polarization to be filtered.
        
        """
                  
        str1="FILT"
        if P==0:
            str2='(H)'
        else:
            str2='(V)'
        strf=str1+str2
        self.circ.pol_filter(ch,P)
        self.disp.one_gate(ch,strf,False)
        
    #---------------------------------------------------------------------------      
    # Adds to the device a half waveplate
    #---------------------------------------------------------------------------      
    def half(self, ch, alpha):
        """
            
        Adds a half waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double): Rotation angle in degrees.
        
        """
        str1="Half \u03B1="
        str2=str(round(alpha,2))
        str3="º"
        strf=str1+str2+str3
        self.circ.half(ch,alpha)
        self.disp.one_gate(ch,strf)

    #---------------------------------------------------------------------------      
    # Adds to the device a quarter waveplate
    #---------------------------------------------------------------------------      
    def quarter(self, ch, alpha):
        """
            
        Adds a quarter waveplate attached to channel ch.
        
        :ch (int): Waveplate input channel.
        :alpha (double):  Rotation angle in degrees.
        
        """
        str1="Quart \u03B1="
        str2=str(round(alpha,2))
        str3="º"
        strf=str1+str2+str3
        self.circ.quarter(ch,alpha)
        self.disp.one_gate(ch,strf)        

    #---------------------------------------------------------------------------              
    # Adds to the device a detector
    #---------------------------------------------------------------------------      
    def detector(self, ch, cond=-1, pol=-1, mpi=-1, mpo=-1, eff=1.0, blnk=0.0, gamma=0.0):
        """
            
        Adds a detector to a channel of the device.
        
        :ch   (int): Detector channel.
        :cond (optional[int]): Detection condition. |br|
         cond>=0:  Readings in the remaining channels are considered only by calculations in probability bins and density matrices if the number of photons in this channel is equal to cond. |br|
         cond=-1:  There is no condition and works as a normal detector. |br|
         cond=-2:  The channel is ignored by outcome calculations in probability bins and density matrices. |br|
        :pol (optional[int]): Polarization condition. If cond>=0, pol determines the polarization of the photons to fulfill the condition. Note that if pol=-1 no assumption about the polarization of those photons is made.
        :mpi (optional[int]): Initial period of the detection window (if -1 takes the first one as default).
        :mpo (optional[int]): Final period of the detection window (if -1 takes the first one as default).
        :eff   (optional[float]): Efficiency of the detector.
        :blnk  (optional[float]): Ratio of time in which the detector is inactive due other detections.
        :gamma (optional[float]): Average rate of dark counts in this channel.
        
        """
        self.circ.detector(ch,cond,pol,mpi,mpo,eff,blnk,gamma)
        self.disp.final_dec(ch,cond,pol)
        self.outlist.append(ch)
        

    #---------------------------------------------------------------------------      
    # Ignore a channel
    #---------------------------------------------------------------------------      
    def ignore(self, ch):
        """
            
        Flags the channel to be ignored in the outcome calculations. No detector is placed in this channel.
        
        :ch (int): Ignored channel.
        
        """
        self.circ.detector(ch,-2,-1,-1,-1,1.0,0.0,0.0)
        self.disp.final_dec(ch,-2,-1)
        self.outlist.append(ch)

    #---------------------------------------------------------------------------              
    # Adds to the output a random Gaussian whote noise of stdev^2 dispersion
    #---------------------------------------------------------------------------      
    def noise(self, stdev2):
        """
            
        Adds Gaussian white noise to the output.
        
        :stdev2 (float): Dispersion of the Gaussian noise.
        
        """
        self.circ.noise(stdev2)

    #---------------------------------------------------------------------------      
    #  Apply a *single ket* post-selection condition defined by the detectors 
    # ( only valid for ideal circuits).
    #---------------------------------------------------------------------------         
    def apply_condition(self, inputst, ignore=True):
        """
        
        Apply the post-selection condition defined by the detectors in an ideal circuit to a state. |br|
        **Warning!** This can be only applied to ideal devices ns=1 where detector naturally define 
        a single ket projector. Otherwise the conditions lead to a density matrix output and there are 
        other tools available in SOQCS for that purpose. Note that ignored channels are not transcribed 
        but that may lead to collisions between otherwise different terms.           
    
        :inputst(state):  Input state to be post-selected applying the detector conditions.
        :ignore(optional[bool]): Consider ignored channels (False=No/True=Yes)
        :return (state): Post-selected state
    
        """
        return self.circ.apply_condition(inputst, ignore)

    #---------------------------------------------------------------------------              
    # Plots the device
    #---------------------------------------------------------------------------      
    def show(self,depth=10,sizexy=100,font=18, slin=0, rows=1):
        """
        
        Plots the quantum device.
        
        :depth (optional[int]): Number of plot layers. The longer the circuit the more needed.
        :sizexy(optional[int]): How many pixels by layer. Controls the size of the plot.
        :font (optional[int]):  Font size of the labels.
        :slin(optional(int)): Starting lines 0='No'/1='Yes'. If starting lines is 0 no lines are shown until the channel is initialized.
        :rows(int): Numbers of printing lines to be used by the plotter. |br| 
                   If rows=1, the plot is printed in one line unless hard separators are used to manually 
                   point out the change of printing line. |br| 
                   If rows>1, the plot is drawn in a fixed quantity of lines. In this case the change of line is done automatically. |br|

        
        """    
        self.disp.show(self.nch,depth,sizexy,font, slin, rows)
        
#------------------------------------------------------------------------------#                    
# Plots a probability function vs dt                                           #
#------------------------------------------------------------------------------#      
def plot(pfunc,cnvx,cnvy,textx,x0,x1,nxt,texty,y0,y1,nyt,N,argsv=[{0:0}],colorv=['b'],padb=0.0,padt=0.0):
    """
        
    Plots a multivalued function.
    
    :pfunc (method): Method with the function to be plotted pfunc(x,args) where x is a float and args is a dictionary of parameters.
    :cnvx (int): Size in x of the figure.
    :cnvy (int): Size in y of the figure.
    :textx (string): Label x axis.
    :x0 (float): Start x value.
    :x1 (float): End x value.
    :nxt (int): Number of x ticks.
    :texty (string): Label y axis.
    :y0 (float): Start y value.
    :y1 (float): End y value.
    :nyt (int): Number of y ticks.
    :N (int): Number of points.
    :argsv(dicc): Dictionary with the parameters for each sweep to pfunc.
    :colorv(list): Colors assigned to each sweep of pfunc.
    :padb (float): Pad bottom. Some extra space at the bottom to better visualization.
    :padt (float): Pad top. Some extra space at the top to better visualization.
        
    """
    Nv=len(argsv)
    Nc=len(colorv)
    if(Nv!=Nc):
        print("Plot warning!: Number of arguments different of number of colors to print")

    dt = np.linspace(x0, x1, N)    
    vpfunc=np.vectorize(pfunc)
    probs=[0]*Nv
    for i in range(0,Nv):
        probs[i]=vpfunc(dt,argsv[i])
    
    f=plt.figure(figsize=(cnvx,cnvy))
    ax1=f.add_subplot(111)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')

    plt.xlim(x0,x1)
    plt.ylim(y0-padb,y1+padt)
    tx=np.linspace(x0, x1, nxt)   
    ty=np.linspace(y0, y1, nyt)   
    plt.xticks(tx, size = 18)
    plt.yticks(ty, size = 18)
    plt.xlabel(textx, fontsize=18, labelpad=0)
    plt.ylabel(texty,  fontsize=18, labelpad=0)
    plt.subplots_adjust(left = 0.14)
    plt.subplots_adjust(bottom = 0.15)
    plt.tick_params(which='both', width=1)
    plt.tick_params(which='major', length=-7, pad=15)
    plt.tick_params(which='minor', length=-4, pad=15)
    for i in range(0,Nv):
        plt.plot(dt, probs[i],color=colorv[i], linewidth=3.0)
    plt.show()
    
#------------------------------------------------------------------------------#                    
# Scatter plot of (X,Y) pairs                                                 #
#------------------------------------------------------------------------------#      
def plot_data(probs,X,cnvx,cnvy,textx,x0,x1,nxt,texty,y0,y1,nyt,N,colorv=['b'],padb=0.0,padt=0.0):
    """
        
    Scatter plot of (X,Y) pairts  
    
    :Y (list[double]): Y-Data to be plotted of a (X,Y) pair
    :X (list[double]): X-Data to be plotted of a (X,Y) pair
    :cnvx (int): Size in x of the figure.
    :cnvy (int): Size in y of the figure.
    :textx (string): Label x axis.
    :x0 (float): Start x value.
    :x1 (float): End x value.
    :nxt (int): Number of x ticks.
    :texty (string): Label y axis.
    :y0 (float): Start y value.
    :y1 (float): End y value.
    :nyt (int): Number of y ticks.
    :N (int): Number of points.
    :color: Color of the plot.
    :padb (float): Pad bottom. Some extra space at the bottom to better visualization.
    :padt (float): Pad top. Some extra space at the top to better visualization.
        
    """

         
    Nv=len(probs)
    Nc=len(colorv)
    if(Nv!=Nc):
        print("Plot data warning!: Number of data arguments different of number of colors to print")
        
    f=plt.figure(figsize=(cnvx,cnvy))
    ax1=f.add_subplot(111)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')

    plt.xlim(x0,x1)
    plt.ylim(y0-padb,y1+padt)
    tx=np.linspace(x0, x1, nxt)   
    ty=np.linspace(y0, y1, nyt)   
    plt.xticks(tx, size = 18)
    plt.yticks(ty, size = 18)
    plt.xlabel(textx, fontsize=18, labelpad=0)
    plt.ylabel(texty,  fontsize=18, labelpad=0)
    plt.subplots_adjust(left = 0.14)
    plt.subplots_adjust(bottom = 0.15)
    plt.tick_params(which='both', width=1)
    plt.tick_params(which='major', length=-7, pad=15)
    plt.tick_params(which='minor', length=-4, pad=15)
    for i in range(0,Nv):
        plt.plot(X, probs[i],color=colorv[i], linewidth=3.0)
    plt.show()    