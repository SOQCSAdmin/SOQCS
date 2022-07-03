"""@package pysoqcs
SOQCS Python port

This library is a wrapper that provides the functionality of C++ SOQCS library in Python 
"""

##########################################################################################################
#                                                                                                        #
#  SOQCS interface with Python. Python side.                                                             #
#                                                                                                        #
#                                                                                                        #
# Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.    #
# The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the    #
# root directory of this source tree. Use of the source code in this file is only permitted under the    #
# terms of the licence in LICENCE.TXT.                                                                   #
##########################################################################################################

import math 
import cmath 
import sys
from numpy import *
from ctypes import *
from ctypes import cdll,c_char,c_int,c_long,c_double,POINTER
import matplotlib.pyplot as plt
soqcs = cdll.LoadLibrary('./libSOQCS.so')


# Memory management auxiliary methods.
# Free the memory of an object created in C++
# No valid for classes
def free_ptr(array):
        array.restype=POINTER(c_char)
        soqcs.free_ptr(array)

# Converter MTX->int*
def to_int_ptr(mtx):
    n=len(mtx)
    m=len(mtx[0])
    length=n*m


    send=(c_int*length)()
    for i in range(0,n):
        for j in range(0,m):            
            send[i*m+j]=mtx[i][j]

    return send,n,m

# Wrapper for C++ SOQCS configuration method
# Configuration of the maximum number of photons of the simulation
def cfg_soqcs(nph):
    soqcs.all_cfg_soqcs(nph);   
    
# Wrapper for C++ SOQCS class qocircit
# Definition of a quantum optical circuit
class qocircuit(object):
    # Create circuit
    def __init__(self, *args):
        if len(args)==0:
            print("pyQocircuit : This object can not be created with 0 arguments")
            
        if len(args)==1:
            self.init1(args)
            
        if len(args)==2:
            print("pyQocircuit : This object can not be created with 2 arguments")

        if len(args)==3:
            self.init3(args)

        if len(args)==4:
            self.init4(args)

        if len(args)==5:
            print("pyQocircuit : This object can not be created with 5 arguments")

        if len(args)==6:
            self.init6(args)
            
        if len(args)>6:
            print("pyQocircuit : This object can not be created with more than 6 arguments")
            
    # Create circuit 1 argument        
    def init1(self,args):
        func=soqcs.qoc_new_qocircuit
        func.restype=c_long
        self.obj = func(args[0],1,1,0,0,False)
        
    # Create circuit 3 arguments        
    def init3(self,args):
        func=soqcs.qoc_new_qocircuit
        func.restype=c_long
        self.obj = func(args[0],args[1],args[2],0,0,False)

    # Create circuit 4 arguments        
    def init4(self,args):
        func=soqcs.qoc_new_qocircuit
        func.restype=c_long
        self.obj = func(args[0],args[1],args[2],args[3],0,False)

    # Create circuit 6 arguments        
    def init6(self,args):
        func=soqcs.qoc_new_qocircuit
        func.restype=c_long    
        self.obj = func(args[0],args[1],args[2],args[3],args[4],args[5])

    # Delete circuit        
    def __del__(self):
        soqcs.qoc_destroy_qocircuit(c_long(self.obj))

    #Return nunmber of levels          
    def num_levels(self):
        return soqcs.qoc_num_levels(c_long(self.obj))

    # Adds to the circuit a random circuit
    def random_circuit(self):
        soqcs.qoc_random_circuit(c_long(self.obj))

    # Adds to the circuit a NSX subcircuit
    def NSX(self, ch1, ch2, ch3):
        soqcs.qoc_NSX(c_long(self.obj),ch1,ch2,ch3)

    # Adds to the circuit an ideal beamsplitter
    def beamsplitter(self, ch1, ch2, theta, phi):
        soqcs.qoc_beamsplitter(c_long(self.obj),ch1,ch2,c_double(theta),c_double(phi))

    # Adds to the circuit an ideal beamsplitter
    def dielectric(self, ch1, ch2, t, r):
       soqcs.qoc_dielectric(c_long(self.obj),ch1,ch2,c_double(t.real),c_double(t.imag),c_double(r.real),c_double(r.imag))

    # Adds to the circuit a MMI 2x2 beamsplitter
    def MMI2(self, ch1, ch2):
        soqcs.qoc_MMI2(c_long(self.obj),ch1,ch2)
   
    # Adds to the circuit a general homogeneous medium
    def general_medium(self, ch, t):
        soqcs.qoc_phase_shifter(c_long(self.obj),ch,c_double(t.real),c_double(t.imag))
        
    # Adds to the circuit a phase shifter
    def phase_shifter(self, ch, phi):
        phir=math.pi*phi/180
        t=cmath.exp(1j*phir)
        self.general_medium(ch,t)
    
    # Adds to the circuit a lossy medium
    def loss(self, ch, l):
        self.general_medium(ch, math.sqrt(1.0-l))

     # Adds to the circuit a detector
    def detector(self, *args):
        if len(args)==0:
            print("pyDetector : This device can not be created with 0 arguments")
            
        if len(args)==1:
            self.detector1(args)

        if len(args)==2:
            self.detector2(args)
            
        if len(args)==3:
            print("pyDetector : This device can not be created with 3 arguments")

        if len(args)==4:
            print("pyDetector : This device can not be created with 4 arguments")

        if len(args)==5:
            self.detector5(args)
                        
        if len(args)>5:
            print("pyDetector : This device can not be created with more than 5 arguments")
            
            
    # Adds to the circuit a detector 1 argument
    def detector1(self, args):
        soqcs.qoc_detector(c_long(self.obj),args[0],-1,c_double(1.0),c_double(0.0),c_double(0.0))
        
    # Adds to the circuit a detector 2 arguments
    def detector2(self, args):
        soqcs.qoc_detector(c_long(self.obj),args[0],args[1],c_double(1.0),c_double(0.0),c_double(0.0))
        
    # Adds to the circuit a detector 5 arguments
    def detector5(self, args):
        soqcs.qoc_detector(c_long(self.obj),args[0],args[1],c_double(args[2]),c_double(args[3]),c_double(args[4]))

    # Adds to the output a random Gaussian whote noise of stdev^2 dispersion
    def noise(self, stdev2):
        soqcs.qoc_noise(c_long(self.obj),c_double(stdev2))
        
    # Adds to the circuit a packet definition
    def def_packet(self, n,t,f,w):
        soqcs.qoc_def_packet(c_long(self.obj),n,c_double(t),c_double(f),c_double(w))        

    # Calculates the emitted packet visibility
    def emitted_vis(self, i,j):
        soqcs.qoc_emitted_vis(c_long(self.obj),i,j)                

    # Adds to the circuit an emitter
    def emitter(self, ckind, rand):
        if(ckind=='G'): 
            ikind=0
        else:
            ikind=1
            
        soqcs.qoc_emitter(c_long(self.obj),ikind,rand)     
    
    #  Adds to the circuit a delay
    def delay(self, ch, dt):
        soqcs.qoc_delay(c_long(self.obj),ch,c_double(dt))     

    # Print circuit
    def prnt(self):
        soqcs.qoc_prnt(c_long(self.obj))
        sys.stdout.flush()

    # Set the print format
    def prnt_format(self, flag):
        soqcs.qoc_set_prnt_flag(c_long(self.obj),flag)


# Wrapper for C++ SOQCS class state
# Definition of a quantum photonic state
class state(object):
    
    # Create a state
    def __init__(self, *args):
        if len(args)==0:
            print("pyState : This object can not be created with 0 arguments")
            
        if len(args)==1:
            self.init1(args)

        if len(args)==2:
            self.init2(args)

        if len(args)>2:
            print("pyState : This object can not be created with more than 2 arguments")

            
    # Create state with a C++ object               
    def init1(self,args):
        self.obj = args[0]
        
    # Create a state
    def init2(self,args):
        func=soqcs.st_new_state
        func.restype=c_long
        self.obj = func(args[0],args[1])
        
    #Delete a state
    def __del__(self):
        soqcs.st_destroy_state(c_long(self.obj))

    # Add a new term (ket + amplitude) to a state
    def add_term(self, ampl, mtx,qoc):
            param=to_int_ptr(mtx)
            soqcs.st_add_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), param[0],param[1],param[2],c_long(qoc.obj)) 

    # Post-selection by a projector
    def post_selection(self, prj):           
            func=soqcs.st_post_selection
            func.restype=c_long
            obj=func(c_long(self.obj), c_long(prj.obj))
            newstate=state(obj) 
            return newstate

     # Prnt_state
    def prnt_state(self, qoc, loss, column):
            soqcs.st_prnt_state(c_long(self.obj), c_long(qoc.obj),loss,column)
            sys.stdout.flush()

    # Print bins
    def prnt_state(self, *args):
        if len(args)==0:
            self.prnt_state0()
            
        if len(args)==1:
            print("pyState : prnt state is not defined with 1 argument")

        if len(args)==2:
            self.prnt_state2(args)

        if len(args)==3:
            self.prnt_state3(args)
            
        if len(args)>3:
            print("pyState : prnt state is not defined with more than 3 argumenta")


     # Prnt_state 0 arguments
    def prnt_state0(self):
            soqcs.st_prnt_state(c_long(self.obj), c_long(0),False, 0)
            sys.stdout.flush()


     # Prnt_state 2 arguments
    def prnt_state2(self, args):
            qoc=args[0];
            qoc.restype=qocircuit
            soqcs.st_prnt_state(c_long(self.obj), c_long(qoc.obj),False,args[1])
            sys.stdout.flush()

     # Prnt_state 3 arguments
    def prnt_state3(self, args):
            qoc=args[0];
            qoc.restype=qocircuit
            soqcs.st_prnt_state(c_long(self.obj), c_long(qoc.obj),args[1], args[2])
            sys.stdout.flush()


# Wrapper for C++ SOQCS class projector
# Definition of a projector to perform post-selection
class projector(object):
   
    # Create a proejctor
    def __init__(self,nlevel,maxket):
        func=soqcs.prj_new_projector
        func.restype=c_long    
        self.obj = func(nlevel,maxket)
    
    # Delete a projector
    def __del__(self):
        soqcs.prj_destroy_projector(c_long(self.obj))
          
    # Add a new term (ket + amplitude) to a projector    
    def add_term(self, ampl, mtx,qoc):
            param=to_int_ptr(mtx)
            soqcs.prj_add_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), param[0],param[1],param[2],c_long(qoc.obj)) 


# Wrapper for C++ SOQCS class p_bin
# Definition of a set of probability bins
class p_bin(object):
    # Create a photon bunch
    def __init__(self, *args):
        if len(args)==0:
            print("pyPBin : This object can not be created with 0 arguments")
            
        if len(args)==1:
            self.init1(args)

        if len(args)==2:
            self.init2(args)
            
        if len(args)>2:
            print("pyPBin : This object can not be created with more than 2 arguments")

            
    # Create photon bunch with a C++ object        
    def init1(self,args):
        self.obj = args[0]
        
    # Create a state
    def init2(self,args):
        func=soqcs.pb_new_pbin
        func.restype=c_long
        self.obj = func(args[0],args[1])
        
    # Delete a p_bin
    def __del__(self):
        soqcs.pb_destroy_pbin(c_long(self.obj))
          
    # Add the statistics of a state to a pbin
    def add_state(self, state):
        soqcs.pb_add_state(c_long(self.obj),c_long(state.obj)) 

    # Calculate the measurement
    def calc_measure(self, qoc):
        soqcs.pb_calc_measure(c_long(self.obj),c_long(qoc.obj)) 

    # Return probability of the bin
    def nbins(self):
        return soqcs.pb_nbins(c_long(self.obj))
    
    #Return nunmber of levels          
    def num_levels(self):
        return soqcs.pb_num_levels(c_long(self.obj))

        # Return tag of the bin "bin name"
    def tag(self, index):
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

    # Return probability of the bin
    def prob(self, index):
        func=soqcs.pb_prob
        func.restype=c_double
        return func(c_long(self.obj),index)
    
    # Return probability of the bin from its definition
    def prob_def(self, mtx,qoc):
        func=soqcs.pb_prob_def
        func.restype=c_double
        param=to_int_ptr(mtx)
        return func(c_long(self.obj), param[0],param[1],param[2],c_long(qoc.obj))    

    # Print bins
    def prnt_bins(self, *args):
        if len(args)==0:
            self.prnt_bins0()
            
        if len(args)==1:
            print("pyPBin : prnt bins is not defined with 1 argument")

        if len(args)==2:
            self.prnt_bins2(args)

        if len(args)==3:
            self.prnt_bins3(args)
            
        if len(args)>3:
            print("pyPBin : prnt bins is not defined with more than 3 arguments")

    # Print bins 0 argument
    def prnt_bins0(self):        
        soqcs.pb_prnt_bins(c_long(self.obj),c_long(0), c_double(0.0),False)    
        sys.stdout.flush()
                
    # Print bins 2 argument
    def prnt_bins2(self, args):
        qoc=args[0];
        qoc.restype=qocircuit
        soqcs.pb_prnt_bins(c_long(self.obj),c_long(qoc.obj), c_double(args[1]),False)    
        sys.stdout.flush()
        
    # Print bins 3 arguments
    def prnt_bins3(self, args):
        qoc=args[0];
        qoc.restype=qocircuit
        soqcs.pb_prnt_bins(c_long(self.obj),c_long(qoc.obj), c_double(args[1]),args[2])    
        sys.stdout.flush()
    
    def show(self):
        taglist=list()
        problist=list()
        for i in range(0,self.nbins()):
            taglist.append(self.tag(i))
            problist.append(self.prob(i))
  
        zipped_lists = zip(taglist, problist)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        staglist, sproblist = [ list(tuple) for tuple in  tuples]
        
        fig, ax = plt.subplots(1, 1)
        y_pos = arange(len(staglist))
        ax.bar(y_pos, sproblist, align='center', alpha=0.5)
        ax.set_axisbelow(True)
        plt.grid(visible=True, which='major', axis='y', color='0.8', linestyle='--', linewidth=2)
        plt.xticks(y_pos, staglist)
        plt.ylabel('Frequency')
        plt.title('Output')
        
        plt.show()

# Wrapper for C++ SOQCS class ph_bunch
# Definition of a set of photon bunches
class ph_bunch(object):
    # Create a photon bunch
    def __init__(self, *args):
        if len(args)==0:
            print("pyPHBunch : This object can not be created with 0 arguments")

        if len(args)==1:
            self.init1(args)

        if len(args)==2:
            self.init2(args)
            
        if len(args)>2:
            print("pyPHBunch : This object can not be created with more than 2 arguments")

            
    # Create photoon bunch with a C++ object        
    def init1(self,args):
        self.obj = args[0]
        
    # Create a state
    def init2(self,args):
        func=soqcs.ph_new_bunch
        func.restype=c_long
        self.obj = func(args[0],args[1])
        
    # Delete a photon bunch
    def __del__(self):
        soqcs.ph_destroy_bunch(c_long(self.obj))

    # Add photons to the photon bunch
    def add_photons(self, *args):
        if len(args)<3:
            print("pyPHBunch : Whe can not add photons with less than 3 arguments")
            
        if len(args)==3:
            self.add_photons3(args)
            
        if len(args)>3 and len(args)<7:
            print("pyPHBunch : Whe can not add photons with", len(args),"arguments")
            print("It has to be either 3 or 7 arguments")

        if len(args)==7:
            self.add_photons7(args)
            
        if len(args)>7:
            print("pyPHBunch : Whe can not add photons with more than 7 arguments")
             
            
    # Add photons to the photon bunch 3 arguments
    def add_photons3(self, args):
        qoc=args[2];
        qoc.restype=qocircuit
        soqcs.ph_add_photons(c_long(self.obj),args[0],args[1],0,c_double(0.0),c_double(0.0),c_double(0.0),c_long(qoc.obj)) 
          
    # Add photons to the photon bunch
    def add_photons7(self, args):
        qoc=args[6];
        qoc.restype=qocircuit
        soqcs.ph_add_photons(c_long(self.obj),args[0],args[1],args[2],c_double(args[3]),c_double(args[4]),c_double(args[5]),c_long(qoc.obj)) 

    # Send the photons to the circuit
    def send2circuit(self, ckind,rand,qoc):
        if(ckind=='G'): 
            ikind=0
        else:
            ikind=1
        soqcs.ph_send2circuit(c_long(self.obj),ikind,rand,c_long(qoc.obj)) 
    
    # Weight of the photons particular configuration (amplitude of the state)
    def weight(self, A):
        soqcs.ph_weight(c_long(self.obj), c_double(A.real),c_double(A.imag))
    
    # Consider an alternative. ( Start a new ket )
    def alternative(self):
        soqcs.ph_alternative(c_long(self.obj))

    # Return the state of the bunch of photons
    def bunch_state(self):  ## Look at memory management here too. Maybe easier to put flag and do not create
        func=soqcs.ph_bunch_state
        func.restype=c_long
        obj=func(c_long(self.obj))
        newstate=state(obj)
        return newstate        

    # Return print the state of the bunch of photons
    def prnt_state(self, qoc,loss,column):
        soqcs.ph_prnt_state(c_long(self.obj),c_long(qoc.obj),loss,column)    
        sys.stdout.flush()

        
# Wrapper for C++ SOQCS class ph_bunch
# Definition of a set of photon bunches        
class simulator(object):
    # Create a simulator
    def __init__(self, *args):
        if len(args)==0:
            self.init0()

        if len(args)==1:
            print("pySimulator : This object can not be created with only 1 argument")

        if len(args)==2:
            self.init2(args)
            
        if len(args)>2:
            print("pySimulator : This object can not be created with more than 2 arguments")

            
    # Create a simulator 0 arguments
    def init0(self):
        func=soqcs.sim_new_simulator
        func.restype=c_long    
        self.obj = func(0,1000)
        
    # Create a simulator 2 arguments
    def init2(self,args):
        func=soqcs.sim_new_simulator
        func.restype=c_long    
        self.obj = func(args[0],args[1])

    # Delete a simulator
    def __del__(self):
        soqcs.sim_destroy_simulator(c_long(self.obj))
          
    # Run a simulator for a buch of photons
    def run(self, photons,qoc):
        func=soqcs.sim_run
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(photons.obj),c_long(qoc.obj)) 
        newoutcome=p_bin(obj)
        return newoutcome

    # Run a simulator for a input state
    def run_state(self, istate,qoc):
        func=soqcs.sim_run_state
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj)) 
        newstate=state(obj)
        return newstate

    # Run a simulator for a buch of photons
    def sample(self, photons,qoc,N):
        func=soqcs.sim_sample
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(photons.obj),c_long(qoc.obj),N) 
        newoutcome=p_bin(obj)
        return newoutcome

    # Run a simulator for a input state
    def sample_state(self, istate,qoc,N):
        func=soqcs.sim_sample_state
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj),N) 
        newstate=state(obj)
        return newstate