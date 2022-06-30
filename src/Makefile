include ../conf.inc

libsocs.a:   dmat.o mthread.o sim.o state.o qocircuit.o util.o
	     $(AR) $(ARFLAGS) libsoqcs.a mthread.o dmat.o sim.o state.o qocircuit.o util.o

dmat.o:      dmat.cpp
	     $(GPP)  $(FLAGS) -c dmat.cpp

mthread.o:   mthread.cpp
	     $(GPP)  $(FLAGS) -c mthread.cpp

sim.o:       sim.cpp 
	     $(GPP)  $(FLAGS) -c sim.cpp  

state.o:     state.cpp
	     $(GPP)  $(FLAGS) -c state.cpp

qocircuit.o: qocircuit.cpp
	     $(GPP)  $(FLAGS) -c qocircuit.cpp

util.o:      util.cpp
	     $(GPP)  $(FLAGS) -c util.cpp


clean:
	     rm *.o *.a
