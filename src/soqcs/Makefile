include ./conf.inc

libpysoqcs.so: interface libsoqcs.a
	$(GPP) $(FLAGS) -shared  pysoqcs.o libsoqcs.a -o libpysoqcs.so

interface:
	$(GPP)  $(FLAGS) -c pysoqcs.cpp -o pysoqcs.o	

libsoqcs.a:  mthread.o sim.o dmat.o pbin.o qodev.o state.o qocircuit.o util.o
	     $(AR) $(ARFLAGS) libsoqcs.a mthread.o sim.o dmat.o pbin.o qodev.o state.o qocircuit.o util.o
     
mthread.o:   mthread.cpp
	     $(GPP)  $(FLAGS) -c mthread.cpp

sim.o:       sim.cpp 
	     $(GPP)  $(FLAGS) -c sim.cpp  

dmat.o:      dmat.cpp
	     $(GPP)  $(FLAGS) -c dmat.cpp

pbin.o:      pbin.cpp
	     $(GPP)  $(FLAGS) -c pbin.cpp
	     
qodev.o:     qodev.cpp
	     $(GPP)  $(FLAGS) -c qodev.cpp
	     
state.o:     state.cpp
	     $(GPP)  $(FLAGS) -c state.cpp

qocircuit.o: qocircuit.cpp
	     $(GPP)  $(FLAGS) -c qocircuit.cpp

util.o:      util.cpp
	     $(GPP)  $(FLAGS) -c util.cpp


clean:
	     rm *.so *.o *.a
