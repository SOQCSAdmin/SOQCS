include ../conf.inc

libsocs.so: interface soqcs
	$(GPP) $(FLAGS) -shared  pySOQCS.o ./src/libsoqcs.a -o libSOQCS.so

interface:
	$(GPP)  $(FLAGS) -c pySOQCS.cpp -o pySOQCS.o	
	
soqcs:
	cd src && make

clean:
	cd src && make clean
	rm *.so *.o	     
