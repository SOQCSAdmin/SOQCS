# include ../conf.inc

libsocs.so: interface soqcs
	g++ -shared  pySOQCS.o ./src/libsoqcs.a -o libSOQCS.so

interface:
	g++ -c -fPIC pySOQCS.cpp -o pySOQCS.o	
	
soqcs:
	cd src && make

clean:
	cd src && make clean
	rm *.so *.o	     
