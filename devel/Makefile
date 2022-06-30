include ../conf.inc
mydir=./devel


main:	
	cd .$(LPATH);\
	make;\
	cd ..;\
	cd $(mydir);\
	$(GPP)  $(FLAGS)  main.cpp $(CFLAGS) $(LDFLAGS) -o main.x

	
cleanall:
	cd .$(LPATH);\
	make clean;\
	cd ..;\
	cd $(mydir);\
	rm *.x
	
clean:
	rm *.x
