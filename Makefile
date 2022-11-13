include ./conf.inc


all:	
	cd $(PPATH);\
	make;\
	cd ..;\
	cd $(EPATH);\
	make;\
	cd ..;  
	     
clean:
	cd $(PPATH);\
	make clean;\
	cd ..;\
	cd $(EPATH);\
	make clean;\
	cd ..;
