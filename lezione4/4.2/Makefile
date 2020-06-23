CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o experiment.o
	$(CC) experiment.o MolDyn_NVE.o -o MolDyn_NVE.exe
experiment.o : experiment.cpp experiment.h 
	$(CC) -c experiment.cpp -o experiment.o $(CFLAGS)
main.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o MolDyn_NVE.exe 
