CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o ising_1D.o
	$(CC) ising_1D.o funzioni.o random.o main.o -o main.exe
funzioni.o: funzioni.h funzioni.cpp 
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
ising_1D.o: ising_1D.h ising_1D.cpp
	$(CC) -c ising_1D.cpp -o ising_1D.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h ising_1D.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out insta_values.0
