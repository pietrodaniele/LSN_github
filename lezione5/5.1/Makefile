CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o hydro_100.o hydro_210.o metropolis.o funzioni.o
	$(CC) funzioni.o metropolis.o hydro_210.o hydro_100.o random.o main.o -o main.exe
funzioni.o: funzioni.h funzioni.cpp
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
hydro_100.o: hydro_100.h hydro_100.cpp
	$(CC) -c hydro_100.cpp -o hydro_100.o $(CFLAGS)
hydro_210.o: hydro_210.h hydro_210.cpp
	$(CC) -c hydro_210.cpp -o hydro_210.o $(CFLAGS)
metropolis.o: metropolis.h metropolis.cpp
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h metropolis.h FunzBase3D.h hydro_100.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
