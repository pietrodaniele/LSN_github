CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funz_try.o metropolis.o funzioni.o
	$(CC) funzioni.o metropolis.o  funz_try.o random.o main.o -o main.exe
funzioni.o: funzioni.h funzioni.cpp
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funz_try.o: funz_try.h funz_try.cpp
	$(CC) -c funz_try.cpp -o funz_try.o $(CFLAGS)
metropolis.o: metropolis.h metropolis.cpp
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h metropolis.h FunzBase.h funz_try.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
