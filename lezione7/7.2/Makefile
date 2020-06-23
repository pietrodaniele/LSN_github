CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o MC_NVT.o
	$(CC) random.o MC_NVT.o funzioni.o main.o -o main.exe
random.o : random.cpp random.h
		$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o: funzioni.h funzioni.cpp
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
MC_NVT.o : MC_NVT.cpp MC_NVT.h
	$(CC) -c MC_NVT.cpp -o MC_NVT.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h MC_NVT.h
		$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe
