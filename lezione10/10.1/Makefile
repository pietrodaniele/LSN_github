CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o SA.o Individual.o funzioni.o
	$(CC) funzioni.o Individual.o  SA.o random.o main.o -o main.exe
funzioni.o: funzioni.h funzioni.cpp
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
SA.o: SA.h SA.cpp
	$(CC) -c SA.cpp -o SA.o $(CFLAGS)
Individual.o: Individual.h Individual.cpp
	$(CC) -c Individual.cpp -o Individual.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h Individual.h SA.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
