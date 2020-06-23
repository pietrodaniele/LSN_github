CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o random_walk.o 
	$(CC) random_walk.o random.o main.o -o main.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
experiment.o : random_walk.cpp random_walk.h funzioni.h 
	$(CC) -c random_walk.cpp -o random_walk.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h random_walk.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
