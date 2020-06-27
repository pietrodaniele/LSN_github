CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : random.o finance.o funzioni.o experiment.o main.o
	$(CC) experiment.o finance.o funzioni.o random.o main.o -o main.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
experiment.o : experiment.cpp experiment.h funzioni.h
	$(CC) -c experiment.cpp -o experiment.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
finance.o : finance.cpp finance.h funzioni.h
		$(CC) -c finance.cpp -o finance.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h experiment.h finance.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
