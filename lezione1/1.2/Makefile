CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o FunzExp.o FunzLnz.o experiment.o 
	$(CC) experiment.o FunzLnz.o FunzExp.o random.o main.o -o main.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
FunzExp.o: FunzExp.h FunzExp.cpp FunzBase.h
	$(CC) -c FunzExp.cpp -o FunzExp.o $(CFLAGS)
FunzLnz.o: FunzLnz.h FunzLnz.cpp FunzBase.h
	$(CC) -c FunzLnz.cpp -o FunzLnz.o $(CFLAGS)
experiment.o : experiment.cpp experiment.h funzioni.h 
	$(CC) -c experiment.cpp -o experiment.o $(CFLAGS)
main.o : main.cpp funzioni.h random.h experiment.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out



