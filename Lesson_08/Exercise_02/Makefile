CC = g++
CFLAGS = -Wall -O3 --std=c++11

metropolis.exe : metropolis.o random.o histogram.o
	$(CC) random.o histogram.o metropolis.o -o metropolis.exe
metropolis.o : metropolis.cpp
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
histogram.o : histogram.cpp histogram.h
	$(CC) -c histogram.cpp -o histogram.o $(CFLAGS)
clean :
	rm *.o metropolis.exe seed.out
