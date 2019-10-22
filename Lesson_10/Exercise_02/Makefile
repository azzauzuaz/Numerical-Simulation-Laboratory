MPI = mpic++
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(MPI) random.o main.o -o main.exe
main.o : main.cpp functions.h variables.h
	$(MPI) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(MPI) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
