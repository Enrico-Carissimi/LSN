CC = g++
CFLAGS = -g -Wall -O3 --std=c++11

main.exe : main.o Random/random.o utils.o
	$(CC) Random/random.o utils.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : Random/random.cpp Random/random.h
	$(CC) -c Random/random.cpp -o Random/random.o $(CFLAGS)
utils.o : utils.cpp utils.hpp
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
clean :
	rm *.o Random/*.o main.exe
