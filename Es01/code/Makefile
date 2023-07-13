CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o utils.hpp
	$(CC) Random/random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : Random/random.cpp Random/random.h
	$(CC) -c Random/random.cpp -o Random/random.o $(CFLAGS)
clean :
	rm *.o main.exe Random/*.o Random/seed.out
