CC = mpicxx
CFLAGS = -g -Wall -O3 --std=c++11

TSP.exe : TSP.o Random/random.o utils.o Individual.o Population.o
	$(CC) Random/random.o utils.o Individual.o Population.o TSP.o -o TSP.exe
TSP.o : TSP.cpp
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
random.o : Random/random.cpp Random/random.h
	$(CC) -c Random/random.cpp -o Random/random.o $(CFLAGS)
utils.o : utils.cpp utils.hpp
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
Individual.o : Individual.cpp Individual.hpp
	$(CC) -c Individual.cpp -o Individual.o $(CFLAGS)
Population.o : Population.cpp Population.hpp
	$(CC) -c Population.cpp -o Population.o $(CFLAGS)
clean :
	rm *.o Random/*.o TSP.exe