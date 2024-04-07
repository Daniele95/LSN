CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : ex1.o random.o
	$(CC) random.o ex1.o -o ex1.exe
main.o : ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : ../genRandom/random.cpp ../genRandom/random.h
	$(CC) -c ../genRandom/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ../genRandom/ex1.exe seed.out
