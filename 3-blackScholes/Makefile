CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : ex1.o random.o
	$(CC) random.o ex1.o -o BlackScholes.exe
ex1.o : BlackScholes.cpp
	$(CC) -c BlackScholes.cpp -o ex1.o $(CFLAGS)
clean :
	rm *.o ../genRandom/ex1.exe seed.out
