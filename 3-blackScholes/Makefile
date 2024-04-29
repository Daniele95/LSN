CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : ex1.o random.o ../utils/utils.o
	$(CC) random.o ex1.o ../utils/utils.o -o BlackScholes.exe
ex1.o : BlackScholes.cpp
	$(CC) -c BlackScholes.cpp -o ex1.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)

utils.o : ../utils/utils.cpp ../utils/utils.h 
	$(CC) -c ../utils/utils.cpp -o ../utils/utils.o $(CFLAGS)
clean :
	rm *.o ../random/ex1.exe seed.out
