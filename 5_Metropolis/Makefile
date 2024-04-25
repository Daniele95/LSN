CC = g++
CFLAGS = -Wall -O3 --std=c++11

metropolis.exe : metropolis.o ../random/random.o ../utils/utils.o 
	$(CC) ../random/random.o ../utils/utils.o  metropolis.o -o metropolis.exe

metropolis.o : metropolis.cpp
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)

random.o : ../random/random.cpp ../random/random.h 
	$(CC) -c ../random/random.cpp -o ../random/random.o $(CFLAGS)

utils.o : ../utils/utils.cpp ../utils/utils.h 
	$(CC) -c ../utils/utils.cpp -o ../utils/utils.o $(CFLAGS)

clean : 
	rm *.o metropolis.exe
