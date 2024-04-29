CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : randomWalk.o random.o ../utils/utils.o
	$(CC) random.o randomWalk.o ../utils/utils.o -o randomWalk.exe
main.o : randomWalk.cpp
	$(CC) -c randomWalk.cpp -o randomWalk.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
	

utils.o : ../utils/utils.cpp ../utils/utils.h 
	$(CC) -c ../utils/utils.cpp -o ../utils/utils.o $(CFLAGS)	
clean :
	rm *.o randomWalk.exe seed.out
