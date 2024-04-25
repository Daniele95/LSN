CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o ../random/random.o ../utils/utils.o
	$(CC) ../random/random.o main.o ../utils/utils.o -o main.exe
	
	
main.o : variazionale.cpp
	$(CC) -c variazionale.cpp -o main.o $(CFLAGS)
	
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o ../random/random.o $(CFLAGS)
	
utils.o : ../utils/utils.cpp ../utils/utils.h 
	$(CC) -c ../utils/utils.cpp -o ../utils/utils.o $(CFLAGS)
	
clean :
	rm *.o main.exe
