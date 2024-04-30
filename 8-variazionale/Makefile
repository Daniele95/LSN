CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercitazione.exe : esercitazione.o random.o ../utils/utils.o
	$(CC) random.o esercitazione.o ../utils/utils.o -o esercitazione.exe
esercitazione.o : esercitazione.cpp
	$(CC) -c esercitazione.cpp -o esercitazione.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)


utils.o : ../utils/utils.cpp ../utils/utils.h 
	$(CC) -c ../utils/utils.cpp -o ../utils/utils.o $(CFLAGS)	
	
clean :
	rm *.o esercitazione.exe 
