CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

esercitazione.exe : esercitazione.o random.o
	$(CC) random.o esercitazione.o -o esercitazione.exe
esercitazione.o : esercitazione.cpp
	$(CC) -c esercitazione.cpp -o esercitazione.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercitazione.exe
