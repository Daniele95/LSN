CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
Monte_Carlo_NVT.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
	
clean :
	rm *.o MolDyn_NVE.exe
