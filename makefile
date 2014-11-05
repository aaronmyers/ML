CC= icc
OBJS= HW2p2.cpp
INPUT= 20 1 16 ./data/large/ 

all: Runit

omp_als: $(OBJS)
	$(CC) -openmp -o omp_als $(OBJS)

Runit: omp_als
	./omp_als $(INPUT)

