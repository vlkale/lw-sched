#OPTS  = -O3 -lpthread -lrt -lm -openmp
OPTS  = -O3 -lpthread -lrt -lm
MPICC = mpicc
MPICXX = mpicxx
CC  = gcc $(OPTS)
GCC = gcc $(OPTS)
GXX = c++ $(OPTS)
CLANGXX = clang++
DEBUG=-g

# Choose the number of cores on the system for make test
NUMCORES=64 

all: test_vSched testAppTwo_omp-lols-vSched testAppTwo_omp-lols-uds testAppOne_omp-lols-vSched testAppOne_omp-lols-uds

test_vSched: appFor_vSchedSimple.c vSched.h vSched.c
	$(CLANGXX) -fPIC vSched.c appFor_vSchedSimple.c -DCDY_ $(OPTS) -o test_vSched

testAppOne_omp-lols-vSched: appFor_omp-lols.c vSched.h vSched.c
	$(CLANGXX) appFor_omp-lols.c vSched.c -fopenmp -DUSE_VSCHED $(OPTS) -o testAppOne_omp-lols-vSched

testAppOne_omp-lols-uds: appTwoFor_omp-lols.c vSched.h vSched.c
	$(CLANGXX) appFor_omp-lols.c vSched.c -fopenmp $(OPTS) -o testAppOne_omp-lols-uds

testAppTwo_omp-lols-vSched: appTwoFor_omp-lols.c vSched.h vSched.c
	$(CLANGXX) appTwoFor_omp-lols.c vSched.c -fopenmp -DUSE_VSCHED $(OPTS) -o testAppTwo_omp-lols-vSched

testAppTwo_omp-lols-uds: appTwoFor_omp-lols.c vSched.h vSched.c
	$(CLANGXX) appTwoFor_omp-lols.c vSched.c -fopenmp $(OPTS) -o testAppTwo_omp-lols-uds

test_tQueue: appFor_vSchedSimple.C
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_threadedQueueSimple.C -DCDY_ $(OPTS) -o test_vSched

test_vSchedforMac: appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c
	$(GXX) -fPIC -g vSched.h vSched.c appFor_vSchedSimple.c pthBarrierforOSX.c -DCDY_ $(OPTS) -o test_vSchedforMac

test_tqueue_forMac: threadedQueue.h threadedQueue.C pthBarrierforOSX.c
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_tqueueSimple.c pthBarrierforOSX.c $(OPTS) -o test_tqueue_forMac

test_vSchedOpenMP: appFor_vSchedSimpleOpenMP.c vSched.h vSched.c
	$(GXX) -fPIC -g -I. vSched.h vSched.c appFor_vSchedSimpleOpenMP.c $(OPTS) -o test_vSchedOpenMP

test_vSchedomp: appFor_vSched-omp.C vSched.h vSched.c
	$(MPICXX) -fPIC -g $(OPTS) -I. vSched.h vSched.c appFor_vSched-omp.C -o test_vSchedomp

test:
	./test_vSched 65536 10 $(NUMCORES) 64 0.5 0.1
	./testAppOne_omp-lols-vSched 65536 10 $(NUMCORES) 64 0.5 0.1
	./testAppOne_omp-lols-uds 65536 10 $(NUMCORES) 64 0.5 0.1
#	./testAppTwo_omp-lols-vSched 65536 10 $(NUMCORES) 64 0.5 0.1
#	./testAppTwo_omp-lols-uds 65536 10 $(NUMCORES) 64 0.5 0.1
	./testAppTwo_omp-lols-vSched 500 1 64
	./testAppTwo_omp-lols-uds 500 1 64

tgz: appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c 
	tar -cvzf appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c README

clean:
	rm -rf *.o test_vSched test_tqueue_forMac test_vSchedomp  testAppTwo_omp-lols-vSched testAppTwo_omp-lols-uds testAppOne_omp-lols-vSched testAppOne_omp-lols-uds

realclean:
	rm -rf *.o core *.gch test_vSched test_vSchedforMac test_vSchedOpenMP test_vSchedomp  testAppTwo_omp-lols-vSched testAppTwo_omp-lols-uds testAppOne_omp-lols-vSched testAppOne_omp-lols-uds
