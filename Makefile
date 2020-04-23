#OPTS  = -O3 -lpthread -lrt -lm -openmp
OPTS  = -O3 -lpthread -lrt -lm
MPICC = mpicc
MPICXX = mpicxx
CC  = gcc $(OPTS)
GCC = gcc $(OPTS)
GXX = c++ $(OPTS)
CLANGXX = clang++
DEBUG=-g

OPENMPOPT=-fopenmp

# Choose the number of cores on the system for make test
NUMCORES=64

INCLUDE_DIR = include
SRC_DIR = src
LIB_DIR = lib
BIN_DIR = bin
EXAMPLES_DIR = examples
VVTESTS_DIR = tests/vv
PERFTESTS_DIR = tests/perf

# shouldn't be needed here
UKERNELS_DIR = share/ukernels

#vSched and uds could be thought of as two different vendor implementations, though they actually have different interfaces in this case as well (the technique is the same though)

all: testAppTwo_omp-lols-vSched testAppTwo_omp-lols-uds testAppOne_omp-lols-vSched testAppOne_omp-lols-uds

#tODO: need to figure out how to get UKERNELS dir to compile

test_vSched: $(PERFTESTS_DIR)/testOneFor_pthread-lols.c 
	$(CLANGXX) -fPIC $(SRC_DIR)/vSched.c -I/$(INCLUDE_DIR) $(PERFTESTS_DIR)/testOneFor_pthread-lols.c -DCDY_ $(OPTS) -o $(BIN_DIR)/test_vSched

testAppOne_omp-lols-vSched: $(PERFTESTS_DIR)/testOneFor_omp-lols.C
	$(CLANGXX) -fPIC $(SRC_DIR)/vSched.c $(PERFTESTS_DIR)/testOneFor_omp-lols.C $(OPENMPOPT) -DUSE_VSCHED $(OPTS) -o $(BIN_DIR)/testAppOne_omp_lols-vSched

testAppOne_omp-lols-vSched-Mac: $(PERFTESTS_DIR)/testOneFor_omp-lols.C
	$(CLANGXX) -fPIC $(SRC_DIR)/vSched.c $(PERFTESTS_DIR)/testOneFor_omp-lols.C $(OPENMPOPT) -DUSE_VSCHED pthBarrierforOSX.c $(OPTS) -o $(BIN_DIR)/testAppOne_omp_lols-vSched-Mac

testAppOne_omp-lols-uds: $(PERFTESTS_DIR)/testOneFor_omp-lols.C
	$(CLANGXX) -fPIC $(PERFTESTS_DIR)/testOneFor_omp-lols.C $(OPTS) $(OPENMPOPT) -o $(BIN_DIR)/testAppOne_omp_lols-uds

testAppTwo_omp-lols-vSched: $(PERFTESTS_DIR)/testTwoFor_omp-lols.c
	$(CLANGXX) -fPIC $(SRC_DIR)/vSched.c $(PERFTESTS_DIR)/testTwoFor_omp-lols.C -DUSE_VSCHED $(OPENMPOPT) $(OPTS) -o $(BIN_DIR)/testAppTwo_omp_lols-vSched

testAppTwo_omp-lols-uds: $(PERFTESTS_DIR)/testTwoFor_omp-lols.c
	$(CLANGXX) $(PERFTESTS_DIR)/testTwoFor_omp-lols.c $(OPENMPOPT) $(OPTS) -o $(BIN_DIR)/testAppTwo_omp-lols-uds

test_tQueue: appFor_vSchedSimple.C
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_threadedQueueSimple.C -DCDY_ $(OPTS) -o test_vSched

test_vSchedforMac: appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c
	$(GXX) -fPIC -g vSched.h vSched.c appFor_vSchedSimple.c pthBarrierforOSX.c -DCDY_ $(OPTS) -o test_vSchedforMac

test_tqueue_forMac: threadedQueue.h threadedQueue.C pthBarrierforOSX.c
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_tqueueSimple.c pthBarrierforOSX.c $(OPTS) -o test_tqueue_forMac

test_vSchedOpenMP: appFor_vSchedSimpleOpenMP.c vSched.h vSched.c
	$(GXX) -fPIC -g -I. vSched.h vSched.c appFor_vSchedSimpleOpenMP.c $(OPENMPOPT) $(OPTS) -o test_vSchedOpenMP

test_vSchedomp: appFor_vSched-omp.C vSched.h vSched.c
	$(MPICXX) -fPIC -g $(OPTS) -I. vSched.h vSched.c appFor_vSched-omp.C $(OPENMPOPT) -o test_vSchedomp

test:
	$(BIN_DIR)/test_vSched 65536 10 $(NUMCORES) 64 0.5 0.1
	$(BIN_DIR)/testAppOne_omp-lols-vSched 65536 10 $(NUMCORES) 64 0.5 0.1
	$(BIN_DIR)/testAppOne_omp-lols-uds 65536 10 $(NUMCORES) 64 0.5 0.1
#	./testAppTwo_omp-lols-vSched 65536 10 $(NUMCORES) 64 0.5 0.1
#	./testAppTwo_omp-lols-uds 65536 10 $(NUMCORES) 64 0.5 0.1
	$(BIN_DIR)/testAppTwo_omp-lols-vSched 500 1 64
	$(BIN_DIR)/testAppTwo_omp-lols-uds 500 1 64

tgz: appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c 
	tar -cvzf appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c README

clean:
	rm -rf *.o $(BIN_DIR)/test_v* $(BIN_DIR)/test_o*

realclean:
	rm -rf *.o core *.gch $(BIN_DIR)/test_v* $(BIN_DIR)/test_o* 
