OPTS  = -O3 -lpthread -lrt -lm
MPICC   = mpicc $(OPTS)
CC  = gcc $(OPTS)
GCC = gcc $(OPTS)
GXX = c++ $(OPTS)

all: test_vSched

test_vSched: appFor_vSchedSimple.c vSched.h vSched.c
	$(GXX) -fPIC -g vSched.h vSched.c appFor_vSchedSimple.c -DCDY_ $(OPTS) -o test_vSched

test_vSchedforMac: appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c
	$(GXX) -fPIC -g vSched.h vSched.c appFor_vSchedSimple.c pthBarrierforOSX.c -DCDY_ $(OPTS) -o test_vSchedforMac

test_vSchedOpenMP: appFor_vSchedSimpleOpenMP.c vSched.h vSched.c
	$(GXX) -fPIC -g -I. vSched.h vSched.c appFor_vSchedSimpleOpenMP.c $(OPTS) -fopenmp -o test_vSchedOpenMP

tgz: appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c 
	tar -cvzf appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c README

clean:
	rm -rf *.o test_vSched test_vSchedforMac test_vSchedOpenMP

realclean:
	rm -rf *.o core *.gch test_vSched test_vSchedforMac test_vSchedOpenMP