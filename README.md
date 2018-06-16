
---
# Description of the Library

This library allows users of an OpenMP program to use different loop scheduling strategies than the provided loop schedules. The library that is built from the code files in this directory
lets a user insert annotations for OpenMP loops in an application program to use one of the several efficient loop schedules provided in this library. The existing set of loop schedules in this library can be augmented through open-source code development. The schedules are intended to provides an effective tradeoff among factors including dynamic load imbalance, synchronization overhead and cache misses due to loss of locality. Specifically, the code of this loop scheduling library contains the key strategies for a static/dynamic scheduling scheme for MPI+OpenMP programs. Since Oct 2017, the library has served to guide the addition of a user-defined schedule to an upcoming specification of OpenMP (see more information about the user-defined schedule at bit.ly/udsomp). This is the library for lightweight scheduling which includes a constraint for reduction of cache misses and locking overhead.

The code file vSched.c contains the functions for initializing the scheduler, loop_start_<strat>, and the dequeueing functions
(loop_next_<strat>), where <strat> is the particular strategy. The implementor can add other strategies, by looking at the 
example functions in vSched.c shown.

There are two example application code files included in this folder: one is using pthreads, called appFor_vSched.c, and the other is using OpenMP, named appFor_vSchedOpenMP.c. To use this library's functionality in those application programs or to your MPI+OpenMP application program, you need to link this library with a library that predicts MPI communication time for an upcoming collective call or MPI_Isend() / MPI_Irecv() / MPI_Waitall() during execution of an application program. The library is in the repository named 'slack-trace' on my github page at https://github.com/vivek224.

The code of this loop scheduling library contains the key idea for the within-node solution proposed and studied in my dissertation. See my dissertation by visiting my googlepages webpage https://vivek112.googlepages.com, clicking on the link Research at the top of the page, scrolling to the bottom of the page where you see a list of publications, and clicking on the publication titled 'Low-overhead Scheduling to Improve Performance of Scientific Applications'. Other relevant papers can be found in that list by searching for the keyword 'low-overhead', 'locality-sensitive' or 'load balancing' in that list.

---

# Installing and Configuring the Library

The code doesn't have a configure file and I've just left out the configure file for now. I'm modifying the configure file to ensure its correctness before adding it to this open-source version. I intend to add it in the next month. I believe the
pthBarrierOSX.h should work for any Mac OS version. If you'd like to run on a Mac, you'll be using the pthBarrierOSX.h file when running the code.


---

# Compiling and Using the Library

1. Pthreads example
 
    a. To compile the pthreads example, type:
 
      make clean; make

    b. To run the pthreads example on Mac, you may need the pthBarrierOSX.c file and the pthBarrierOSX.h file, included in this directory. This is taken care of in the Makefile, as below:

      make clean; make test_vSchedOnMac; 
  
  
Include the pthBarrierforOSX.h in the appFor_vSchedSimple.c file.


2. To run the OpenMP example:


   set OMP_NUM_THREADS={number_of_cores_on_node} make clean; make test_vSched_OMP; 
   

   This can also be done also using my library.

---
