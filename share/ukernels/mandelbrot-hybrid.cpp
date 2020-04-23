#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

#ifdef USE_MPI
#include <mpi.h>
#endif 

#define MAX_ITERATIONS 1000
#define WIDTH 1536
#define HEIGHT 1024
#define BLOCK_HEIGHT 4
#define BLOCK_WIDTH 1536
#define NUM_JOBS 256

double When()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void fileWrite(int* pixels) {
	
	int i, j;
	FILE *fp;
	fp = fopen("test.pgm", "w");
	
	if (fp == NULL) {
		perror ( "Unable to open file" );
		exit (EXIT_FAILURE);
	}
	
	fprintf(fp,"P2\n");
	fprintf(fp,"%d %d\n",WIDTH,HEIGHT);
	fprintf(fp,"%d\n",MAX_ITERATIONS);
	
	
	for (j = 0; j < HEIGHT; j++) {
		for (i = 0; i < WIDTH; i++) {
			fprintf(fp,"%d ",pixels[i + j * WIDTH]);
		}
		fprintf(fp,"\n");
	}
	
	
	//fprintf(fp, "Testing...\n");
	
	fclose(fp);
}

int computePoint(int _x, int _y) {
	
	int iteration, color;
	double xtemp;
	double x0, y0, x, y;
		
	x0 = (((double)_x - 1024) / ((double)1024 / (double)2));
	y0 = (((double)_y - 512) / ((double)HEIGHT / (double)2));
	
	/*
	if ((_x == 0) && (_y == 1023)) {
		fprintf(stderr,"X: %d, Y: %d, ScaledX: %lf, ScaledY: %lf\n",_x,_y,x0,y0);
	}
	if ((_x == 1023) && (_y == 1023)) {
		fprintf(stderr,"X: %d, Y: %d, ScaledX: %lf, ScaledY: %lf\n",_x,_y,x0,y0);
	}
	if ((_x == 1024) && (_y == 1023)) {
		fprintf(stderr,"X: %d, Y: %d, ScaledX: %lf, ScaledY: %lf\n",_x,_y,x0,y0);
	}
	if ((_x == 1535) && (_y == 1023)) {
		fprintf(stderr,"X: %d, Y: %d, ScaledX: %lf, ScaledY: %lf\n",_x,_y,x0,y0);
	}
	*/
	
	iteration = 0;
	x = 0;
	y = 0;
	
	while (((x*x + y*y) < 4) && (iteration < MAX_ITERATIONS)) //Remember: 4 == (2*2)
	{
		xtemp = x*x - y*y + x0;
		y = 2*x*y + y0;
		
		x = xtemp;
		
		iteration++;
	}
	
	color = MAX_ITERATIONS - iteration;
	return color;
	
}



int main(int argc, char** argv) {

	// mandelbrot
	int* pixels;
	int i, j;
	int color;
	
	// mpi
	int nproc, iproc;
	int* mySendArr;
	int* myRecvArr;
	#ifdef USE_MPI
	MPI_Status status;
	#endif
	
	// miscellaneous
	int loop; // used as boolean
	int loopCount;

	
	#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	#endif 
	
	// everyone needs a send and a recv array
	// the extra 1 is for the start location
	mySendArr = (int*)malloc((BLOCK_WIDTH * BLOCK_HEIGHT + 1) * sizeof(int));
	myRecvArr = (int*)malloc((BLOCK_WIDTH * BLOCK_HEIGHT + 1) * sizeof(int));
	int iter;
	
	if (iproc == 0) {	// master code
		int numJobs;
		int jobCount;
		int jobStart;
		double timestart, timefinish, timetaken;
		
		numJobs = NUM_JOBS;
		jobCount = 0;
		
		//fprintf(stderr,"(%d) I'm the master\n",iproc);
		pixels = (int*)malloc(WIDTH * HEIGHT * sizeof(int));
		
		timestart = When();
		
	/*
		for (j = 0; j < HEIGHT; j++) {
			for (i = 0; i < WIDTH; i++) {
				pixels[i + j * WIDTH] = computePoint(i, j);
			}
		}
	*/	
		//loop = 1;
		for (loopCount = 0; loopCount < (NUM_JOBS * 2 + nproc - 1); loopCount++) {
			//fprintf(stderr,"(%d) I'm waiting\n",iproc);
		  #ifdef USE_MPI
			
			MPI_Recv(myRecvArr,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
		//	fprintf(stderr,"(%d) %d sent me a message: %lf\n",iproc,status.MPI_SOURCE,myRecvArr[0]);
			
			if (myRecvArr[0] == -1) {	// worker wants more work
				//fprintf(stderr,"(%d) %d wants more work\n",iproc,status.MPI_SOURCE);
				
				if (numJobs > 0) {
					// tell worker there is a job starting at numJobs - 1 * BLOCK_WIDTH
					mySendArr[0] = jobCount * BLOCK_WIDTH;
					jobCount++;
					MPI_Send(mySendArr,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
					numJobs--;
				}
				else {
					// tell worker there isn't any more
					mySendArr[0] = -1;
					MPI_Send(mySendArr,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
				}
				
				
			}
			else if (myRecvArr[0] == -2) {	// worker wants to send finished work
				//fprintf(stderr,"(%d) %d wants to send me their work\n",iproc,status.MPI_SOURCE);
			
				MPI_Recv(myRecvArr,BLOCK_WIDTH * BLOCK_HEIGHT + 1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD,&status);
				
				// worker sent work with start number
				//fprintf(stderr,"(%d) %d sent me their work starting at %d\n",iproc,status.MPI_SOURCE,myRecvArr[0]);
				
				// copy work into pixel array
				jobStart = myRecvArr[0];
				for (i = 1; i < BLOCK_WIDTH * BLOCK_HEIGHT + 1; i++) {
					pixels[i-1 + jobStart] = myRecvArr[i];
					//fprintf(stderr,"%d ",pixels[i-1 + jobStart]);
				}
			}
#endif
			
			//loop = 0;
		}
		
		
		timefinish = When();
		timetaken = timefinish - timestart;
		//fprintf(stderr,"(%d) I'm finished\n",iproc);
		fprintf(stdout,"(%d) Time taken: %lf\n",iproc,timetaken);
	
		fileWrite(pixels);
		free(pixels);
	
	
	}
	else {	// worker code
		int myJobStart;
		//fprintf(stderr,"\t(%d) I'm a worker\n",iproc);
		pixels = (int*)malloc(BLOCK_WIDTH * BLOCK_HEIGHT * sizeof(int));
		
		loop = 1;
		while ((loop) && (iter < 1000) ) {
			// ask master for work
			mySendArr[0] = -1;
			#ifdef USE_MPI
			MPI_Send(mySendArr,1,MPI_INT,0,0,MPI_COMM_WORLD);
			#endif
			
			// recv response (starting number or -1)
			#ifdef USE_MPI
			MPI_Recv(myRecvArr,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
			#endif
			
			if (myRecvArr[0] == -1) {	// -1 means no more
				//fprintf(stderr,"\t(%d) Master says no more work\n",iproc);

                        #ifdef USE_MPI
			   loop = 0;
                         #else
			iter++;
			#endif 
			 
			}
			
			else {
				//fprintf(stderr,"\t(%d) Master gave me start of %d\n",iproc,myRecvArr[0]);
				myJobStart = myRecvArr[0];
				
				
				// do work
				#pragma omp parallel for private(i, j)
				for (j = 0; j < BLOCK_HEIGHT; j++) {
					for (i = 0; i < BLOCK_WIDTH; i++) {
						pixels[i + j * BLOCK_WIDTH] = computePoint(i, j + myJobStart / 1536);
						//fprintf(stderr,"%d ",pixels[i + j * BLOCK_WIDTH]);
					}
				}				

				
				// tell master work is done and ready to send
				mySendArr[0] = -2;
				#ifdef USE_MPI
				MPI_Send(mySendArr,1,MPI_INT,0,0,MPI_COMM_WORLD);
				#endif 
				// send work
				mySendArr[0] = myJobStart;
				for (i = 1; i < BLOCK_WIDTH * BLOCK_HEIGHT + 1; i++) {
					mySendArr[i] = pixels[i-1];
				}
				#ifdef USE_MPI
				MPI_Send(mySendArr,BLOCK_WIDTH * BLOCK_HEIGHT + 1,MPI_INT,0,0,MPI_COMM_WORLD);
				#endif 
				
			}	
		}
			
		//fprintf(stderr,"\t(%d) I'm finished\n",iproc);
		
	}
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
	
	return 0;
}

