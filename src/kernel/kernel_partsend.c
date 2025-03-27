/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include <omp.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdlib>
#include <random>
#include <ctime>
#include <getopt.h>
#include <string>

#include "kernel_partsend.h"



double doTest21( int rank, int, char*, char*, int, size_t, double, double );
int kernel_partsend_doTest21_func( int thread_num, int num_threads, void* args);

struct Args {
	int threadModel;
	int test;
	int numIterations;
	int numThreads;
	int minNumThreads;
	int maxNumThreads;
	size_t maxBuffSize;
	size_t minBuffSize;
	double compTime;
	double noise;
}; 

struct kernel_partsend_doTest21_args {
    int rank;
    int numIterations;
    long sleep;
    long sleepPlus;
    gran_request* gran_req;
#ifdef MPI_PARALLEL
    char* recvBuf;
    size_t threadPart;
    MPIX_Request myReq;
#endif
};

typedef struct
{
	Args* args;
} kernel_partsend_data;


int kernel_partsend_argc_temp;
char** kernel_partsend_argv_temp;
char** gran_send_bufs;
char** gran_recv_bufs;

std::string getThreadModelStr(int model)
{
	switch(model) {
		case MPI_THREAD_MULTIPLE: return "Multiple";
		case MPI_THREAD_SERIALIZED: return "Serial";
	}
	assert(0);
}

int parseArgs( int argc, char* argv[], Args* args )
{
        #define OPTIONAL 2
        static struct option long_options[] = {
                {"i", OPTIONAL, NULL, 0 },
                {"minT", OPTIONAL, NULL, 1 },
                {"maxT", OPTIONAL, NULL, 2 },
                {"minB", OPTIONAL, NULL, 3 },
                {"maxB", OPTIONAL, NULL, 4 },
                {"test", OPTIONAL, NULL, 5 },
                {"numT", OPTIONAL, NULL, 6 },
                {"tm", OPTIONAL, NULL, 7 },
                {"ct", OPTIONAL, NULL, 8 },
                {"n", OPTIONAL, NULL, 9 },
                {0,0,0,0}
        };

        while (1) {
                int option_index = 0;
                int c =  getopt_long( argc, argv, "", long_options, &option_index );
                if ( c == -1) break;
                switch( c ) {
                        case 0:
                                args->numIterations = atoi(optarg);
                                break;
                        case 1:
                                args->minNumThreads = atoi(optarg);
                                break;
                        case 2:
                                args->maxNumThreads = atoi(optarg);
                                break;
                        case 3:
                                args->minBuffSize = atoi(optarg);
                                break;
                        case 4:
                                args->maxBuffSize = atoi(optarg);
                                break;
                        case 5:
                                args->test = atoi(optarg);
                                break;
                        case 6:
                                args->numThreads = atoi(optarg);
                                break;
                        case 7:
                                args->threadModel = MPI_THREAD_MULTIPLE;
                                break;
                        case 8:
                                args->compTime = atof(optarg);
                                break;
                        case 9:
                                args->noise = atof(optarg);
                                break;
                        default:
                                return -1;
                }
        }
        return 0;
}

int kernel_partsend_parse_args(void* _data, int argc, char** argv)
{
	kernel_partsend_argc_temp = argc;

	kernel_partsend_argv_temp = (char**)malloc((argc + 1) * sizeof(char*));
	for(int i = 0; i < argc; i++)
	{
		size_t length = strlen(argv[i])+1;
		kernel_partsend_argv_temp[i] = (char*)malloc(length);
		memcpy(kernel_partsend_argv_temp[i], argv[i], length);
	}
	kernel_partsend_argv_temp[argc] = NULL;
	return(0);
}

int kernel_partsend_init(void** kernel_data, int id, int job_size)
{
        int rc;

	kernel_partsend_data* data = (kernel_partsend_data*) malloc(sizeof(kernel_partsend_data));
	Args* args = (Args*) malloc(sizeof(Args));
	args->threadModel = MPI_THREAD_MULTIPLE;
	args->test = 21;
	args->numThreads = 0;
	args->numIterations = 20;
	args->minNumThreads = 1;
	args->maxNumThreads = 2;

	args->minBuffSize = 1;

	args->maxBuffSize = 1048576;
	args->compTime = 0.050;

	args->noise = 0.04;

	
	
	
	

	if ( args->numThreads != 0 ) {
		args->minNumThreads = args->maxNumThreads = args->numThreads;
	}

	data->args = args;
	*kernel_data = (void*)data;
	return 0;
}

int kernel_partsend_run(void* _kernel_data)
{
    int provided;
    int rank, size;
    int rc;

    kernel_partsend_data* data = (kernel_partsend_data*) _kernel_data;
    Args* args = data->args;




    




    comm_get_id(&rank);

	if ( 0 == rank ) {
        printf("# OMP_NUM_THREADS=%s\n",getenv("OMP_NUM_THREADS")?getenv("OMP_NUM_THREADS"):" not set");

		printf("# test=%d iterations=%d minThreads=%d maxThreads=%d minBuf=%lu maxBuf=%lu\n",args->test,args->numIterations,
							args->minNumThreads,args->maxNumThreads,args->minBuffSize,args->maxBuffSize);
		printf("compTime=%f noise=%f\n",args->compTime,args->noise);
	}
	char hostname[100];
	gethostname(hostname,100);
	printf("# rank=%d hostname=%s\n",rank,hostname);



    comm_get_job_size(&size);

	assert( size == 2 ); 
 
	char* sendBuf = (char*) malloc( args->maxBuffSize );
	char* recvBuf = (char*) malloc( args->maxBuffSize );

	for ( size_t bufSize = args->minBuffSize; bufSize <= args->maxBuffSize; bufSize *= 2 ) {
		for ( int numThreads = args->minNumThreads; numThreads <= args->maxNumThreads; numThreads*=2 ) { 
			size_t threadPart = bufSize/numThreads;
			

			double duration; 
			double start = MPI_Wtime();

			comm_barrier();
                        duration = doTest21( rank, args->numIterations, sendBuf, recvBuf, numThreads, threadPart, args->compTime, args->noise );
			double stop = MPI_Wtime();
			if ( 0 == rank ) {
				double bandwidth = bufSize * args->numIterations / duration; 
				bandwidth /= 1000000.0;

				double msgRate = args->numIterations / duration;

				printf("bufSize %lu, numThreads %d, threadPart %lu, time %f, bandwidth %.3f MB/s, msgRate %.3f msg/sec, %lf\n",
					bufSize,numThreads,threadPart, duration, bandwidth, msgRate, (stop-start)/args->numIterations );
			}
		}
		int numThreads = 56;
		size_t threadPart = bufSize/numThreads;
		

		double duration; 
		double start = MPI_Wtime();

		comm_barrier();
                duration = doTest21( rank, args->numIterations, sendBuf, recvBuf, numThreads, threadPart, args->compTime, args->noise );
		double stop = MPI_Wtime();
		if ( 0 == rank ) {
			double bandwidth = bufSize * args->numIterations / duration; 
			bandwidth /= 1000000.0;

			double msgRate = args->numIterations / duration;

			printf("bufSize %lu, numThreads %d, threadPart %lu, time %f, bandwidth %.3f MB/s, msgRate %.3f msg/sec, %lf\n",
				bufSize,numThreads,threadPart, duration, bandwidth, msgRate, (stop-start)/args->numIterations );
		}
	}




    printf("# %d: rank=%d exiting\n",__LINE__,rank);

    return 0;
}

int kernel_partsend_clean(void* _kernel_data)
{
  return 0; 
}

int kernel_partsend_load(kernelfuncs* funcs)
{
  funcs->parse_args = &kernel_partsend_parse_args;
  funcs->init       = &kernel_partsend_init;
  funcs->run        = &kernel_partsend_run;
  funcs->clean      = &kernel_partsend_clean;
  return 0;
}

double doTest21( int rank, int numIterations, char* sendBuf, char* recvBuf, int numThreads, size_t threadPart, double compTime, double noise )
{
#ifdef MPI_PARALLEL
    MPIX_Request myReq;  
    int TAG = 0xdead;
#endif
    int other = (rank + 1) % 2;
    double start;
    int rc;
    int err;

#ifdef VERIFY
    for ( int i = 0; i < (threadPart * numThreads) / 8; i++ ) 
    {
        ((uint64_t*)sendBuf)[i] = i;
    }

    bzero( recvBuf, threadPart * numThreads );
#endif

    
    gran_request* gran_req;

    
    gran_send_bufs = (char**) malloc(sizeof(char*));
    gran_recv_bufs = (char**) malloc(sizeof(char*));

    
    gran_send_bufs[0] = (char*) malloc(threadPart*numThreads);

    
    void*** recv_bufs = (void***)malloc(sizeof(void**) * 1);
    void** send_bufs = (void**)malloc(sizeof(void*) * 1);
    int* recv_ids = (int*)malloc(sizeof(int) * 1);
    int* send_ids = (int*)malloc(sizeof(int) * 1);
    int* num_entries = (int*)malloc(sizeof(int) * 1);
    int* entry_sizes = (int*)malloc(sizeof(int) * 1);

    
    if ( rank == 0 ) 
    {
#ifdef MPI_PARALLEL
        rc = MPIX_Psend_init(sendBuf,
                             numThreads,
                             threadPart,
                             MPI_CHAR,
                             other,
                             TAG,
                             MPI_COMM_WORLD,
                             &myReq);

        assert( rc == MPI_SUCCESS );
#endif
        recv_bufs[0] = (void**)&gran_recv_bufs[0];
        send_bufs[0] = (void*)gran_send_bufs[0];
        recv_ids[0] = -1;
        send_ids[0] = other;
        num_entries[0] = numThreads;
        entry_sizes[0] = threadPart;
    } 
    else 
    {
#ifdef MPI_PARALLEL
        rc = MPIX_Precv_init(recvBuf,
                             numThreads, 
                             threadPart, 
                             MPI_CHAR,
                             other,
                             TAG,
                             MPI_COMM_WORLD,
                             &myReq);
        assert( rc == MPI_SUCCESS );
#endif
        recv_bufs[0] = (void**)&gran_recv_bufs[0];
        send_bufs[0] = (void*)gran_send_bufs[0];
        recv_ids[0] = other;
        send_ids[0] = -1;
        num_entries[0] = numThreads;
        entry_sizes[0] = threadPart;
    }

    
    err = gran_stencil_init(&gran_req, recv_bufs, send_bufs, recv_ids, send_ids, num_entries, entry_sizes, 1);
    if (err) { 
        printf("error in stencil_init:%d\n",err);
        assert(err == 0);
    }

    
    free(recv_bufs);
    free(send_bufs);
    free(recv_ids);
    free(send_ids);
    free(num_entries);
    free(entry_sizes);

    
    threads_request* t_req;
    threads_init(&t_req, numThreads);

    start = MPI_Wtime();
    srand(time(NULL));

    long sleep = compTime * 1000000000;
    long sleepPlus = (compTime + ( compTime * noise)) * 1000000000 ;

    
    kernel_partsend_doTest21_args* args = (kernel_partsend_doTest21_args*) malloc(sizeof(kernel_partsend_doTest21_args));
    args->rank = rank;
    args->numIterations = numIterations;
    args->sleep = sleep;
    args->sleepPlus = sleepPlus;
    args->gran_req = gran_req;
#ifdef MPI_PARALLEL
    args->recvBuf = recvBuf;
    args->threadPart = threadPart;
    args->myReq = myReq;
#endif
    threads_run(t_req, &kernel_partsend_doTest21_func, args);

#ifdef SINGLE_THREADED

    {
        int rc;
        int tid = omp_get_thread_num();
        int iteration;
        struct timespec req,rem;
        req.tv_sec = 0;
        if ( numThreads > 1 && tid == numThreads - 1 ) 
        {
            req.tv_nsec = sleepPlus;
        } 
        else 
        {
            req.tv_nsec = sleep;
    	}
	
        for ( iteration = 0; iteration < numIterations; iteration++ ) 
        {

#ifdef MPI_PARALLEL

            {
                rc = MPIX_Start(&myReq);
                
                assert( rc == MPI_SUCCESS );
                MPI_Barrier( MPI_COMM_WORLD );
            }
#endif
            gran_stencil_start(gran_req);

            threads_synch(t_req);


            if ( 0 == rank ) 
            {
                rc = clock_nanosleep(CLOCK_REALTIME,0,&req, &rem);
                if ( 0  !=  rc ) 
                {
                    printf("rc=%s rem %li\n",strerror(rc),rem.tv_nsec);
                }
#ifdef MPI_PARALLEL
                rc = MPIX_Pready(tid, &myReq); 
                assert( rc  == MPI_SUCCESS );
#endif
                err = gran_stencil_ready(gran_req, 0, tid);
                assert( err == 0 );
            }


          threads_synch(t_req);


#ifdef MPI_PARALLEL

            {
                rc = MPIX_Wait(&myReq, NULL);
                assert( rc == MPI_SUCCESS );
            }
#endif
            err = gran_stencil_end(gran_req);
            assert( err == 0 );

#ifdef VERIFY


            threads_synch(t_req);


            if( 1 == rank )
            {
                for ( int i = 0; i < (threadPart * numThreads) / 8; i++ )
                {
                    assert( ((uint64_t*)recvBuf)[i] == i );
                }
                bzero( recvBuf, threadPart * numThreads );
           }
#endif

        }
    }
#endif

    
    double duration = MPI_Wtime() - start;

    if ( numThreads > 1 ) 
    {
        duration -=  sleepPlus / 1000000000.0 * numIterations;
    } 
    else 
    {
        duration -=  sleep / 1000000000.0 * numIterations;
    }
    return duration;  
}

int kernel_partsend_doTest21_func( int thread_num, int num_threads, void* args)
{
    kernel_partsend_doTest21_args* doTest21_args = (kernel_partsend_doTest21_args*) args;
    int rank = doTest21_args->rank;
    int numIterations = doTest21_args->numIterations;
    long sleep = doTest21_args->sleep;
    long sleepPlus = doTest21_args->sleepPlus;
    gran_request* gran_req = doTest21_args->gran_req;
#ifdef MPI_PARALLEL
    char* recvBuf = doTest21_args->recvBuf;
    size_t threadPart = doTest21_args->threadPart;
    MPIX_Request myReq = doTest21_args->myReq;
#endif

    int rc;
    int err;
    int tid = thread_num;
    int numThreads = num_threads;
    int iteration;
    struct timespec req,rem;
    req.tv_sec = 0;
    if ( numThreads > 1 && tid == numThreads - 1 ) 
    {
        req.tv_nsec = sleepPlus;
    } 
    else 
    {
        req.tv_nsec = sleep;
    }
    
    for ( iteration = 0; iteration < numIterations; iteration++ ) 
    {

#ifdef MPI_PARALLEL

        {
            rc = MPIX_Start(&myReq);
            
            assert( rc == MPI_SUCCESS );
            MPI_Barrier( MPI_COMM_WORLD );
        }
#endif
        if (tid == 0)
        {
            err = gran_stencil_start(gran_req);
            assert( err == 0 );
        }


        threads_synch(NULL);

        if ( 0 == rank ) 
        {
            rc = clock_nanosleep(CLOCK_REALTIME,0,&req, &rem);
            if ( 0  !=  rc ) 
            {
                printf("rc=%s rem %li\n",strerror(rc),rem.tv_nsec);
            }
#ifdef MPI_PARALLEL
            rc = MPIX_Pready(tid, &myReq); 
            assert( rc  == MPI_SUCCESS );
#endif
            err = gran_stencil_ready(gran_req, 0, tid);
            assert( err == 0 );
        }


        threads_synch(NULL);

#ifdef MPI_PARALLEL

        {
            rc = MPIX_Wait(&myReq, NULL);
            assert( rc == MPI_SUCCESS );
        }
#endif
        if (tid == 0)
        {
            err = gran_stencil_end(gran_req);
            assert( err == 0 );
        }

#ifdef VERIFY

        threads_synch(NULL);

        if( 1 == rank )
        {
            for ( int i = 0; i < (threadPart * numThreads) / 8; i++ )
            {
                assert( ((uint64_t*)recvBuf)[i] == i );
            }
            bzero( recvBuf, threadPart * numThreads );
        }
#endif

    }
    return 0;
}
