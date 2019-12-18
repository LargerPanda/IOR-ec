#if !defined(_EC_H)
#define _EC_H

#include "iordef.h"

/*jerasure include*/
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <gf_rand.h>
#include <unistd.h>
#include <jerasure.h>
#include <reed_sol.h>
#include <cauchy.h>
#include <liberation.h>
#include <pthread.h>
//#include <timing.h>

#define N 10
//define k*m
#define K 2
#define M 2
#define W 8
#define P 8
#define TOTAL_STRIPE_NUM (K+M)
#define IMIDIATE_EC 1
//word count
// #define W 8  
// //0 means chosen adaptively 
// #define EC_BUFFER_SIZE 0       
// //size of k+m files
// #define EC_BLOCK_SIZE (1048576/K)
// //ignore in reed_sol's
// #define EC_PACKET_SIZE 8


#define info_flag 0

typedef enum Coding_Technique
{
    Reed_Sol_Van,
    Reed_Sol_R6_Op,
    Cauchy_Orig,
    Cauchy_Good,
    Liberation,
    Blaum_Roth,
    Liber8tion,
    RDP,
    EVENODD,
    No_Coding
}Coding_Technique;

typedef enum EC_Strategy{
    Immediate
}EC_Strategy;

//char *Methods[N] = {"reed_sol_van", "reed_sol_r6_op", "cauchy_orig", "cauchy_good", "liberation", "blaum_roth", "liber8tion", "no_coding"};

// enum Coding_Technique ec_method = Cauchy_Good; //here to point the coding method

/*a global structure to configure target servers*/
typedef struct FileList
{
    /* data */
    char *file[TOTAL_STRIPE_NUM];
}FileList;

typedef struct FdList
{
    /* data */
    void *fd[TOTAL_STRIPE_NUM];
}FdList;


/*need to allocate space in every process*/
typedef struct ec_info
{
    /* data */
    double *writeVal[2];/* array to write results */
    double *readVal[2];
}ec_info;

typedef struct ec_read_timer
{
    /* data */
    double readTotalTime;
}ec_read_timer;

typedef struct ec_read_thread_args
{
    /* data */
    int id;
    IOR_param_t *test;
    FdList *fds;
    int access;
    ec_read_timer *ec_timer;
    IOR_offset_t transfer;
    
    /* ec decode parameters */
    char **ec_data;
    char **ec_coding;
    int *ec_matrix;
    int *ec_bitmatrix;
    enum Coding_Technique method;

    /*ec control parameters*/
} ec_read_thread_args;

#endif // _EC_H
