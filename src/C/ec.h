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
#include <timing.h>

#define N 10
//define k*m
#define K 2
#define M 2
//word count
// #define W 8  
// //0 means chosen adaptively 
// #define EC_BUFFER_SIZE 0       
// //size of k+m files
// #define EC_BLOCK_SIZE (1048576/K)
// //ignore in reed_sol's
// #define EC_PACKET_SIZE 8      


enum Coding_Technique
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
};

char *Methods[N] = {"reed_sol_van", "reed_sol_r6_op", "cauchy_orig", "cauchy_good", "liberation", "blaum_roth", "liber8tion", "no_coding"};

// enum Coding_Technique ec_method = Cauchy_Good; //here to point the coding method

/*a global structure to configure target servers*/
typedef struct FileList
{
    /* data */
    char file1[MAX_STR] = "/data/data1/ec_testfile.part1";
    char file2[MAX_STR] = "/data/data2/ec_testfile.part2";
    char file3[MAX_STR] = "/data/data3/ec_testfile.parity1";
    char file4[MAX_STR] = "/data/data4/ec_testfile.parity2";
}FileList;

typedef struct FdList
{
    /* data */
    void *fd1 = NULL;
    void *fd2 = NULL;
    void *fd3 = NULL;
    void *fd4 = NULL;
}FdList;


/*need to allocate space in every process*/
typedef struct ec_info
{
    /* data */
    double *writeVal[2]; /* array to write results */
    double *readVal[2];
}*ec_info;



#endif // _EC_H
