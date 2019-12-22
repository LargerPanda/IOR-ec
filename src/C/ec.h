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
#define IMMEDIATE_EC 0

typedef struct ec_read_thread_args
{
    /* data */
    int id;
    IOR_param_t *test;
    void **fds;
    int access;
    IOR_offset_t *offSetArray;
    /* ec decode parameters */
    char **ec_data;
    char **ec_coding;
    int *ec_matrix;
    int *ec_bitmatrix;
    enum Coding_Technique method;

    /*ec control parameters*/
} ec_read_thread_args;


#endif // _EC_H
