#if !defined(_EC_H)
#define _EC_H

#include "iordef.h"

//define k*m
#define K 2
#define M 2
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
