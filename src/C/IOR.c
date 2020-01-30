/******************************************************************************\
*                                                                              *
*        Copyright (c) 2003, The Regents of the University of California       *
*      See the file COPYRIGHT for a complete copyright notice and license.     *
*                                                                              *
********************************************************************************
*
* CVS info:
*   $RCSfile: IOR.c,v $
*   $Revision: 1.4 $
*   $Date: 2008/12/03 16:16:02 $
*   $Author: rklundt $
*
* Purpose:
*       Test file system I/O.
*
* Settings and Usage:
*       View DisplayUsage() for settings.
*       Usage is with either by commandline or with an input script
*       file of the form shown in DisplayUsage().
*
* History (see CVS log for detailed history):
*       2001.11.21  wel
*                   Started initial implementation of code.
*       2002.02.07  wel
*                   Added abstract IOR interface for I/O.
*       2002.03.29  wel
*                   Added MPI synchronization.
*       2002.04.15  wel
*                   Added MPIIO.
*       2002.08.07  wel
*                   Added MPI file hints, collective calls, file views, etc.
*       2002.11.06  wel
*                   Added HDF5.
*       2003.10.03  wel
*                   Added NCMPI.
*
* Known problems and limitations:
*       None known.
*
\******************************************************************************/
/********************* Modifications to IOR-2.10.1 ****************************
* Hodson, 8/18/2008:                                                          *
* Documentation updated for the following new option:                         *
* The modifications to IOR-2.10.1 extend existing random I/O capabilities and *
* enhance performance output statistics, such as IOPs.                        *
*                                                                             *
*  cli        script                  Description                             *
* -----      -----------------        ----------------------------------------*
* 1)  -A N   testnum                - test reference number for easier test   *
*                                     identification in log files             *
* 2)  -Q N   taskpernodeoffset      - for read tests. Use with -C & -Z options*
*                                     If -C (reordertasks) specified,         *
*                                     then node offset read by CONSTANT  N.   *
*                                     If -Z (reordertasksrandom) specified,   *
*                                     then node offset read by RANDOM >= N.   *
* 3)  -Z     reordertasksrandom     - random node task ordering for read tests*
*                                     In this case, processes read            *
*                                     data written by other processes with    *
*                                     node offsets specified by the -Q option *
*                                     and -X option.                          *
* 4)  -X N   reordertasksrandomseed - random seed for -Z (reordertasksrandom) *
*                                     If N>=0, use same seed for all iters    *
*                                     If N< 0, use different seed for ea. iter*
* 5)  -Y     fsyncperwrite          - perform fsync after every POSIX write   *
\*****************************************************************************/

#include "aiori.h"                                    /* IOR I/O interfaces */
#include "IOR.h"                                      /* IOR definitions
                                                         and prototypes */
#include "IOR-aiori.h"                                /* IOR abstract
                                                         interfaces */
#include <ctype.h>                                    /* tolower() */
#include <errno.h>                                    /* sys_errlist */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>                                 /* struct stat */
#include <time.h>
#include "ec.h"                                      /* ec header */
#include <omp.h>
#include <sys/types.h>
#include <assert.h>
#include "mylist.h"
#ifndef _WIN32
#   include <sys/time.h>                              /* gettimeofday() */
#   include <sys/utsname.h>                           /* uname() */
#endif

/************************** D E C L A R A T I O N S ***************************/

extern IOR_param_t defaultParameters,
                   initialTestParams;
extern int     errno;                                   /* error number */
extern char ** environ;
int            totalErrorCount       = 0;
int            numTasksWorld         = 0;
int            rank                  = 0;
int            rankOffset            = 0;
int            tasksPerNode          = 0;               /* tasks per node */
int            verbose               = VERBOSE_0;       /* verbose output */
double         wall_clock_delta      = 0;
double         wall_clock_deviation;
MPI_Comm       testComm;



/********************thread_pool*********************/
//任务队列的结构体定义
typedef struct worker
{
    void *(*process)(void *arg); //一个返回值为void*的参数为void*的函数指针
    void *arg;                   //函数参数
    struct worker *next;         //指向下一个结构体
} CThread_worker;

//线程池的结构体定义
typedef struct
{
    pthread_mutex_t queue_lock; //互斥量
    pthread_cond_t queue_ready; //条件变量
    CThread_worker *queue_head; //线程任务队列的链表头指针
    int shutdown;               //是否销毁这个线程池
    pthread_t *threadid;        //指向线程id的指针
    int max_thread_num;         //最大线程数量
    int cur_queue_size;         //任务队列中的任务数
} CThread_pool;

//线程池相关的全局变量
int pool_add_worker(void *(*process)(void *arg), void *arg); //增加一个任务到任务队列
void *thread_routine(void *arg);                             //用于初始化线程的函数指针
static CThread_pool *pool = NULL;                            //指向线程池的指针

//线程池初始化函数
void pool_init(int max_thread_num)
{
    pool = (CThread_pool *)malloc(sizeof(CThread_pool));                      //初始化线程池指针
    pthread_mutex_init(&(pool->queue_lock), NULL);                            //初始化互斥量
    pthread_cond_init(&(pool->queue_ready), NULL);                            //初始化条件变量
    pool->queue_head = NULL;                                                  //初始化任务队列头指针
    pool->max_thread_num = max_thread_num;                                    //初始化最大线程数
    pool->cur_queue_size = 0;                                                 //初始化任务队列当前数量
    pool->shutdown = 0;                                                       //初始化是否销毁线程池的标志
    pool->threadid = (pthread_t *)malloc(max_thread_num * sizeof(pthread_t)); //初始化线程id
    int i = 0;
    for (i = 0; i < max_thread_num; i++)
    {
        pthread_create(&(pool->threadid[i]), NULL, thread_routine, NULL); //创建线程
    }
}

//执行任务队列的任务的函数，用于初始化线程
/*
在调用pthread_cond_wait()前必须由本线程加锁（pthread_mutex_lock()），而在更新条件等待队列以前，mutex保持锁定状态，
并在线程挂起进入等待前解锁。在条件满足从而离开pthread_cond_wait()之前，mutex将被重新加锁，
以与进入pthread_cond_wait()前的加锁动作对应。
*/
void *thread_routine(void *arg)
{
    printf("start thread 0x%x\n", pthread_self());
    while (1)
    {
        pthread_mutex_lock(&(pool->queue_lock)); //互斥量加锁
        while (pool->cur_queue_size == 0 && !pool->shutdown)
        { //等待条件发生
            //printf("thread 0x%x is waiting\n", pthread_self());
            pthread_cond_wait(&(pool->queue_ready), &(pool->queue_lock));
        }
        if (pool->shutdown)
        { //如果shutdown为1则线程退出
            pthread_mutex_unlock(&(pool->queue_lock));
            //printf("thread 0x%x will exit\n", pthread_self());
            pthread_exit(NULL);
        }
        //printf("thread 0x%x is starting to work\n", pthread_self());
        assert(pool->cur_queue_size != 0);
        assert(pool->queue_head != NULL);
        pool->cur_queue_size--;
        CThread_worker *worker = pool->queue_head; //将任务加入任务队列
        pool->queue_head = worker->next;
        pthread_mutex_unlock(&(pool->queue_lock)); //互斥量解锁
        (*(worker->process))(worker->arg);         //执行这个任务
        free(worker);                              //释放任务指针
        worker = NULL;
    }
    pthread_exit(NULL);
}

//添加一个任务到任务队列
int pool_add_worker(void *(*process)(void *arg), void *arg)
{
    CThread_worker *newworker = (CThread_worker *)malloc(sizeof(CThread_worker)); //创建队列的一个节点
    newworker->process = process;                                                 //对节点进行初始化
    newworker->arg = arg;
    newworker->next = NULL;
    pthread_mutex_lock(&(pool->queue_lock));   //互斥量上锁
    CThread_worker *member = pool->queue_head; //将这个节点插入任务队列中
    if (member != NULL)
    {
        while (member->next != NULL)
            member = member->next;
        member->next = newworker;
    }
    else
        pool->queue_head = newworker;
    assert(pool->queue_head != NULL);
    pool->cur_queue_size++;                    //任务队列个数加一
    pthread_mutex_unlock(&(pool->queue_lock)); //互斥量解锁
    pthread_cond_signal(&(pool->queue_ready)); //给条件变量发信号
    return 0;
}

//销毁线程
int pool_destroy()
{
    if (pool->shutdown)
        return -1;
    pool->shutdown = 1;
    pthread_cond_broadcast(&(pool->queue_ready)); //唤醒所有线程
    int i;
    for (i = 0; i < pool->max_thread_num; i++)
        pthread_join(pool->threadid[i], NULL); //以参与的方式回收线程资源
    free(pool->threadid);                      //释放线程id的指针
    CThread_worker *head = NULL;
    while (pool->queue_head != NULL)
    { //释放队列中所有的节点
        head = pool->queue_head;
        pool->queue_head = pool->queue_head->next;
        free(head);
    }
    pthread_mutex_destroy(&(pool->queue_lock)); //销毁互斥量
    pthread_cond_destroy(&(pool->queue_ready)); //销毁条件变量
    free(pool);                                 //释放线程池指针
    pool = NULL;                                //线程池指针置空
    return 0;
}
/********************************** M A I N ***********************************/

int
main(int     argc,
     char ** argv)
{
    int           i;
    IOR_queue_t * tests;

    /*
     * check -h option from commandline without starting MPI;
     * if the help option is requested in a script file (showHelp=TRUE),
     * the help output will be displayed in the MPI job
     */
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            DisplayUsage(argv);
            return(0);
        }
    }

    /* start the MPI code */
    MPI_CHECK(MPI_Init(&argc, &argv), "cannot initialize MPI");
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &numTasksWorld),
              "cannot get number of tasks");
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "cannot get rank");
    /* set error-handling */
    /*MPI_CHECK(MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN),
      "cannot set errhandler");*/
    
    /* setup tests before verifying test validity */
    tests = SetupTests(argc, argv);
    verbose = tests->testParameters.verbose;
    tests->testParameters.testComm = MPI_COMM_WORLD;
    
    /* check for commandline usage */
    if (rank == 0 && tests->testParameters.showHelp == TRUE) {
      DisplayUsage(argv);
    }
    
    /* perform each test */

    while (tests != NULL) {
        verbose = tests->testParameters.verbose;
        if (rank == 0 && verbose >= VERBOSE_0) {
            ShowInfo(argc, argv, &tests->testParameters);
        }
        if (rank == 0 && verbose >= VERBOSE_3) {
            ShowTest(&tests->testParameters);
        }
        TestIoSys(&tests->testParameters);
        tests = tests->nextTest;
    }
    

    /* display finish time */
    if (rank == 0 && verbose >= VERBOSE_0) {
      fprintf(stdout, "Run finished: %s", CurrentTimeString());
    }

    MPI_CHECK(MPI_Finalize(), "cannot finalize MPI");

    return(totalErrorCount);

} /* main() */


/***************************** F U N C T I O N S ******************************/

/******************************************************************************/
/*
 * Bind abstract I/O function pointers to API-specific functions.
 */

void
AioriBind(char * api)
{
    if (strcmp(api, "POSIX") == 0) {
        IOR_Create      = IOR_Create_POSIX;
        IOR_Open        = IOR_Open_POSIX;
        IOR_Xfer        = IOR_Xfer_POSIX;
        IOR_Close       = IOR_Close_POSIX;
        IOR_Delete      = IOR_Delete_POSIX;
        IOR_SetVersion  = IOR_SetVersion_POSIX;
        IOR_Fsync       = IOR_Fsync_POSIX;
        IOR_GetFileSize = IOR_GetFileSize_POSIX;
    } else if (strcmp(api, "MPIIO") == 0) {
        IOR_Create      = IOR_Create_MPIIO;
        IOR_Open        = IOR_Open_MPIIO;
        IOR_Xfer        = IOR_Xfer_MPIIO;
        IOR_Close       = IOR_Close_MPIIO;
        IOR_Delete      = IOR_Delete_MPIIO;
        IOR_SetVersion  = IOR_SetVersion_MPIIO;
        IOR_Fsync       = IOR_Fsync_MPIIO;
        IOR_GetFileSize = IOR_GetFileSize_MPIIO;
    } else if (strcmp(api, "HDF5") == 0) {
        IOR_Create      = IOR_Create_HDF5;
        IOR_Open        = IOR_Open_HDF5;
        IOR_Xfer        = IOR_Xfer_HDF5;
        IOR_Close       = IOR_Close_HDF5;
        IOR_Delete      = IOR_Delete_HDF5;
        IOR_SetVersion  = IOR_SetVersion_HDF5;
        IOR_Fsync       = IOR_Fsync_HDF5;
        IOR_GetFileSize = IOR_GetFileSize_HDF5;
    } else if (strcmp(api, "NCMPI") == 0) {
        IOR_Create      = IOR_Create_NCMPI;
        IOR_Open        = IOR_Open_NCMPI;
        IOR_Xfer        = IOR_Xfer_NCMPI;
        IOR_Close       = IOR_Close_NCMPI;
        IOR_Delete      = IOR_Delete_NCMPI;
        IOR_SetVersion  = IOR_SetVersion_NCMPI;
        IOR_Fsync       = IOR_Fsync_NCMPI;
        IOR_GetFileSize = IOR_GetFileSize_NCMPI;
    } else {
        WARN("unrecognized IO API");
    }
} /* AioriBind() */


/******************************************************************************/
/*
 * Display outliers found.
 */

void
DisplayOutliers(int    numTasks,
                double timerVal,
                char * timeString,
                int    access,
                int    outlierThreshold)
{
    char accessString[MAX_STR];
    double sum, mean, sqrDiff, var, sd;

    /* for local timerVal, don't compensate for wall clock delta */
    timerVal += wall_clock_delta;

    MPI_CHECK(MPI_Allreduce(&timerVal, &sum, 1, MPI_DOUBLE, MPI_SUM, testComm),
              "MPI_Allreduce()");
    mean = sum / numTasks;
    sqrDiff = pow((mean - timerVal), 2);
    MPI_CHECK(MPI_Allreduce(&sqrDiff, &var, 1, MPI_DOUBLE, MPI_SUM, testComm),
              "MPI_Allreduce()");
    var = var / numTasks;
    sd = sqrt(var);

    if (access == WRITE) {
        strcpy(accessString, "write");
    } else { /* READ */
        strcpy(accessString, "read");
    }
    if (fabs(timerVal - mean) > (double)outlierThreshold){
        fprintf(stdout, "WARNING: for task %d, %s %s is %f\n",
                rank, accessString, timeString, timerVal);
        fprintf(stdout, "         (mean=%f, stddev=%f)\n", mean, sd);
        fflush(stdout);
    }
} /* DisplayOutliers() */


/******************************************************************************/
/*
 * Check for outliers in start/end times and elapsed create/xfer/close times.
 */

void
CheckForOutliers(IOR_param_t  * test,
                 double      ** timer,
                 int            rep,
                 int            access)
{
    int shift;

    if (access == WRITE) {
        shift = 0;
    } else { /* READ */
        shift = 6;
    }

    DisplayOutliers(test->numTasks, timer[shift+0][rep],
                    "start time", access, test->outlierThreshold);
    DisplayOutliers(test->numTasks, timer[shift+1][rep]-timer[shift+0][rep],
                    "elapsed create time", access, test->outlierThreshold);
    DisplayOutliers(test->numTasks, timer[shift+3][rep]-timer[shift+2][rep],
                    "elapsed transfer time", access, test->outlierThreshold);
    DisplayOutliers(test->numTasks, timer[shift+5][rep]-timer[shift+4][rep],
                    "elapsed close time", access, test->outlierThreshold);
    DisplayOutliers(test->numTasks, timer[shift+5][rep],
                    "end time", access, test->outlierThreshold);

} /* CheckForOutliers() */


/******************************************************************************/
/*
 * Check if actual file size equals expected size; if not use actual for
 * calculating performance rate.
 */

void
CheckFileSize(IOR_param_t *test,
              IOR_offset_t dataMoved,
              int          rep)
{
    MPI_CHECK(MPI_Allreduce(&dataMoved, &test->aggFileSizeFromXfer[rep],
                            1, MPI_LONG_LONG_INT, MPI_SUM, testComm),
              "cannot total data moved");

    if (strcmp(test->api, "HDF5") != 0 && strcmp(test->api, "NCMPI") != 0) {
        if (verbose >= VERBOSE_0 && rank == 0) {
            if ((test->aggFileSizeFromCalc[rep]
                 != test->aggFileSizeFromXfer[rep])
                ||
                (test->aggFileSizeFromStat[rep]
                 != test->aggFileSizeFromXfer[rep])) {
                fprintf(stdout,
                    "WARNING: Expected aggregate file size       = %lld.\n",
                    (long long)test->aggFileSizeFromCalc[rep]);
                fprintf(stdout,
                    "WARNING: Stat() of aggregate file size      = %lld.\n",
                    (long long)test->aggFileSizeFromStat[rep]);
                fprintf(stdout,
                    "WARNING: Using actual aggregate bytes moved = %lld.\n",
                    (long long)test->aggFileSizeFromXfer[rep]);
            }
        }
    }
    test->aggFileSizeForBW[rep] = test->aggFileSizeFromXfer[rep];

} /* CheckFileSize() */


/******************************************************************************/
/*
 * Check if string is true or false.
 */

char *
CheckTorF(char * string)
{
    string = LowerCase(string);
    if (strcmp(string, "false") == 0) {
        strcpy(string, "0");
    } else if (strcmp(string, "true") == 0) {
        strcpy(string, "1");
    }
    return(string);
} /* CheckTorF() */


/******************************************************************************/
/*
 * Compare buffers after reading/writing each transfer.  Displays only first
 * difference in buffers and returns total errors counted.
 */

size_t
CompareBuffers(void          * expectedBuffer,
               void          * unknownBuffer,
               size_t          size,
               IOR_offset_t    transferCount,
               IOR_param_t   * test,
               int             access)
{
    char                testFileName[MAXPATHLEN],
                        bufferLabel1[MAX_STR],
                        bufferLabel2[MAX_STR];
    size_t              i, j, length, first, last;
    size_t              errorCount  = 0;
    int                 inError = 0;
    unsigned long long *goodbuf = (unsigned long long *)expectedBuffer;
    unsigned long long *testbuf = (unsigned long long *)unknownBuffer;

    if (access == WRITECHECK) {
        strcpy(bufferLabel1, "Expected: ");
        strcpy(bufferLabel2, "Actual:   ");
    } else if (access == READCHECK) {
        strcpy(bufferLabel1, "1st Read: ");
        strcpy(bufferLabel2, "2nd Read: ");
    } else {
        ERR("incorrect argument for CompareBuffers()");
    }

    length = size / sizeof(IOR_size_t);
    first = -1;
    if (verbose >= VERBOSE_3) {
        fprintf(stdout,
                "[%d] At file byte offset %lld, comparing %llu-byte transfer\n",
                rank, test->offset, (long long)size);
    }
    for (i = 0; i < length; i++) {
        if (testbuf[i] != goodbuf[i]) {
            errorCount++;
            if (verbose >= VERBOSE_2) {
                fprintf(stdout,
          "[%d] At transfer buffer #%lld, index #%lld (file byte offset %lld):\n",
                        rank, transferCount-1, (long long)i,
                        test->offset + (IOR_size_t)(i * sizeof(IOR_size_t)));
                fprintf(stdout, "[%d] %s0x", rank, bufferLabel1);
                fprintf(stdout, "%016llx\n", goodbuf[i]);
                fprintf(stdout, "[%d] %s0x", rank, bufferLabel2);
                fprintf(stdout, "%016llx\n", testbuf[i]);
            }
            if (!inError) {
                inError = 1;
                first = i;
                last = i;
            } else {
                last = i;
            }
        } else if (verbose >= VERBOSE_5 && i % 4 == 0) {
            fprintf(stdout, "[%d] PASSED offset = %lld bytes, transfer %lld\n",
                    rank, ((i * sizeof(unsigned long long)) + test->offset),
                    transferCount);
            fprintf(stdout, "[%d] GOOD %s0x", rank, bufferLabel1);
            for (j = 0; j < 4; j++)
                fprintf(stdout, "%016llx ", goodbuf[i + j]);
            fprintf(stdout, "\n[%d] GOOD %s0x", rank, bufferLabel2);
            for (j = 0; j < 4; j++)
                fprintf(stdout, "%016llx ", testbuf[i + j]);
            fprintf(stdout, "\n");
        }
    }
    if (inError) {
        inError = 0;
        GetTestFileName(testFileName, test);
        fprintf(stdout,
            "[%d] FAILED comparison of buffer containing %d-byte ints:\n",
            rank, (int)sizeof(unsigned long long int));
        fprintf(stdout, "[%d]   File name = %s\n", rank, testFileName);
        fprintf(stdout, "[%d]   In transfer %lld, ",
                rank, transferCount);
        fprintf(stdout, "%lld errors between buffer indices %lld and %lld.\n",
                (long long)errorCount, (long long)first, (long long)last);
        fprintf(stdout, "[%d]   File byte offset = %lld:\n",
                rank, ((first * sizeof(unsigned long long)) + test->offset));

        fprintf(stdout, "[%d]     %s0x", rank, bufferLabel1);
        for (j = first; j < length && j < first + 4; j++)
            fprintf(stdout, "%016llx ", goodbuf[j]);
        if (j == length) fprintf(stdout, "[end of buffer]");
        fprintf(stdout, "\n[%d]     %s0x", rank, bufferLabel2);
        for (j = first; j < length && j < first + 4; j++)
            fprintf(stdout, "%016llx ", testbuf[j]);
        if (j == length) fprintf(stdout, "[end of buffer]");
        fprintf(stdout, "\n");
        if (test->quitOnError == TRUE)
            ERR("data check error, aborting execution");
    }
    return(errorCount);
} /* CompareBuffers() */


/******************************************************************************//*
 * Count all errors across all tasks; report errors found.
 */

int
CountErrors(IOR_param_t *test,
            int          access,
            int          errors)
{
    int allErrors = 0;

    if (test->checkWrite || test->checkRead) {
        MPI_CHECK(MPI_Reduce(&errors, &allErrors, 1, MPI_INT, MPI_SUM,
                             0, testComm), "cannot reduce errors");
        MPI_CHECK(MPI_Bcast(&allErrors, 1, MPI_INT, 0, testComm),
                  "cannot broadcast allErrors value");
        if (allErrors != 0) {
            totalErrorCount += allErrors;
            test->errorFound = TRUE;
        }
        if (rank == 0 && allErrors != 0) {
            if (allErrors < 0) {
                WARN("overflow in errors counted");
                allErrors = -1;
            }
            if (access == WRITECHECK) {
                fprintf(stdout, "WARNING: incorrect data on write.\n");
                fprintf(stdout,
                    "         %d errors found on write check.\n", allErrors);
            } else {
                fprintf(stdout, "WARNING: incorrect data on read.\n");
                fprintf(stdout,
                    "         %d errors found on read check.\n", allErrors);
            }
            fprintf(stdout, "Used Time Stamp %u (0x%x) for Data Signature\n",
                test->timeStampSignatureValue, test->timeStampSignatureValue);
        }
    }
    return(allErrors);
} /* CountErrors() */


/******************************************************************************/
/*
 * Compares hostnames to determine the number of tasks per node
 */

int
CountTasksPerNode(int numTasks, MPI_Comm comm)
{
    char       localhost[MAX_STR],
               hostname[MAX_STR],
               taskOnNode[MAX_STR];
    int        count               = 1,
               resultsLen          = MAX_STR,
               i;
    static int firstPass           = TRUE;
    MPI_Status status;

    MPI_CHECK(MPI_Get_processor_name(localhost, &resultsLen),
              "cannot get processor name");

    if (verbose >= VERBOSE_2 && firstPass) {
        sprintf(taskOnNode, "task %d on %s", rank, localhost);
        OutputToRoot(numTasks, comm, taskOnNode);
        firstPass = FALSE;
    }

    if (numTasks > 1) {
        if (rank == 0) {
            /* MPI_receive all hostnames, and compare to local hostname */
            for (i = 0; i < numTasks-1; i++) {
                MPI_CHECK(MPI_Recv(hostname, MAX_STR, MPI_CHAR, MPI_ANY_SOURCE,
                                   MPI_ANY_TAG, comm, &status),
                          "cannot receive hostnames");
                if (strcmp(hostname, localhost) == 0)
                    count++;
            }
        } else {
            /* MPI_send hostname to root node */
            MPI_CHECK(MPI_Send(localhost, MAX_STR, MPI_CHAR, 0, 0,
                               comm), "cannot send hostname");
        }
        MPI_CHECK(MPI_Bcast(&count, 1, MPI_INT, 0, comm),
                  "cannot broadcast tasks-per-node value");
    }

    return(count);
} /* CountTasksPerNode() */


/******************************************************************************/
/*
 * Allocate a page-aligned (required by O_DIRECT) buffer.
 */

void *
CreateBuffer(size_t size)
{
    size_t pageSize;
    size_t pageMask;
    char *buf, *tmp;
    char *aligned;

    pageSize = getpagesize();
    pageMask = pageSize-1;
    buf = malloc(size + pageSize + sizeof(void *));
    if (buf == NULL) ERR("out of memory");
    /* find the alinged buffer */
    tmp = buf + sizeof(char *);
    aligned = tmp + pageSize - ((size_t)tmp & pageMask);
    /* write a pointer to the original malloc()ed buffer into the bytes
       preceding "aligned", so that the aligned buffer can later be free()ed */
    tmp = aligned - sizeof(void *);
    *(void **)tmp = buf;
    
    return((void *)aligned);
} /* CreateBuffer() */


/******************************************************************************/
/*
 * Create new test for list of tests.
 */

IOR_queue_t *
CreateNewTest(int test_num)
{
    IOR_queue_t * newTest         = NULL;

    newTest = (IOR_queue_t *)malloc(sizeof(IOR_queue_t));
    if (newTest == NULL) ERR("out of memory");
    newTest->testParameters = initialTestParams;
    GetPlatformName(newTest->testParameters.platform);
    newTest->testParameters.nodes = initialTestParams.numTasks / tasksPerNode;
    newTest->testParameters.tasksPerNode = tasksPerNode;
    newTest->testParameters.id = test_num;
    newTest->nextTest = NULL;
    return(newTest);
} /* CreateNewTest() */


/******************************************************************************/
/*
 * Sleep for 'delay' seconds.
 */

void
DelaySecs(int delay)
{
    if (rank == 0 && delay > 0 ) {
        if (verbose >= VERBOSE_1)
            fprintf(stdout, "delaying %d seconds . . .\n", delay);
        sleep(delay);
    }
} /* DelaySecs() */


/******************************************************************************/
/*
 * Display freespace (df).
 */
void
DisplayFreespace(IOR_param_t  * test)
{
    char fileName[MAX_STR]   = {0};
    int  i;
    int  directoryFound      = FALSE;

    /* get outfile name */
    GetTestFileName(fileName, test);

    /* get directory for outfile */
    i = strlen(fileName);
    while (i-- > 0) {
        if (fileName[i] == '/') {
            fileName[i] = '\0';
            directoryFound = TRUE;
            break;
        }
    }

    /* if no directory/, use '.' */
    if (directoryFound == FALSE) {
        strcpy(fileName, ".");
    }

#if USE_UNDOC_OPT /* NFS */
    if (test->NFS_serverCount) {
        strcpy(fileName, test->NFS_rootPath);
    }
#endif /* USE_UNDOC_OPT - NFS */

    ShowFileSystemSize(fileName);

    return;
} /* DisplayFreespace() */


/******************************************************************************/
/*
 * Display usage of script file.
 */

void
DisplayUsage(char ** argv)
{
    char * opts[] = {
"OPTIONS:",
" -A N  testNum -- test number for reference in some output",
" -a S  api --  API for I/O [POSIX|MPIIO|HDF5|NCMPI]",
" -b N  blockSize -- contiguous bytes to write per task  (e.g.: 8, 4k, 2m, 1g)",
" -B    useO_DIRECT -- uses O_DIRECT for POSIX, bypassing I/O buffers",
" -c    collective -- collective I/O",
" -C    reorderTasks -- changes task ordering to n+1 ordering for readback",
" -Q N  taskPerNodeOffset for read tests use with -C & -Z options (-C constant N, -Z at least N)",
" -Z    reorderTasksRandom -- changes task ordering to random ordering for readback",
" -X N  reorderTasksRandomSeed -- random seed for -Z option",
" -d N  interTestDelay -- delay between reps in seconds",
" -D N  deadlineForStonewalling -- seconds before stopping write or read phase",
" -Y    fsyncPerWrite -- perform fsync after each POSIX write",
" -e    fsync -- perform fsync upon POSIX write close",
" -E    useExistingTestFile -- do not remove test file before write access",
" -f S  scriptFile -- test script name",
" -F    filePerProc -- file-per-process",
" -g    intraTestBarriers -- use barriers between open, write/read, and close",
" -G N  setTimeStampSignature -- set value for time stamp signature",
" -h    showHelp -- displays options and help",
" -H    showHints -- show hints",
" -i N  repetitions -- number of repetitions of test",
" -I    individualDataSets -- datasets not shared by all procs [not working]",
" -j N  outlierThreshold -- warn on outlier N seconds from mean",
" -J N  setAlignment -- HDF5 alignment in bytes (e.g.: 8, 4k, 2m, 1g)",
" -k    keepFile -- don't remove the test file(s) on program exit",
" -K    keepFileWithError  -- keep error-filled file(s) after data-checking",
" -l    storeFileOffset -- use file offset as stored signature",
" -m    multiFile -- use number of reps (-i) for multiple file count",
" -n    noFill -- no fill in HDF5 file creation",
" -N N  numTasks -- number of tasks that should participate in the test",
" -o S  testFile -- full name for test",
" -O S  string of IOR directives (e.g. -O checkRead=1,lustreStripeCount=32)",
" -p    preallocate -- preallocate file size",
" -P    useSharedFilePointer -- use shared file pointer [not working]",
" -q    quitOnError -- during file error-checking, abort on error",
" -r    readFile -- read existing file",
" -R    checkRead -- check read after read",
" -s N  segmentCount -- number of segments",
" -S    useStridedDatatype -- put strided access into datatype [not working]",
" -t N  transferSize -- size of transfer in bytes (e.g.: 8, 4k, 2m, 1g)",
" -T N  maxTimeDuration -- max time in minutes to run tests",
" -u    uniqueDir -- use unique directory name for each file-per-process",
" -U S  hintsFileName -- full name for hints file",
" -v    verbose -- output information (repeating flag increases level)",
" -V    useFileView -- use MPI_File_set_view",
" -w    writeFile -- write file",
" -W    checkWrite -- check read after write",
" -x    singleXferAttempt -- do not retry transfer if incomplete",
" -z    randomOffset -- access is to random, not sequential, offsets within a file",
" ",
"         NOTE: S is a string, N is an integer number.",
" ",
"" };
    int i = 0;

    fprintf(stdout, "Usage: %s [OPTIONS]\n\n", * argv);
    for (i=0; strlen(opts[i]) > 0; i++)
        fprintf(stdout, "%s\n", opts[i]);

    return;
} /* DisplayUsage() */


/******************************************************************************/
/*
 * Distribute IOR_HINTs to all tasks' environments.
 */

void
DistributeHints(void)
{
    char   hint[MAX_HINTS][MAX_STR],
           fullHint[MAX_STR],
           hintVariable[MAX_STR];
    int    hintCount                 = 0,
           i;

    if (rank == 0) {
        for (i = 0; environ[i] != NULL; i++) {
            if (strncmp(environ[i], "IOR_HINT", strlen("IOR_HINT")) == 0) {
                hintCount++;
                if (hintCount == MAX_HINTS) {
                    WARN("exceeded max hints; reset MAX_HINTS and recompile");
                    hintCount = MAX_HINTS;
                    break;
                }
                /* assume no IOR_HINT is greater than MAX_STR in length */
                strncpy(hint[hintCount-1], environ[i], MAX_STR-1);
            }
        }
    }

    MPI_CHECK(MPI_Bcast(&hintCount, sizeof(hintCount), MPI_BYTE,
                        0, MPI_COMM_WORLD), "cannot broadcast hints");
    for (i = 0; i < hintCount; i++) {
        MPI_CHECK(MPI_Bcast(&hint[i], MAX_STR, MPI_BYTE,
                            0, MPI_COMM_WORLD), "cannot broadcast hints");
        strcpy(fullHint, hint[i]);
        strcpy(hintVariable, strtok(fullHint, "="));
        if (getenv(hintVariable) == NULL) {
            /* doesn't exist in this task's environment; better set it */
            if (putenv(hint[i]) != 0) WARN("cannot set environment variable");
        }
    }
} /* DistributeHints() */


/******************************************************************************/
/*
 * Fill buffer, which is transfer size bytes long, with known 8-byte long long
 * int values.  In even-numbered 8-byte long long ints, store MPI task in high
 * bits and timestamp signature in low bits.  In odd-numbered 8-byte long long
 * ints, store transfer offset.  If storeFileOffset option is used, the file
 * (not transfer) offset is stored instead.
 */

void
FillBuffer(void               * buffer,
           IOR_param_t        * test,
           unsigned long long   offset,
           int                  fillrank)
{
    size_t              i;
    unsigned long long  hi, lo;
    unsigned long long *buf = (unsigned long long *)buffer;

    hi = ((unsigned long long)fillrank) << 32;
    lo = (unsigned long long)test->timeStampSignatureValue;
    for (i = 0; i < test->transferSize / sizeof(unsigned long long); i++) {
        if ((i%2) == 0) {
            /* evens contain MPI rank and time in seconds */
            buf[i] = hi | lo;
        } else {
            /* odds contain offset */
            buf[i] = offset + (i * sizeof(unsigned long long));
        }
    }
} /* FillBuffer() */

/******************************************************************************//*
 * Free transfer buffers.
 */

void
FreeBuffers(int           access,
            void         *checkBuffer,
            void         *readCheckBuffer,
            void         *buffer,
            IOR_offset_t *offsetArray)
{
    /* free() aligned buffers */
    if (access == WRITECHECK || access == READCHECK) {
        free( *(void **)((char *)checkBuffer - sizeof(char *)) );
    }
    if (access == READCHECK) {
        free( *(void **)((char *)readCheckBuffer - sizeof(char *)) );
    }
    free( *(void **)((char *)buffer - sizeof(char *)) );

    /* nothing special needed to free() this unaligned buffer */
    free(offsetArray);
    return;
} /* FreeBuffers() */


/******************************************************************************/
/*
 * Return string describing machine name and type.
 */
     
void
GetPlatformName(char * platformName)
{
    char             nodeName[MAX_STR],
                   * p,
                   * start,
                     sysName[MAX_STR];
    struct utsname   name;

    if (uname(&name) != 0) {
        WARN("cannot get platform name");
        sprintf(sysName, "%s", "Unknown");
        sprintf(nodeName, "%s", "Unknown");
    } else {
        sprintf(sysName, "%s", name.sysname);
        sprintf(nodeName, "%s", name.nodename);
    }

    start = nodeName;
    if (strlen(nodeName) == 0) {
        p = start;
    } else {
        /* point to one character back from '\0' */
        p = start + strlen(nodeName) - 1;
    }
    /*
     * to cut off trailing node number, search backwards
     * for the first non-numeric character
     */
    while (p != start) {
        if (* p < '0' || * p > '9') {
            * (p+1) = '\0';
            break;
        } else {
            p--;
        }
    }

    sprintf(platformName, "%s(%s)", nodeName, sysName);
} /* GetPlatformName() */


/******************************************************************************/
/*
 * Return test file name to access.
 * for single shared file, fileNames[0] is returned in testFileName
 */
/*this is a function to choose the start ost#*/
int ChooseStartOST(){
    return (rank*4);
}

void GetTestFileName_ec(char **ec_testFileNames, IOR_param_t *test)
{
    //format is /data/data(ost#)/filename.process#.strip# or parity#
    // if(test->ec_verbose >= VERBOSE_3){
    //     fprintf(stdout, "in GetTestFileName_ec\n");
    // }
    char initialTestFileName[MAX_STR],
         targetDirectoryRoot[MAX_STR] = "/data/ost";

    int total_stripe_num = test->ec_k+test->ec_m;
    int start_ost = ChooseStartOST();
    strcpy(initialTestFileName, test->testFileName);

    int i;
    if(test->filePerProc){
        for(i=0;i<total_stripe_num;i++){
            if(i<(test->ec_k)){
                sprintf(ec_testFileNames[i], "%s%d/%s.process%d.stripe%d", targetDirectoryRoot, (start_ost + i) % test->ec_num_ost, initialTestFileName, rank, i);
            }else{
                sprintf(ec_testFileNames[i], "%s%d/%s.process%d.parity%d", targetDirectoryRoot, (start_ost + i) % test->ec_num_ost, initialTestFileName, rank, i - (test->ec_k));
            }
        }
    }else{
        for(i=0;i<total_stripe_num;i++){
            if(i<(test->ec_k)){
                sprintf(ec_testFileNames[i], "%s%d/%s.stripe%d", targetDirectoryRoot, (start_ost + i) % test->ec_num_ost, initialTestFileName, i);
            }else{
                sprintf(ec_testFileNames[i], "%s%d/%s.parity%d", targetDirectoryRoot, (start_ost + i) % test->ec_num_ost, initialTestFileName, i - (test->ec_k));
            }
        }
    }
    
} /* GetTestFileName() */

void
GetTestFileName(char * testFileName, IOR_param_t * test)
{
    char ** fileNames,
            initialTestFileName[MAXPATHLEN],
            testFileNameRoot[MAX_STR],
            tmpString[MAX_STR];
    int     count;

    /* parse filename for multiple file systems */
    strcpy(initialTestFileName, test->testFileName);
    fileNames = ParseFileName(initialTestFileName, &count);
    if (count > 1 && test->uniqueDir == TRUE)
        ERR("cannot use multiple file names with unique directories");
    if (test->filePerProc) {
        strcpy(testFileNameRoot,
               fileNames[((rank+rankOffset)%test->numTasks) % count]);
    } else {
        strcpy(testFileNameRoot, fileNames[0]);
    }

    /* give unique name if using multiple files */
    if (test->filePerProc) {
        /*
         * prepend rank subdirectory before filename
         * e.g., /dir/file => /dir/<rank>/file
         */
        if (test->uniqueDir == TRUE) {
            strcpy(testFileNameRoot, PrependDir(test, testFileNameRoot));
        }
#if USE_UNDOC_OPT /* NFS */
        if (test->NFS_serverCount) {
            sprintf(tmpString, "%s/%s%d/%s", test->NFS_rootPath,
                    test->NFS_serverName, rank%(test->NFS_serverCount),
                    testFileNameRoot);
            strcpy(testFileNameRoot, tmpString);
        }
#endif /* USE_UNDOC_OPT - NFS */
        sprintf(testFileName, "%s.%08d", testFileNameRoot,
                (rank+rankOffset)%test->numTasks); 
    } else {
        strcpy(testFileName, testFileNameRoot);
    }

    /* add suffix for multiple files */
    if (test->repCounter > -1) {
        sprintf(tmpString, ".%d", test->repCounter);
        strcat(testFileName, tmpString);
    }
} /* GetTestFileName() */


/******************************************************************************/
/*
 * Get time stamp.  Use MPI_Timer() unless _NO_MPI_TIMER is defined,
 * in which case use gettimeofday().
 */

double
GetTimeStamp(void)
{
    double         timeVal;
#ifdef _NO_MPI_TIMER
    struct timeval timer;

    if (gettimeofday(&timer, (struct timezone *)NULL) != 0)
        ERR("cannot use gettimeofday()");
    timeVal = (double)timer.tv_sec + ((double)timer.tv_usec/1000000);
#else /* not _NO_MPI_TIMER */
    timeVal = MPI_Wtime(); /* no MPI_CHECK(), just check return value */
    if (timeVal < 0) ERR("cannot use MPI_Wtime()");
#endif /* _NO_MPI_TIMER */

    /* wall_clock_delta is difference from root node's time */
    timeVal -= wall_clock_delta;

    return(timeVal);
} /* GetTimeStamp() */


/******************************************************************************/
/*
 * Convert IOR_offset_t value to human readable string.
 */

char *
HumanReadable(IOR_offset_t value, int base)
{
    char * valueStr;
    int m = 0, g = 0;
    char m_str[8], g_str[8];

    valueStr = (char *)malloc(MAX_STR);
    if (valueStr == NULL) ERR("out of memory");

    if (base == BASE_TWO) {
        m = MEBIBYTE;
        g = GIBIBYTE;
        strcpy(m_str, "MiB");
        strcpy(g_str, "GiB");
    } else if (base == BASE_TEN) {
        m = MEGABYTE;
        g = GIGABYTE;
        strcpy(m_str, "MB");
        strcpy(g_str, "GB");
    }

    if (value >= g) {
        if (value % (IOR_offset_t)g) {
            sprintf(valueStr, "%.2f %s", (double)((double)value/g), g_str);
        } else {
            sprintf(valueStr, "%d %s", (int)(value/g), g_str);
        }
    } else if (value >= m) {
        if (value % (IOR_offset_t)m) {
            sprintf(valueStr, "%.2f %s", (double)((double)value/m), m_str);
        } else {
            sprintf(valueStr, "%d %s", (int)(value/m), m_str);
        }
    } else if (value >= 0) {
        sprintf(valueStr, "%d bytes", (int)value);
    } else {
        sprintf(valueStr, "-");
    }
    return valueStr;
} /* HumanReadable() */


/******************************************************************************/
/*
 * Change string to lower case.
 */

char *
LowerCase(char * string)
{
    char * nextChar = string;

    while (* nextChar != '\0') {
        * nextChar = (char)tolower((int)* nextChar);
        nextChar++;
    }
    return(string);
} /* LowerCase() */


/******************************************************************************/
/*
 * Parse file name.
 */

char **
ParseFileName(char * name, int * count)
{
    char ** fileNames,
          * tmp,
          * token;
    char    delimiterString[3] = { FILENAME_DELIMITER, '\n', '\0' };
    int     i                  = 0;

    * count = 0;
    tmp = name;

    /* pass one */
    /* if something there, count the first item */
    if (* tmp != '\0') {
        (* count)++;
    }
    /* count the rest of the filenames */
    while (* tmp != '\0') {
        if (* tmp == FILENAME_DELIMITER) {
             (* count)++;
        }
        tmp++;
    }

    fileNames = (char **)malloc((* count) * sizeof(char **));
    if (fileNames == NULL) ERR("out of memory");

    /* pass two */
    token = strtok(name, delimiterString);
    while (token != NULL) {
        fileNames[i] = token;
        token = strtok(NULL, delimiterString);
        i++;
    }
    return(fileNames);
} /* ParseFileName() */


/******************************************************************************/
/*
 * Pretty Print a Double.  The First parameter is a flag determining if left
 * justification should be used.  The third parameter a null-terminated string
 * that should be appended to the number field.
 */

void
PPDouble(int      leftjustify,
         double   number,
         char   * append)
{
    if (number < 0) {
        fprintf(stdout, "   -      %s", append);
    } else {
        if (leftjustify) {
            if (number < 1)
                fprintf(stdout, "%-10.6f%s", number, append);
            else if (number < 3600)
                fprintf(stdout, "%-10.2f%s", number, append);
            else
                fprintf(stdout, "%-10.0f%s", number, append);
        } else {
            if (number < 1)
                fprintf(stdout, "%10.6f%s", number, append);
            else if (number < 3600)
                fprintf(stdout, "%10.2f%s", number, append);
            else
                fprintf(stdout, "%10.0f%s", number, append);
        }
    }
} /* PPDouble() */


/******************************************************************************/
/*
 * From absolute directory, insert rank as subdirectory.  Allows each task
 * to write to its own directory.  E.g., /dir/file => /dir/<rank>/file.
 * 
 */

char *
PrependDir(IOR_param_t * test, char * rootDir)
{
    char * dir;
    char   fname[MAX_STR+1];
    char * p;
    int    i;

    dir = (char *)malloc(MAX_STR+1);
    if (dir == NULL) ERR("out of memory");

    /* get dir name */
    strcpy(dir, rootDir);
    i = strlen(dir)-1;
    while (i > 0) {
        if (dir[i] == '\0' || dir[i] == '/') {
            dir[i] = '/';
            dir[i+1] = '\0';
            break;
        }
        i--;
    }

    /* get file name */
    strcpy(fname, rootDir);
    p = fname;
    while (i > 0) {
        if (fname[i] == '\0' || fname[i] == '/') {
            p = fname + (i+1);
            break;
        }
        i--;
    }

    /* create directory with rank as subdirectory */
        sprintf(dir, "%s%d", dir, (rank+rankOffset)%test->numTasks);

    /* dir doesn't exist, so create */
    if (access(dir, F_OK) != 0) {
        if (mkdir(dir, S_IRWXU) < 0 ) {
            ERR("cannot create directory");
        }

    /* check if correct permissions */
    } else if (access(dir, R_OK) != 0 || access(dir, W_OK) != 0 ||
               access(dir, X_OK) != 0) {
        ERR("invalid directory permissions");
    }

    /* concatenate dir and file names */
    strcat(dir, "/");
    strcat(dir, p);

    return dir;
} /* PrependDir() */


/******************************************************************************//*
 * Read and then reread buffer to confirm data read twice matches.
 */

void
ReadCheck(void         *fd,
          void         *buffer,
          void         *checkBuffer,
          void         *readCheckBuffer,
          IOR_param_t  *test,
          IOR_offset_t  transfer,
          IOR_offset_t  blockSize,
          IOR_offset_t *amtXferred,
          IOR_offset_t *transferCount,
          int           access,
          int          *errors)
{
    int          readCheckToRank, readCheckFromRank;
    MPI_Status   status;
    IOR_offset_t tmpOffset;
    IOR_offset_t segmentSize, segmentNum;

    memset(buffer, 'a', transfer);
    *amtXferred = IOR_Xfer(access, fd, buffer, transfer, test);
    tmpOffset = test->offset;
    if (test->filePerProc == FALSE) {
        /* offset changes for shared file, not for file-per-proc */
        segmentSize = test->numTasks * blockSize;
        segmentNum = test->offset / segmentSize;

                       /* work in current segment */
        test->offset = (((test->offset % segmentSize)
                       /* offset to neighbor's data */
                       + ((test->reorderTasks ?
                       test->tasksPerNode : 0) * blockSize))
                       /* stay within current segment */
                       % segmentSize)
                       /* return segment to actual file offset */
                       + (segmentNum * segmentSize);
    }
    if (*amtXferred != transfer)
        ERR("cannot read from file on read check");
    memset(checkBuffer, 'a', transfer);        /* empty buffer */
#if USE_UNDOC_OPT /* corruptFile */
    MPI_CHECK(MPI_Barrier(testComm), "barrier error");
    /* intentionally corrupt file to determine if check works */
    if (test->corruptFile) {
        CorruptFile(test->testFileName, test, 0, READCHECK);
    }
#endif /* USE_UNDOC_OPT - corruptFile */
    MPI_CHECK(MPI_Barrier(testComm), "barrier error");
    if (test->filePerProc) {
        *amtXferred = IOR_Xfer(access, test->fd_fppReadCheck,
                               checkBuffer, transfer, test);
    } else {
        *amtXferred = IOR_Xfer(access, fd, checkBuffer, transfer, test);
    }
    test->offset = tmpOffset;
    if (*amtXferred != transfer)
        ERR("cannot reread from file read check");
    (*transferCount)++;
    /* exchange buffers */
    memset(readCheckBuffer, 'a', transfer);
    readCheckToRank = (rank + (test->reorderTasks ? test->tasksPerNode : 0))
                      % test->numTasks;
    readCheckFromRank = (rank + (test->numTasks
                        - (test->reorderTasks ? test->tasksPerNode : 0)))
                        % test->numTasks;
    MPI_CHECK(MPI_Barrier(testComm), "barrier error");
    MPI_Sendrecv(checkBuffer, transfer, MPI_CHAR, readCheckToRank, 1,
                 readCheckBuffer, transfer, MPI_CHAR, readCheckFromRank, 1,
                 testComm, &status);
    MPI_CHECK(MPI_Barrier(testComm), "barrier error");
    *errors += CompareBuffers(buffer, readCheckBuffer, transfer,
                              *transferCount, test, READCHECK);
    return;
} /* ReadCheck() */


/******************************************************************************/
/*
 * Reduce test results, and show if verbose set.
 */

void
ReduceIterResults(IOR_param_t  * test,
                  double      ** timer,
                  int            rep,
                  int            access)
{
    double       reduced[12]    = { 0 },
                 diff[6],
                 currentWrite,
                 currentRead,
                 totalWriteTime,
                 totalReadTime;
    enum         {RIGHT, LEFT};
    int          i;
    static int   firstIteration = TRUE;
    IOR_offset_t aggFileSizeForBW;
    MPI_Op       op;

    /* Find the minimum start time of the even numbered timers, and the
       maximum finish time for the odd numbered timers */
    for (i = 0; i < 12; i++) {
        op = i % 2 ? MPI_MAX : MPI_MIN;
        MPI_CHECK(MPI_Reduce(&timer[i][rep], &reduced[i], 1, MPI_DOUBLE,
                             op, 0, testComm),
                  "MPI_Reduce()");
    }

    if (rank == 0) {
        /* Calculate elapsed times and throughput numbers */
        for (i = 0; i < 6; i++)
            diff[i] = reduced[2*i+1] - reduced[2*i];
        totalReadTime  = reduced[11] - reduced[6];
        totalWriteTime = reduced[5] - reduced[0];
        aggFileSizeForBW = test->aggFileSizeForBW[rep];
        if (access == WRITE) {
            currentWrite = (double)((double)aggFileSizeForBW / totalWriteTime);
            test->writeVal[0][rep] = currentWrite;
            test->writeVal[1][rep] = totalWriteTime;
        }
        if (access == READ) {
            currentRead = (double)((double)aggFileSizeForBW / totalReadTime);
            test->readVal[0][rep] = currentRead;
            test->readVal[1][rep] = totalReadTime;
        }
    }


#if USE_UNDOC_OPT /* fillTheFileSystem */
    if (test->fillTheFileSystem && rank == 0 && firstIteration
        && rep == 0 && verbose >= VERBOSE_1) {
        fprintf(stdout, " . . . skipping iteration results output . . .\n");
        fflush(stdout);
    }
#endif /* USE_UNDOC_OPT - fillTheFileSystem */

    if (rank == 0 && verbose >= VERBOSE_2
#if USE_UNDOC_OPT /* fillTheFileSystem */
        && test->fillTheFileSystem == FALSE
#endif /* USE_UNDOC_OPT - fillTheFileSystem */
        ) {
        /* print out the results */
        if (firstIteration && rep == 0) {
            fprintf(stdout, "access    bw(MiB/s)  block(KiB) xfer(KiB)");
            fprintf(stdout, "  open(s)    wr/rd(s)   close(s) total(s)  iter\n");
            fprintf(stdout, "------    ---------  ---------- ---------");
            fprintf(stdout, "  --------   --------   --------  --------   ----\n");
        }
        if (access == WRITE) {
            fprintf(stdout, "write     ");
            PPDouble(LEFT, (currentWrite/MEBIBYTE), " \0");
            PPDouble(LEFT, (double)test->blockSize/KIBIBYTE, " \0");
            PPDouble(LEFT, (double)test->transferSize/KIBIBYTE, " \0");
            if (test->writeFile) {
                PPDouble(LEFT, diff[0], " \0");
                PPDouble(LEFT, diff[1], " \0");
                PPDouble(LEFT, diff[2], " \0");
                PPDouble(LEFT, totalWriteTime, " \0");
            }
            fprintf(stdout, "%-4d XXCEL\n", rep);
        }
        if (access == READ) {
            fprintf(stdout, "read      ");
            PPDouble(LEFT, (currentRead/MEBIBYTE), " \0");
            PPDouble(LEFT, (double)test->blockSize/KIBIBYTE, " \0");
            PPDouble(LEFT, (double)test->transferSize/KIBIBYTE, " \0");
            if (test->readFile) {
                PPDouble(LEFT, diff[3], " \0");
                PPDouble(LEFT, diff[4], " \0");
                PPDouble(LEFT, diff[5], " \0");
                PPDouble(LEFT, totalReadTime, " \0");
            }
            fprintf(stdout, "%-4d XXCEL\n", rep);
        }
        fflush(stdout);
    }
    firstIteration = FALSE;  /* set to TRUE to repeat this header */
} /* ReduceIterResults() */


/******************************************************************************/
/*
 * Check for file(s), then remove all files if file-per-proc, else single file.
 */
     
void
RemoveFile(char        * testFileName,
           int           filePerProc,
           IOR_param_t * test)
{
    int tmpRankOffset;
    if (filePerProc) 
    {
       /* in random tasks, delete own file */
       if (test->reorderTasksRandom == TRUE) 
       {  
         tmpRankOffset = rankOffset;
         rankOffset = 0;
         GetTestFileName(testFileName, test);
       }
       if (access(testFileName, F_OK) == 0)
       {
         IOR_Delete(testFileName, test);
       }
       if (test->reorderTasksRandom == TRUE) 
       {  
         rankOffset = tmpRankOffset;
         GetTestFileName(testFileName, test);
       }
    } else {
        if ((rank == 0) && (access(testFileName, F_OK) == 0)) {
            IOR_Delete(testFileName, test);
        }
    }
} /* RemoveFile() */


/******************************************************************************/
/*
 * Setup tests by parsing commandline and creating test script.
 */

IOR_queue_t *
SetupTests(int argc, char ** argv)
{
    
    //fprintf(stdout, "in SetupTests\n");
    
    
    IOR_queue_t * tests,
                * testsHead;

    /* count the tasks per node */
    tasksPerNode = CountTasksPerNode(numTasksWorld, MPI_COMM_WORLD);

    testsHead = tests = ParseCommandLine(argc, argv);
    /*
     * Since there is no guarantee that anyone other than
     * task 0 has the environment settings for the hints, pass
     * the hint=value pair to everyone else in MPI_COMM_WORLD
     */
    DistributeHints();

    /* check validity of tests and create test queue */
    while (tests != NULL) {
        ValidTests(&tests->testParameters);
        tests = tests->nextTest;
    }

    /* check for skew between tasks' start times */
    wall_clock_deviation = TimeDeviation();

    /* seed random number generator */
    SeedRandGen(MPI_COMM_WORLD);

    return(testsHead);
} /* SetupTests() */


/******************************************************************************//*
 * Setup transfer buffers, creating and filling as needed.
 */

void
SetupXferBuffers(void        **buffer,
                 void        **checkBuffer,
                 void        **readCheckBuffer,
                 IOR_param_t  *test,
                 int           pretendRank,
                 int           access)
{
    /* create buffer of filled data */
    *buffer = CreateBuffer(test->transferSize);
    FillBuffer(*buffer, test, 0, pretendRank);
    if (access == WRITECHECK || access == READCHECK) {
        /* this allocates buffer only */
        *checkBuffer = CreateBuffer(test->transferSize);
        if (access == READCHECK) {
            *readCheckBuffer = CreateBuffer(test->transferSize);
        }
    }
    return;
} /* SetupXferBuffers() */


/******************************************************************************/
/*
 * Print header information for test output.
 */

void
ShowInfo(int argc, char **argv, IOR_param_t *test)
{
    int    i;
    struct utsname unamebuf;

    fprintf(stdout, "%s: MPI Coordinated Test of Parallel I/O\n\n",
            IOR_RELEASE);

    fprintf(stdout, "Run began: %s", CurrentTimeString());
#if USE_UNDOC_OPT /* NFS */
    if (test->NFS_serverCount) {
        fprintf(stdout, "NFS path: %s%s[0..%d]\n", test->NFS_rootPath,
                test->NFS_serverName, test->NFS_serverCount-1);
    }
#endif /* USE_UNDOC_OPT - NFS */
    fprintf(stdout, "Command line used:");
    for (i = 0; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
    }
    fprintf(stdout, "\n");
    if (uname(&unamebuf) != 0) {
        WARN("uname failed");
        fprintf(stdout, "Machine: Unknown");
    } else {
        fprintf(stdout, "Machine: %s %s", unamebuf.sysname, unamebuf.nodename);
        if (verbose >= VERBOSE_2) {
             fprintf(stdout, " %s %s %s", unamebuf.release, unamebuf.version,
                    unamebuf.machine);
        }
    }
    fprintf(stdout, "\n");
#ifdef _NO_MPI_TIMER
    if (verbose >= VERBOSE_2)
        fprintf(stdout, "Using unsynchronized POSIX timer\n");
#else /* not _NO_MPI_TIMER */
    if (MPI_WTIME_IS_GLOBAL) {
        if (verbose >= VERBOSE_2)
            fprintf(stdout, "Using synchronized MPI timer\n");
    } else {
        if (verbose >= VERBOSE_2)
            fprintf(stdout, "Using unsynchronized MPI timer\n");
    }
#endif /* _NO_MPI_TIMER */
    if (verbose >= VERBOSE_1) {
        fprintf(stdout, "Start time skew across all tasks: %.02f sec\n",
            wall_clock_deviation);
        /* if pvfs2:, then skip */
        if (Regex(test->testFileName, "^[a-z][a-z].*:") == 0) {
            DisplayFreespace(test);
        }
    }
    if (verbose >= VERBOSE_3) {       /* show env */
            fprintf(stdout, "STARTING ENVIRON LOOP\n");
        for (i = 0; environ[i] != NULL; i++) {
            fprintf(stdout, "%s\n", environ[i]);
        }
            fprintf(stdout, "ENDING ENVIRON LOOP\n");
    }
    fflush(stdout);
} /* ShowInfo() */


/******************************************************************************/
/*
 * Show simple test output with max results for iterations.
 */

void
ShowSetup(IOR_param_t * test)
{
    IOR_offset_t aggFileSizeForBW;

    /* use expected file size of iteration #0 for display */
    aggFileSizeForBW = test->aggFileSizeFromCalc[0];

    if (strcmp(test->debug, "") != 0) {
        fprintf(stdout, "\n*** DEBUG MODE ***\n");
        fprintf(stdout, "*** %s ***\n\n", test->debug);
    }
    fprintf(stdout, "\nSummary:\n");
    fprintf(stdout, "\tapi                = %s\n",
            test->apiVersion);
    fprintf(stdout, "\ttest filename      = %s\n", test->testFileName);
    fprintf(stdout, "\taccess             = ");
    if (test->filePerProc) {
        fprintf(stdout, "file-per-process");
    } else {
        fprintf(stdout, "single-shared-file");
    }
    if (verbose >= VERBOSE_1 && strcmp(test->api, "POSIX") != 0) {
        if (test->collective == FALSE) {
            fprintf(stdout, ", independent");
        } else {
            fprintf(stdout, ", collective");
        }
    }
    fprintf(stdout, "\n");
    if (verbose >= VERBOSE_1) {
        if (test->segmentCount > 1) {
            fprintf(stdout, "\tpattern            = strided (%d segments)\n",
            (int)test->segmentCount);
        } else {
            fprintf(stdout, "\tpattern            = segmented (1 segment)\n");
        }
    }
    fprintf(stdout, "\tordering in a file =");
    if (test->randomOffset == FALSE) {
        fprintf(stdout, " sequential offsets\n");
    } else {
        fprintf(stdout, " random offsets\n");
    }
    fprintf(stdout, "\tordering inter file=");
    if (test->reorderTasks == FALSE && test->reorderTasksRandom == FALSE){
        fprintf(stdout, " no tasks offsets\n");
    } 
    if (test->reorderTasks == TRUE)
    {
        fprintf(stdout, "constant task offsets = %d\n",test->taskPerNodeOffset);
    }
    if (test->reorderTasksRandom == TRUE) 
    {
        fprintf(stdout, "random task offsets >= %d, seed=%d\n",test->taskPerNodeOffset,test->reorderTasksRandomSeed);
    }
    fprintf(stdout, "\tclients            = %d (%d per node)\n",
            test->numTasks, test->tasksPerNode);
    fprintf(stdout, "\trepetitions        = %d\n",
            test->repetitions);
    fprintf(stdout, "\txfersize           = %s\n",
            HumanReadable(test->transferSize, BASE_TWO));
    fprintf(stdout, "\tblocksize          = %s\n",
            HumanReadable(test->blockSize, BASE_TWO));
    fprintf(stdout, "\taggregate filesize = %s\n",
            HumanReadable(aggFileSizeForBW, BASE_TWO));
#ifdef _MANUALLY_SET_LUSTRE_STRIPING
    fprintf(stdout, "\tLustre stripe size = %s\n",
            ((test->lustre_stripe_size == 0) ? "Use default" :
            HumanReadable(test->lustre_stripe_size, BASE_TWO)));
    if (test->lustre_stripe_count == 0) {
        fprintf(stdout, "\t      stripe count = %s\n", "Use default");
    } else {
        fprintf(stdout, "\t      stripe count = %d\n",
                test->lustre_stripe_count);
    }
#endif /* _MANUALLY_SET_LUSTRE_STRIPING */
    if (test->deadlineForStonewalling > 0) {
        fprintf(stdout, "\tUsing stonewalling = %d second(s)\n",
                test->deadlineForStonewalling);
    }
    fprintf(stdout, "\n");
    /*fprintf(stdout, "\n  ================================\n\n");*/
    fflush(stdout);
} /* ShowSetup() */


/******************************************************************************/
/*
 * Show test description.
 */

void
ShowTest(IOR_param_t * test)
{
    fprintf(stdout, "\n\n  ================================\n\n");
    fprintf(stdout, "TEST:\t%s=%d\n", "id", test->id);
    fprintf(stdout, "\t%s=%d\n", "testnum", test->TestNum);
    fprintf(stdout, "\t%s=%s\n", "api", test->api);
    fprintf(stdout, "\t%s=%s\n", "platform", test->platform);
    fprintf(stdout, "\t%s=%s\n", "testFileName", test->testFileName);
    fprintf(stdout, "\t%s=%s\n", "hintsFileName", test->hintsFileName);
    fprintf(stdout, "\t%s=%d\n", "deadlineForStonewall",
            test->deadlineForStonewalling);
    fprintf(stdout, "\t%s=%d\n", "maxTimeDuration", test->maxTimeDuration);
    fprintf(stdout, "\t%s=%d\n", "outlierThreshold", test->outlierThreshold);
    fprintf(stdout, "\t%s=%s\n", "options", test->options);
    fprintf(stdout, "\t%s=%d\n", "nodes", test->nodes);
    fprintf(stdout, "\t%s=%d\n", "tasksPerNode", tasksPerNode);
    fprintf(stdout, "\t%s=%d\n", "repetitions", test->repetitions);
    fprintf(stdout, "\t%s=%d\n", "multiFile", test->multiFile);
    fprintf(stdout, "\t%s=%d\n", "interTestDelay", test->interTestDelay);
    fprintf(stdout, "\t%s=%d\n", "fsync", test->fsync);
    fprintf(stdout, "\t%s=%d\n", "fsYncperwrite", test->fsyncPerWrite);
    fprintf(stdout, "\t%s=%d\n", "useExistingTestFile",
            test->useExistingTestFile);
    fprintf(stdout, "\t%s=%d\n", "showHints", test->showHints);
    fprintf(stdout, "\t%s=%d\n", "uniqueDir", test->uniqueDir);
    fprintf(stdout, "\t%s=%d\n", "showHelp", test->showHelp);
    fprintf(stdout, "\t%s=%d\n", "individualDataSets",test->individualDataSets);
    fprintf(stdout, "\t%s=%d\n", "singleXferAttempt", test->singleXferAttempt);
    fprintf(stdout, "\t%s=%d\n", "readFile", test->readFile);
    fprintf(stdout, "\t%s=%d\n", "writeFile", test->writeFile);
    fprintf(stdout, "\t%s=%d\n", "filePerProc", test->filePerProc);
    fprintf(stdout, "\t%s=%d\n", "reorderTasks", test->reorderTasks);
    fprintf(stdout, "\t%s=%d\n", "reorderTasksRandom", test->reorderTasksRandom);
    fprintf(stdout, "\t%s=%d\n", "reorderTasksRandomSeed", test->reorderTasksRandomSeed);
    fprintf(stdout, "\t%s=%d\n", "randomOffset", test->randomOffset);
    fprintf(stdout, "\t%s=%d\n", "checkWrite", test->checkWrite);
    fprintf(stdout, "\t%s=%d\n", "checkRead", test->checkRead);
    fprintf(stdout, "\t%s=%d\n", "preallocate", test->preallocate);
    fprintf(stdout, "\t%s=%d\n", "useFileView", test->useFileView);
    fprintf(stdout, "\t%s=%lld\n", "setAlignment", test->setAlignment);
    fprintf(stdout, "\t%s=%d\n", "storeFileOffset", test->storeFileOffset);
    fprintf(stdout, "\t%s=%d\n", "useSharedFilePointer",
            test->useSharedFilePointer);
    fprintf(stdout, "\t%s=%d\n", "useO_DIRECT", test->useO_DIRECT);
    fprintf(stdout, "\t%s=%d\n", "useStridedDatatype",
            test->useStridedDatatype);
    fprintf(stdout, "\t%s=%d\n", "keepFile", test->keepFile);
    fprintf(stdout, "\t%s=%d\n", "keepFileWithError", test->keepFileWithError);
    fprintf(stdout, "\t%s=%d\n", "quitOnError", test->quitOnError);
    fprintf(stdout, "\t%s=%d\n", "verbose", verbose);
    fprintf(stdout, "\t%s=%d\n", "setTimeStampSignature",
            test->setTimeStampSignature);
    fprintf(stdout, "\t%s=%d\n", "collective", test->collective);
    fprintf(stdout, "\t%s=%lld", "segmentCount", test->segmentCount);
    if (strcmp(test->api, "HDF5") == 0) {
        fprintf(stdout, " (datasets)");
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "\t%s=%lld\n", "transferSize", test->transferSize);
    fprintf(stdout, "\t%s=%lld\n", "blockSize", test->blockSize);
} /* ShowTest() */


/******************************************************************************/
/*
 * Summarize results, showing max rates (and min, mean, stddev if verbose)
 */

void
SummarizeResults(IOR_param_t * test)
{
    int    rep, ival;
    double maxWrite[2], minWrite[2], maxRead[2], minRead[2],
           meanWrite[2], meanRead[2],
           varWrite[2], varRead[2],
           sdWrite[2], sdRead[2],
           sumWrite[2], sumRead[2],
           writeTimeSum, readTimeSum;

   for (ival=0; ival<2; ival++)
   {
      varWrite[ival] = 0; 
      varRead[ival] = 0;
      sdWrite[ival] = 0; 
      sdRead[ival] = 0;
      sumWrite[ival] = 0;
      sumRead[ival] = 0;
      writeTimeSum = 0;
      readTimeSum = 0;

      if (ival == 1)
      {
         for (rep = 0; rep < test->repetitions; rep++) 
        {
         writeTimeSum += test->writeVal[ival][rep];
         readTimeSum  += test->readVal [ival][rep];
         test->writeVal[ival][rep] = (double)test->numTasks*((double)test->blockSize/(double)test->transferSize)/test->writeVal[ival][rep];
         test->readVal [ival][rep] = (double)test->numTasks*((double)test->blockSize/(double)test->transferSize)/test->readVal [ival][rep];
        }
      }

    maxWrite[ival] = minWrite[ival] = test->writeVal[ival][0];
    maxRead[ival] = minRead[ival] = test->readVal[ival][0];

    for (rep = 0; rep < test->repetitions; rep++) {
        if (maxWrite[ival] < test->writeVal[ival][rep]) {
            maxWrite[ival] = test->writeVal[ival][rep];
        }
        if (maxRead[ival] < test->readVal[ival][rep]) {
            maxRead[ival] = test->readVal[ival][rep];
        }
        if (minWrite[ival] > test->writeVal[ival][rep]) {
            minWrite[ival] = test->writeVal[ival][rep];
        }
        if (minRead[ival] > test->readVal[ival][rep]) {
            minRead[ival] = test->readVal[ival][rep];
        }
        sumWrite[ival] += test->writeVal[ival][rep];
        sumRead[ival] += test->readVal[ival][rep];
    }

    meanWrite[ival] = sumWrite[ival] / test->repetitions;
    meanRead[ival] = sumRead[ival] / test->repetitions;

    for (rep = 0; rep < test->repetitions; rep++) {
        varWrite[ival] += pow((meanWrite[ival] - test->writeVal[ival][rep]), 2);
        varRead[ival] += pow((meanRead[ival] - test->readVal[ival][rep]), 2);
    }
    varWrite[ival] = varWrite[ival] / test->repetitions;
    varRead[ival] = varRead[ival] / test->repetitions;
    sdWrite[ival] = sqrt(varWrite[ival]);
    sdRead[ival] = sqrt(varRead[ival]);
   }

    if (rank == 0 && verbose >= VERBOSE_0) 
    {
        fprintf(stdout, "Operation  Max (MiB)  Min (MiB)  Mean (MiB)   Std Dev  Max (OPs)  Min (OPs)  Mean (OPs)   Std Dev  Mean (s)  ");
        if (verbose >= VERBOSE_1)
          fprintf(stdout, "Op grep #Tasks tPN reps  fPP reord reordoff reordrand seed segcnt blksiz xsize aggsize\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "---------  ---------  ---------  ----------   -------  ---------  ---------  ----------   -------  --------\n");
        if (maxWrite[0] > 0.)
        {
            fprintf(stdout,"%s     ", "write");
            fprintf(stdout,"%10.2f ", maxWrite[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f  ", minWrite[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f", meanWrite[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f ", sdWrite[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f ", maxWrite[1]);
            fprintf(stdout,"%10.2f  ", minWrite[1]);
            fprintf(stdout,"%10.2f", meanWrite[1]);
            fprintf(stdout,"%10.2f", sdWrite[1]);
            fprintf(stdout,"%10.5f   ", writeTimeSum/(double)(test->repetitions));
            if (verbose >= VERBOSE_1)
            {
                fprintf(stdout,"%d ", test->numTasks);
                fprintf(stdout,"%d ", test->tasksPerNode);
                fprintf(stdout,"%d ", test->repetitions);
                fprintf(stdout,"%d ", test->filePerProc);
                fprintf(stdout,"%d ", test->reorderTasks);
                fprintf(stdout,"%d ", test->taskPerNodeOffset);
                fprintf(stdout,"%d ", test->reorderTasksRandom);
                fprintf(stdout,"%d ", test->reorderTasksRandomSeed);
                fprintf(stdout,"%lld ", test->segmentCount);
                fprintf(stdout,"%lld ", test->blockSize);
                fprintf(stdout,"%lld ", test->transferSize);
                fprintf(stdout,"%lld ", test->aggFileSizeForBW[0]);
                fprintf(stdout,"%d ", test->TestNum);
                fprintf(stdout,"%s ", test->api);
            }
            fprintf(stdout,"EXCEL\n");
        }
        if (maxRead[0] > 0.)
        {
            fprintf(stdout,"%s      ", "read");
            fprintf(stdout,"%10.2f ", maxRead[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f  ", minRead[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f", meanRead[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f ", sdRead[0]/MEBIBYTE);
            fprintf(stdout,"%10.2f ", maxRead[1]);
            fprintf(stdout,"%10.2f  ", minRead[1]);
            fprintf(stdout,"%10.2f", meanRead[1]);
            fprintf(stdout,"%10.2f", sdRead[1]);
            fprintf(stdout,"%10.5f   ", readTimeSum/(double)(test->repetitions));
            if (verbose >= VERBOSE_1)
            {
                fprintf(stdout,"%d ", test->numTasks);
                fprintf(stdout,"%d ", test->tasksPerNode);
                fprintf(stdout,"%d ", test->repetitions);
                fprintf(stdout,"%d ", test->filePerProc);
                fprintf(stdout,"%d ", test->reorderTasks);
                fprintf(stdout,"%d ", test->taskPerNodeOffset);
                fprintf(stdout,"%d ", test->reorderTasksRandom);
                fprintf(stdout,"%d ", test->reorderTasksRandomSeed);
                fprintf(stdout,"%lld ", test->segmentCount);
                fprintf(stdout,"%lld ", test->blockSize);
                fprintf(stdout,"%lld ", test->transferSize);
                fprintf(stdout,"%lld ", test->aggFileSizeForBW[0]);
                fprintf(stdout,"%d ", test->TestNum);
                fprintf(stdout,"%s ", test->api);
            }
            fprintf(stdout,"EXCEL\n");
        }
        fflush(stdout);
   }

    if (rank == 0 && verbose >= VERBOSE_0) {
        fprintf(stdout, "\n");
        if (test->writeFile) {
            fprintf(stdout, "Max Write: %.2f MiB/sec (%.2f MB/sec)\n",
                    maxWrite[0]/MEBIBYTE, maxWrite[0]/MEGABYTE);
        }
        if (test->readFile) {
            fprintf(stdout, "Max Read:  %.2f MiB/sec (%.2f MB/sec)\n",
                    maxRead[0]/MEBIBYTE, maxRead[0]/MEGABYTE);
        }
        fprintf(stdout, "\n");
    }
} /* SummarizeResults() */


/******************************************************************************/
/*
 * Using the test parameters, run iteration(s) of single test.
 */

void
TestIoSys(IOR_param_t *test)
{
    // if(test->ec_verbose >= VERBOSE_3){
    //     fprintf(stdout, "in TestIoSys\n");
    // }
    
    char           testFileName[MAX_STR];
    double       * timer[12];
    double         startTime;
    int            i,
                   rep,
                   maxTimeDuration;
    void         * fd;//origin mark
    /********************ec fd*******************/
    char        ** ec_testFileNames;
    void        ** ec_fds;
    /********************ec fd*******************/

    MPI_Group      orig_group, new_group;
    int            range[3];
    IOR_offset_t   dataMoved;             /* for data rate calculation */

    /*ec dataMoved*/
    // IOR_offset_t ec_dataMoved1;
    // IOR_offset_t ec_dataMoved2;
    // IOR_offset_t ec_dataMoved3;
    // IOR_offset_t ec_dataMoved4;
    /*space for ec-structrue*/
    //ec_info my_ec_info;
    /*space for ec-structrue*/

    /* set up communicator for test */
    if (test->numTasks > numTasksWorld) {
        if (rank == 0) {
            fprintf(stdout,
                    "WARNING: More tasks requested (%d) than available (%d),",
                    test->numTasks, numTasksWorld);
            fprintf(stdout,"         running on %d tasks.\n", numTasksWorld);
        }
        test->numTasks = numTasksWorld;
    }
    MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &orig_group),
              "MPI_Comm_group() error");
    range[0] = 0; /* first rank */
    range[1] = test->numTasks - 1; /* last rank */
    range[2] = 1; /* stride */
    MPI_CHECK(MPI_Group_range_incl(orig_group, 1, &range, &new_group),
              "MPI_Group_range_incl() error");
    MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, new_group, &testComm),
              "MPI_Comm_create() error");
    test->testComm = testComm;
    if (testComm == MPI_COMM_NULL) {
        /* tasks not in the group do not participate in this test */
        MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD), "barrier error");
        return;
    }
    if (rank == 0 && verbose >= VERBOSE_1) {
        fprintf(stdout, "Participating tasks: %d\n", test->numTasks);
        fflush(stdout);
    }
    if (rank == 0 && test->reorderTasks == TRUE && verbose >= VERBOSE_1) {
        fprintf(stdout,
    "Using reorderTasks '-C' (expecting block, not cyclic, task assignment)\n");
        fflush(stdout);
    }
    test->tasksPerNode = CountTasksPerNode(test->numTasks, testComm);

    /* setup timers */
    for (i = 0; i < 12; i++) {
        timer[i] = (double *)malloc(test->repetitions * sizeof(double));
        if (timer[i] == NULL) ERR("out of memory");
    }
    for (i = 0; i < 2; i++) 
    {   
    test->writeVal[i] = (double *)malloc(test->repetitions * sizeof(double));
    if (test->writeVal[i] == NULL) ERR("out of memory");
    test->readVal[i] = (double *)malloc(test->repetitions * sizeof(double));
    if (test->readVal[i] == NULL) ERR("out of memory");
    for (rep = 0; rep < test->repetitions; rep++) {
        test->writeVal[i][rep] = test->readVal[i][rep] = 0.0;
    }
    }
    test->aggFileSizeFromCalc =
        (IOR_offset_t *)malloc(test->repetitions * sizeof(IOR_offset_t));
    if (test->aggFileSizeFromCalc == NULL) ERR("out of memory");
    test->aggFileSizeFromStat =
        (IOR_offset_t *)malloc(test->repetitions * sizeof(IOR_offset_t));
    if (test->aggFileSizeFromStat == NULL) ERR("out of memory");
    test->aggFileSizeFromXfer =
        (IOR_offset_t *)malloc(test->repetitions * sizeof(IOR_offset_t));
    if (test->aggFileSizeFromXfer == NULL) ERR("out of memory");
    test->aggFileSizeForBW =
        (IOR_offset_t *)malloc(test->repetitions * sizeof(IOR_offset_t));
    if (test->aggFileSizeForBW == NULL) ERR("out of memory");

    /* bind I/O calls to specific API */
    AioriBind(test->api);

    /* file size is first calculated on for comparison */
    for (rep = 0; rep < test->repetitions; rep++) {
        test->aggFileSizeFromCalc[rep] = test->blockSize * test->segmentCount *
                                         test->numTasks;
    }

    /* show test setup */
    if (rank == 0 && verbose >= VERBOSE_0) ShowSetup(test);
    //fprintf(stdout, "break at 2000\n");

    startTime = GetTimeStamp();
    maxTimeDuration = test->maxTimeDuration * 60;	/* convert to seconds */

#if USE_UNDOC_OPT /* fillTheFileSystem */
    if (rank == 0 && test->fillTheFileSystem && verbose >= VERBOSE_0) {
        fprintf(stdout, "Run started: %s", CurrentTimeString());
    }
#endif /* USE_UNDOC_OPT - fillTheFileSystem */

    /* loop over test iterations */
    for (rep = 0; rep < test->repetitions; rep++) {

        /* Get iteration start time in seconds in task 0 and broadcast to
           all tasks */
        if (rank == 0) {
            if (test->setTimeStampSignature) {
                test->timeStampSignatureValue =
                    (unsigned int)test->setTimeStampSignature;
            } else {
                time_t currentTime;
                if ((currentTime = time(NULL)) == -1) {
                    ERR("cannot get current time");
                }
                test->timeStampSignatureValue = (unsigned int)currentTime;
            }
            if (verbose >= VERBOSE_2) {
                fprintf(stdout,
                    "Using Time Stamp %u (0x%x) for Data Signature\n",
                    test->timeStampSignatureValue,
                    test->timeStampSignatureValue);
            }
        }
        MPI_CHECK(MPI_Bcast(&test->timeStampSignatureValue, 1, MPI_UNSIGNED, 0,
                  testComm), "cannot broadcast start time value");
#if USE_UNDOC_OPT /* fillTheFileSystem */
        if (test->fillTheFileSystem &&
            rep > 0 && rep % (test->fillTheFileSystem / test->numTasks) == 0) {
            if (rank == 0 && verbose >= VERBOSE_0) {
                fprintf(stdout, "at file #%d, time: %s", rep*test->numTasks,
                        CurrentTimeString());
                fflush(stdout);
            }
        }
#endif /* USE_UNDOC_OPT - fillTheFileSystem */
        /* use repetition count for number of multiple files */
        if (test->multiFile) test->repCounter = rep;

        /*
         * write the file(s), getting timing between I/O calls
         */

        if (test->writeFile
#if USE_UNDOC_OPT /* multiReRead */
          && (!test->multiReRead || !rep)
#endif /* USE_UNDOC_OPT - multiReRead */
          && (maxTimeDuration
              ? (GetTimeStamp() - startTime < maxTimeDuration) : 1)) {
            GetTestFileName(testFileName, test);//origin_mark

            /*****************ec_file_init*************************/
            if(test->ec_verbose >= VERBOSE_3){
                fprintf(stdout, "in ec_file_init\n");
            }  
            int total_stripe_num = test->ec_k + test->ec_m;
            ec_testFileNames = (char **)malloc(sizeof(char *) * total_stripe_num);
            //fprintf(stdout, "break at 2067\n");
            for(i = 0;i<total_stripe_num;i++){
                ec_testFileNames[i] = (char *)malloc(sizeof(char) * MAX_STR);
            }
            GetTestFileName_ec(ec_testFileNames, test);
            //fprintf(stdout, "break at 2072\n");
            if(test->ec_verbose >= VERBOSE_1){
                //fprintf(stdout, "break at 2074\n");
                fprintf(stdout, "process %d's target files initialized as:\n", rank);
                for(i = 0;i<total_stripe_num;i++){
                    fprintf(stdout, "%s\n", ec_testFileNames[i]);
                }
            }           
            /*****************ec_file_init*************************/


            if (verbose >= VERBOSE_3) {
                fprintf(stdout, "task %d writing %s\n", rank, testFileName);
            }
            DelaySecs(test->interTestDelay);
            if (test->useExistingTestFile == FALSE) {
                RemoveFile(testFileName, test->filePerProc, test);//origin_mark
                /*******************delete ec files********************/
                for(i=0;i<total_stripe_num;i++){
                    if (access(ec_testFileNames[i], F_OK) == 0){
                        IOR_Delete(ec_testFileNames[i], test);
                    }
                }
                /*******************delete ec files********************/
            }
            //fprintf(stdout, "break at 2096\n");
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            test->open = WRITE;
            timer[0][rep] = GetTimeStamp();
            fd = IOR_Create(testFileName, test);//origin_mark

            /************************create ec fds**********************/
            ec_fds = (void **)malloc(sizeof(void *) * total_stripe_num);
            //fprintf(stdout, "break at 2100\n");
            for(i=0;i<total_stripe_num;i++){
                ec_fds[i] = IOR_Create(ec_testFileNames[i], test);
                //fprintf(stdout, "break at 2103\n");
                if(ec_fds[i] == NULL){
                    ERR("open ec fds failed");
                }
            }
            /************************create ec fds**********************/

            timer[1][rep] = GetTimeStamp();
            if (test->intraTestBarriers)
                MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            if (rank == 0 && verbose >= VERBOSE_1) {
                fprintf(stdout, "Commencing write performance test.\n");
                fprintf(stdout, "%s\n", CurrentTimeString());
            }
            timer[2][rep] = GetTimeStamp();
            //dataMoved = WriteOrRead(test, fd, WRITE); //origin_mark
            // if(rank == (numTasksWorld-1)){
            //     dataMoved = WriteOrRead_CL(test, ec_fds, WRITE);
            // }else{
            //     dataMoved = WriteOrRead_ec(test, ec_fds, WRITE);//ec_version
            // }
            dataMoved = WriteOrRead_ec(test, ec_fds, WRITE);//ec_version
            //fprintf(stdout, "break at 2120\n");
            //fprintf(stdout, "datamoved=%lld\n", dataMoved);
            timer[3][rep] = GetTimeStamp();
            if (test->intraTestBarriers)
                MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            timer[4][rep] = GetTimeStamp();
            //IOR_Close(ec_fds.fd1, test);//origin_mark
            /*************************ec close***************************/
            for (i = 0; i < total_stripe_num; i++)
            {
                IOR_Close(ec_fds[i], test);
            }
            /*************************ec close***************************/

#if USE_UNDOC_OPT /* includeDeleteTime */
            if (test->includeDeleteTime) {
                if (rank == 0 && verbose >= VERBOSE_1) {
                    fprintf(stdout, "** including delete time **\n");
                }
                MPI_CHECK(MPI_Barrier(testComm), "barrier error");
                RemoveFile(testFileName, test->filePerProc, test);
            }
#endif /* USE_UNDOC_OPT - includeDeleteTime */

            timer[5][rep] = GetTimeStamp();
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");

#if USE_UNDOC_OPT /* includeDeleteTime */
            if (test->includeDeleteTime) {
                /* not accurate, but no longer a test file to examine */
                test->aggFileSizeFromStat[rep] = test->aggFileSizeFromCalc[rep];
                test->aggFileSizeFromXfer[rep] = test->aggFileSizeFromCalc[rep];
                test->aggFileSizeForBW[rep] = test->aggFileSizeFromCalc[rep];
            } else {
#endif /* USE_UNDOC_OPT - includeDeleteTime */

            /* get the size of the file just written */
            test->aggFileSizeFromStat[rep]
                = IOR_GetFileSize(test, testComm, testFileName);

            /* check if stat() of file doesn't equal expected file size,
               use actual amount of byte moved */
            CheckFileSize(test, dataMoved, rep);

#if USE_UNDOC_OPT /* includeDeleteTime */
            }
#endif /* USE_UNDOC_OPT - includeDeleteTime */

            if (verbose >= VERBOSE_3) WriteTimes(test, timer, rep, WRITE);
            ReduceIterResults(test, timer, rep, WRITE);
            if (test->outlierThreshold) {
                CheckForOutliers(test, timer, rep, WRITE);
            }
        }

        /*
         * perform a check of data, reading back data and comparing
         * against what was expected to be written
         */
        if (test->checkWrite
#if USE_UNDOC_OPT /* multiReRead */
          && (!test->multiReRead || rep)
#endif /* USE_UNDOC_OPT - multiReRead */
          && (maxTimeDuration
              ? (GetTimeStamp() - startTime < maxTimeDuration) : 1)) {
#if USE_UNDOC_OPT /* corruptFile */
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            /* intentionally corrupt file to determine if check works */
            if (test->corruptFile) {
                CorruptFile(testFileName, test, rep, WRITECHECK);
            }
#endif /* USE_UNDOC_OPT - corruptFile */
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            if (rank == 0 && verbose >= VERBOSE_1) {
                fprintf(stdout,
                        "Verifying contents of the file(s) just written.\n");
                fprintf(stdout, "%s\n", CurrentTimeString());
            }
            if (test->reorderTasks) {
                /* move two nodes away from writing node */
                rankOffset = (2*test->tasksPerNode) % test->numTasks;
            }
            GetTestFileName(testFileName, test);
            test->open = WRITECHECK;
            fd = IOR_Open(testFileName, test);
            dataMoved = WriteOrRead(test, fd, WRITECHECK);
            IOR_Close(fd, test);
            rankOffset = 0;
        }
        /*
         * read the file(s), getting timing between I/O calls
         */
        if (test->readFile
          && (maxTimeDuration ? (GetTimeStamp() - startTime < maxTimeDuration) : 1)) {
            /* Get rankOffset [file offset] for this process to read, based on -C,-Z,-Q,-X options */
            /* Constant process offset reading */
            if (test->reorderTasks) {
                /* move taskPerNodeOffset nodes[1==default] away from writing node */
                rankOffset = (test->taskPerNodeOffset*test->tasksPerNode) % test->numTasks;
            }
            /* random process offset reading */
            if (test->reorderTasksRandom) 
            {
               /* this should not intefere with randomOffset within a file because GetOffsetArrayRandom */
               /* seeds every random() call  */
               int *rankoffs, *filecont, *filehits, ifile, jfile, nodeoffset;
               unsigned int iseed0;
               nodeoffset = test->taskPerNodeOffset;
               nodeoffset = (nodeoffset < test->nodes) ? nodeoffset : test->nodes-1;
               iseed0 = (test->reorderTasksRandomSeed < 0) ? (-1*test->reorderTasksRandomSeed+rep):test->reorderTasksRandomSeed;
               srand(rank+iseed0); 
                                                                    { rankOffset = rand() % test->numTasks; }
               while (rankOffset < (nodeoffset*test->tasksPerNode)) { rankOffset = rand() % test->numTasks; }
               /* Get more detailed stats if requested by verbose level */
               if (verbose >= VERBOSE_2)
               {
                 if (rank == 0)
                 {
                   rankoffs = (int *)malloc(test->numTasks*sizeof(int));
                   filecont = (int *)malloc(test->numTasks*sizeof(int));
                   filehits = (int *)malloc(test->numTasks*sizeof(int));
                 }
                 MPI_CHECK(MPI_Gather(&rankOffset,1,MPI_INT,rankoffs,1,MPI_INT,0,MPI_COMM_WORLD), "MPI_Gather error");
                 /*file hits histogram*/
                 if (rank == 0)
                 {
                   memset((void*)filecont,0,test->numTasks*sizeof(int));
                   for (ifile=0; ifile<test->numTasks; ifile++) { filecont[(ifile+rankoffs[ifile])%test->numTasks]++; }
                   memset((void*)filehits,0,test->numTasks*sizeof(int));
                   for (ifile=0; ifile<test->numTasks; ifile++) 
                    for (jfile=0; jfile<test->numTasks; jfile++) { if (ifile == filecont[jfile]) filehits[ifile]++; }
                   /*                                             fprintf(stdout, "File Contention Dist:");
                   for (ifile=0; ifile<test->numTasks; ifile++) { fprintf(stdout," %d",filecont[ifile]); }
                                                                  fprintf(stdout,"\n"); */
                                                                  fprintf(stdout, "#File Hits Dist:");
                          jfile = 0; ifile = 0;
                   while (jfile<test->numTasks && 
                          ifile<test->numTasks) { fprintf(stdout," %d",filehits[ifile]);jfile+=filehits[ifile],ifile++;}
                                                                  fprintf(stdout," XXCEL\n");
                   free(rankoffs);
                   free(filecont);
                   free(filehits);
                 }
               }
            }
            /* Using globally passed rankOffset, following function generates testFileName to read */
            GetTestFileName(testFileName, test);

            /*****************ec_file_init*************************/
            int total_stripe_num = test->ec_k + test->ec_m;
            ec_testFileNames = (char **)malloc(sizeof(char *) * total_stripe_num);
            for (i = 0; i < total_stripe_num; i++)
            {
                ec_testFileNames[i] = (char *)malloc(sizeof(char) * MAX_STR);
            }
            GetTestFileName_ec(ec_testFileNames, test);
            if (test->ec_verbose >= VERBOSE_1)
            {
                fprintf(stdout, "process %d will read following files:\n", rank);
                for (i = 0; i < total_stripe_num; i++)
                {
                    fprintf(stdout, "%s\n", ec_testFileNames[i]);
                }
            }
            /*****************ec_file_init*************************/

            if (verbose >= VERBOSE_3) {
                fprintf(stdout, "task %d reading %s\n", rank, testFileName);
            }
            DelaySecs(test->interTestDelay);
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            test->open = READ;
            timer[6][rep] = GetTimeStamp();
            fd = IOR_Open(testFileName, test);
            if (rank == 0 && verbose >= VERBOSE_2) {
                fprintf(stdout, "[RANK %03d] open for reading file %s XXCEL\n", rank,testFileName);
            }

            /************************create ec fds**********************/
            ec_fds = (char **)malloc(sizeof(char *) * total_stripe_num);
            for (i = 0; i < total_stripe_num; i++)
            {
                ec_fds[i] = IOR_Create(ec_testFileNames[i], test);
                if (ec_fds[i] == NULL)
                {
                    ERR("open ec fds failed");
                }
            }
            /************************create ec fds**********************/

            timer[7][rep] = GetTimeStamp();
            if (test->intraTestBarriers)
                MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            if (rank == 0 && verbose >= VERBOSE_1) {
                fprintf(stdout, "Commencing read performance test.\n");
                fprintf(stdout, "%s\n", CurrentTimeString());
            }
            timer[8][rep] = GetTimeStamp();
            //dataMoved = WriteOrRead(test, fd, READ);
            // if(rank==(numTasksWorld-1)){
            //     dataMoved = WriteOrRead_CL(test, ec_fds, READ);
            // }else{
            dataMoved = WriteOrRead_ec(test, ec_fds, READ);
            // }
            
            timer[9][rep] = GetTimeStamp();
            if (test->intraTestBarriers)
                MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            timer[10][rep] = GetTimeStamp();
            IOR_Close(fd, test);

            /*************************ec close***************************/
            for (i = 0; i < total_stripe_num; i++)
            {
                IOR_Close(ec_fds[i], test);
            }
            /*************************ec close***************************/

            timer[11][rep] = GetTimeStamp();

            /* get the size of the file just read */
            test->aggFileSizeFromStat[rep] = IOR_GetFileSize(test, testComm,
                                                             testFileName);

            /* check if stat() of file doesn't equal expected file size,
               use actual amount of byte moved */
            CheckFileSize(test, dataMoved, rep);

            if (verbose >= VERBOSE_3) WriteTimes(test, timer, rep, READ);
            ReduceIterResults(test, timer, rep, READ);
            if (test->outlierThreshold) {
                CheckForOutliers(test, timer, rep, READ);
            }
        } /* end readFile test */

        /*
         * perform a check of data, reading back data twice and
         * comparing against what was expected to be read
         */
        if (test->checkRead
          && (maxTimeDuration
          ? (GetTimeStamp() - startTime < maxTimeDuration) : 1)) {
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            if (rank == 0 && verbose >= VERBOSE_1) {
                fprintf(stdout, "Re-reading the file(s) twice to ");
                fprintf(stdout, "verify that reads are consistent.\n");
                fprintf(stdout, "%s\n", CurrentTimeString());
            }
            if (test->reorderTasks) {
                /* move three nodes away from reading node */
                rankOffset = (3*test->tasksPerNode) % test->numTasks;
            }
            GetTestFileName(testFileName, test);
            MPI_CHECK(MPI_Barrier(testComm), "barrier error");
            test->open = READCHECK;
            fd = IOR_Open(testFileName, test);
            if (test->filePerProc) {
                int tmpRankOffset;
                tmpRankOffset = rankOffset;
                /* increment rankOffset to open comparison file on other node */
                if (test->reorderTasks) {
                    /* move four nodes away from reading node */
                    rankOffset = (4*test->tasksPerNode) % test->numTasks;
                }
                GetTestFileName(test->testFileName_fppReadCheck, test);
                rankOffset = tmpRankOffset;
                test->fd_fppReadCheck =
                    IOR_Open(test->testFileName_fppReadCheck, test);
            }
            dataMoved = WriteOrRead(test, fd, READCHECK);
            if (test->filePerProc) {
                IOR_Close(test->fd_fppReadCheck, test);
                test->fd_fppReadCheck = NULL;
            }
            IOR_Close(fd, test);
        }
        /*
         * this final barrier may not be necessary as IOR_Close should
         * be a collective call -- but to make sure that the file has
         * has not be removed by a task before another finishes writing,
         * the MPI_Barrier() call has been included.
         */
        MPI_CHECK(MPI_Barrier(testComm), "barrier error");
        if (!test->keepFile && !(test->keepFileWithError && test->errorFound)) {
            RemoveFile(testFileName, test->filePerProc, test);
        }
        test->errorFound = FALSE;
        MPI_CHECK(MPI_Barrier(testComm), "barrier error");
        rankOffset = 0;
    }

#if USE_UNDOC_OPT /* fillTheFileSystem */
    if (rank == 0 && test->fillTheFileSystem && verbose >= VERBOSE_0) {
        fprintf(stdout, "Run ended: %s", CurrentTimeString());
    }
#endif /* USE_UNDOC_OPT - fillTheFileSystem */

    SummarizeResults(test);

    MPI_CHECK(MPI_Comm_free(&testComm), "MPI_Comm_free() error");
    for (i = 0; i < 2; i++) 
    {
    free(test->writeVal[i]);
    free(test->readVal[i]);
    }
    free(test->aggFileSizeFromCalc);
    free(test->aggFileSizeFromStat);
    free(test->aggFileSizeFromXfer);
    free(test->aggFileSizeForBW);
    for (i = 0; i < 12; i++) {
        free(timer[i]);
    }
    /* Sync with the tasks that did not participate in this test */
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD), "barrier error");

} /* TestIoSys() */


/******************************************************************************/
/*
 * Determine any spread (range) between node times.
 */

double
TimeDeviation(void)
{
    double timestamp, min = 0, max = 0, roottimestamp;

    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD), "barrier error");
    timestamp = GetTimeStamp();
    MPI_CHECK(MPI_Reduce(&timestamp, &min, 1, MPI_DOUBLE,
                         MPI_MIN, 0, MPI_COMM_WORLD),
              "cannot reduce tasks' times");
    MPI_CHECK(MPI_Reduce(&timestamp, &max, 1, MPI_DOUBLE,
                         MPI_MAX, 0, MPI_COMM_WORLD),
              "cannot reduce tasks' times");

    /* delta between individual nodes' time and root node's time */
    roottimestamp = timestamp;
    MPI_CHECK(MPI_Bcast(&roottimestamp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD),
              "cannot broadcast root's time");
    wall_clock_delta = timestamp - roottimestamp;

    return max-min;
} /* TimeDeviation() */


/******************************************************************************/
/*
 * Determine if valid tests from parameters.
 */

void
ValidTests(IOR_param_t * test)
{
    /* get the version of the tests */

    AioriBind(test->api);
    IOR_SetVersion(test);

    if (test->repetitions <= 0)
        WARN_RESET("too few test repetitions", repetitions);
    if (test->numTasks <= 0)
        ERR("too few tasks for testing");
    if (test->interTestDelay < 0)
        WARN_RESET("inter-test delay must be nonnegative value",
                   interTestDelay);
    if (test->readFile != TRUE && test->writeFile != TRUE
        && test->checkRead != TRUE && test->checkWrite != TRUE)
        ERR("test must write, read, or check file");
    if ((test->deadlineForStonewalling > 0)
        && (test->checkWrite == TRUE || test->checkRead == TRUE))
        ERR("can not perform write or read check with stonewalling");
    if (test->segmentCount < 0)
        ERR("segment count must be positive value");
    if ((test->blockSize % sizeof(IOR_size_t)) != 0)
        ERR("block size must be a multiple of access size");
    if (test->blockSize < 0)
        ERR("block size must be non-negative integer");
    if ((test->transferSize % sizeof(IOR_size_t)) != 0)
        ERR("transfer size must be a multiple of access size");
    if (test->setAlignment < 0)
        ERR("alignment must be non-negative integer");
    if (test->transferSize < 0)
        ERR("transfer size must be non-negative integer");
    if (test->transferSize == 0) {
        ERR("test will not complete with zero transfer size");
    } else {
        if ((test->blockSize % test->transferSize) != 0)
            ERR("block size must be a multiple of transfer size");
    }
    if (test->blockSize < test->transferSize)
        ERR("block size must not be smaller than transfer size");
    if ((strcmp(test->api, "MPIIO") == 0)
        && (test->blockSize < sizeof(IOR_size_t)
            || test->transferSize < sizeof(IOR_size_t)))
        ERR("block/transfer size may not be smaller than IOR_size_t for MPIIO");
    if ((strcmp(test->api, "HDF5") == 0)
        && (test->blockSize < sizeof(IOR_size_t)
            || test->transferSize < sizeof(IOR_size_t)))
        ERR("block/transfer size may not be smaller than IOR_size_t for HDF5");
    if ((strcmp(test->api, "NCMPI") == 0)
        && (test->blockSize < sizeof(IOR_size_t)
            || test->transferSize < sizeof(IOR_size_t)))
        ERR("block/transfer size may not be smaller than IOR_size_t for NCMPI");
    if ((strcmp(test->api, "NCMPI") == 0)
        && ((test->numTasks * test->blockSize * test->segmentCount)
            > (2*(IOR_offset_t)GIBIBYTE)))
        ERR("file size must be < 2GiB");
    if((test->useFileView == TRUE)
        && (sizeof(MPI_Aint) < 8) /* used for 64-bit datatypes */
        && ((test->numTasks * test->blockSize) > (2*(IOR_offset_t)GIBIBYTE)))
        ERR("segment size must be < 2GiB");
    if ((strcmp(test->api, "POSIX") != 0) && test->singleXferAttempt)
        WARN_RESET("retry only available in POSIX", singleXferAttempt);
    if ((strcmp(test->api, "POSIX") != 0) && test->fsync)
        WARN_RESET("fsync() only available in POSIX", fsync);
    if ((strcmp(test->api, "MPIIO") != 0) && test->preallocate)
        WARN_RESET("preallocation only available in MPIIO", preallocate);
    if ((strcmp(test->api, "MPIIO") != 0) && test->useFileView)
        WARN_RESET("file view only available in MPIIO", useFileView);
    if ((strcmp(test->api, "MPIIO") != 0) && test->useSharedFilePointer)
        WARN_RESET("shared file pointer only available in MPIIO",
                   useSharedFilePointer);
    if ((strcmp(test->api, "MPIIO") == 0) && test->useSharedFilePointer)
        WARN_RESET("shared file pointer not implemented", useSharedFilePointer);
    if ((strcmp(test->api, "MPIIO") != 0) && test->useStridedDatatype)
        WARN_RESET("strided datatype only available in MPIIO",
                   useStridedDatatype);
    if ((strcmp(test->api, "MPIIO") == 0) && test->useStridedDatatype)
        WARN_RESET("strided datatype not implemented", useStridedDatatype);
    if ((strcmp(test->api, "MPIIO") == 0)
        && test->useStridedDatatype
        && (test->blockSize < sizeof(IOR_size_t)
            || test->transferSize < sizeof(IOR_size_t)))
        ERR("need larger file size for strided datatype in MPIIO");
    if ((strcmp(test->api, "POSIX") == 0) && test->showHints)
        WARN_RESET("hints not available in POSIX", showHints);
    if ((strcmp(test->api, "POSIX") == 0) && test->collective)
        WARN_RESET("collective not available in POSIX", collective);
    if (test->reorderTasks == TRUE && test->reorderTasksRandom == TRUE)
        ERR("Both Constant and Random task re-ordering specified. Choose one and resubmit");
    if (test->randomOffset && test->reorderTasksRandom && test->filePerProc == FALSE)
        ERR("random offset and random reorder tasks specified with single-shared-file. Choose one and resubmit");
    if (test->randomOffset && test->reorderTasks && test->filePerProc == FALSE)
        ERR("random offset and constant reorder tasks specified with single-shared-file. Choose one and resubmit");
    if (test->randomOffset && test->checkRead)
        ERR("random offset not available with read check option (use write check)");
    if (test->randomOffset && test->storeFileOffset)
        ERR("random offset not available with store file offset option)");
    if ((strcmp(test->api, "MPIIO") == 0) && test->randomOffset && test->collective)
        ERR("random offset not available with collective MPIIO");
    if ((strcmp(test->api, "MPIIO") == 0) && test->randomOffset && test->useFileView)
        ERR("random offset not available with MPIIO fileviews");
    if ((strcmp(test->api, "HDF5") == 0) && test->randomOffset)
        ERR("random offset not available with HDF5");
    if ((strcmp(test->api, "NCMPI") == 0) && test->randomOffset)
        ERR("random offset not available with NCMPI");
    if ((strcmp(test->api, "HDF5") != 0) && test->individualDataSets)
        WARN_RESET("individual datasets only available in HDF5",
                   individualDataSets);
    if ((strcmp(test->api, "HDF5") == 0) && test->individualDataSets)
        WARN_RESET("individual data sets not implemented", individualDataSets);
    if ((strcmp(test->api, "NCMPI") == 0) && test->filePerProc)
        ERR("file-per-proc not available in current NCMPI");
    if (test->noFill) {
        if (strcmp(test->api, "HDF5") != 0) {
            ERR("'no fill' option only available in HDF5");
        } else {
            /* check if hdf5 available */
            #if defined (H5_VERS_MAJOR) && defined (H5_VERS_MINOR)
                /* no-fill option not available until hdf5-1.6.x */
                #if (H5_VERS_MAJOR > 0 && H5_VERS_MINOR > 5)
                    ;
                #else
                    char errorString[MAX_STR];
                    sprintf(errorString, "'no fill' option not available in %s",
                            test->apiVersion);
                    ERR(errorString);
                #endif
            #else
                WARN("unable to determine HDF5 version for 'no fill' usage");
            #endif
        }
    }
    if (test->useExistingTestFile
        && (test->lustre_stripe_count != 0
            || test->lustre_stripe_size != 0
            || test->lustre_start_ost != -1))
        ERR("Lustre stripe options are incompatible with useExistingTestFile");
} /* ValidTests() */

IOR_offset_t *
GetOffsetArraySequential(IOR_param_t *test, int pretendRank)
{
    IOR_offset_t i, j, k = 0;
    IOR_offset_t offsets;
    IOR_offset_t *offsetArray;

    /* count needed offsets */
    offsets = (test->blockSize / test->transferSize) * test->segmentCount;

    /* setup empty array */
    offsetArray = (IOR_offset_t *)malloc((offsets+1) * sizeof(IOR_offset_t));
    if (offsetArray == NULL) ERR("out of memory");
    offsetArray[offsets] = -1; /* set last offset with -1 */

    /* fill with offsets */
    for (i = 0; i < test->segmentCount; i++) {
        for (j = 0; j < (test->blockSize / test->transferSize); j++) {
            offsetArray[k] = j * test->transferSize;
            if (test->filePerProc) {
                offsetArray[k] += i * test->blockSize;
            } else {
                offsetArray[k] += (i * test->numTasks * test->blockSize)
                                  + (pretendRank * test->blockSize);
            }
            k++;
        }
    }

    return(offsetArray);
} /* GetOffsetArraySequential() */


IOR_offset_t *
GetOffsetArrayRandom(IOR_param_t *test, int pretendRank, int access)
{
    int seed;
    IOR_offset_t i, value, tmp;
    IOR_offset_t offsets = 0, offsetCnt = 0;
    IOR_offset_t fileSize;
    IOR_offset_t *offsetArray;

    /* set up seed for random() */
    if (access == WRITE || access == READ) {
        test->randomSeed = seed = random();
    } else {
        seed = test->randomSeed;
    }
    srandom(seed);

    fileSize = test->blockSize * test->segmentCount;
    if (test->filePerProc == FALSE) {
        fileSize *= test->numTasks;
    }

    /* count needed offsets (pass 1) */
    for (i = 0; i < fileSize; i += test->transferSize) {
        if (test->filePerProc == FALSE) {
            if ((random() % test->numTasks) == pretendRank) {
                offsets++;
            }
        } else {
            offsets++;
        }
    }

    /* setup empty array */
    offsetArray = (IOR_offset_t *)malloc((offsets+1) * sizeof(IOR_offset_t));
    if (offsetArray == NULL) ERR("out of memory");
    offsetArray[offsets] = -1; /* set last offset with -1 */

    if (test->filePerProc) {
        /* fill array */
        for (i = 0; i < offsets; i++) {
            offsetArray[i] = i * test->transferSize;
        }
    } else {
        /* fill with offsets (pass 2) */
        srandom(seed); /* need same seed */
        for (i = 0; i < fileSize; i += test->transferSize) {
            if ((random() % test->numTasks) == pretendRank) {
                offsetArray[offsetCnt] = i;
                offsetCnt++;
            }
        }
    }
    /* reorder array */
    for (i = 0; i < offsets; i++) {
        value = random() % offsets;
        tmp = offsetArray[value];
        offsetArray[value] = offsetArray[i];
        offsetArray[i] = tmp;
    }
    SeedRandGen(test->testComm); /* synchronize seeds across tasks */

    return(offsetArray);
} /* GetOffsetArrayRandom() */


/******************************************************************************/
/*
 * Write or Read data to file(s).  This loops through the strides, writing
 * out the data to each block in transfer sizes, until the remainder left is 0.
 */

IOR_offset_t
WriteOrRead(IOR_param_t * test,
            void        * fd,
            int           access)
{
    int            errors = 0;
    IOR_offset_t   amtXferred,
                   transfer,
                   transferCount = 0,
                   pairCnt = 0,
                 * offsetArray;
    int            pretendRank;
    void         * buffer = NULL;
    void         * checkBuffer = NULL;
    void         * readCheckBuffer = NULL;
    IOR_offset_t   dataMoved = 0;             /* for data rate calculation */
    double         startForStonewall;
    int            hitStonewall;


    /* initialize values */
    pretendRank = (rank + rankOffset) % test->numTasks;

    if (test->randomOffset) {
        offsetArray = GetOffsetArrayRandom(test, pretendRank, access);
    } else {
        offsetArray = GetOffsetArraySequential(test, pretendRank);
    }

    SetupXferBuffers(&buffer, &checkBuffer, &readCheckBuffer,
                     test, pretendRank, access);

    /* check for stonewall */
    startForStonewall = GetTimeStamp();
    hitStonewall = ((test->deadlineForStonewalling != 0)
                    && ((GetTimeStamp() - startForStonewall)
                        > test->deadlineForStonewalling));

    /* loop over offsets to access */
    while ((offsetArray[pairCnt] != -1) && !hitStonewall) {
        test->offset = offsetArray[pairCnt];
        /*
         * fills each transfer with a unique pattern
         * containing the offset into the file
         */
        if (test->storeFileOffset == TRUE) {
            FillBuffer(buffer, test, test->offset, pretendRank);
        }
        transfer = test->transferSize;
        if (access == WRITE) {
            amtXferred = IOR_Xfer(access, fd, buffer, transfer, test);
            if (amtXferred != transfer) ERR("cannot write to file");
        } else if (access == READ) {
            amtXferred = IOR_Xfer(access, fd, buffer, transfer, test);
            if (amtXferred != transfer) ERR("cannot read from file");
        } else if (access == WRITECHECK) {
            memset(checkBuffer, 'a', transfer);
            amtXferred = IOR_Xfer(access, fd, checkBuffer, transfer, test);
            if (amtXferred != transfer)
                ERR("cannot read from file write check");
            transferCount++;
            errors += CompareBuffers(buffer, checkBuffer, transfer,
                                     transferCount, test, WRITECHECK);
        } else if (access == READCHECK){
            ReadCheck(fd, buffer, checkBuffer, readCheckBuffer, test,
                      transfer, test->blockSize, &amtXferred,
                      &transferCount, access, &errors);
        }
        dataMoved += amtXferred;
        pairCnt++;

        hitStonewall = ((test->deadlineForStonewalling != 0)
                        && ((GetTimeStamp() - startForStonewall)
                            > test->deadlineForStonewalling));
    }

    totalErrorCount += CountErrors(test, access, errors);

    FreeBuffers(access, checkBuffer, readCheckBuffer, buffer, offsetArray);

    if (access == WRITE && test->fsync == TRUE) {
        IOR_Fsync(fd, test); /*fsync after all accesses*/
    }
    return(dataMoved);
} /* WriteOrRead() */

/****************************ec functions*******************************************/


/*****************experiment timers***************/
double *ec_timers;
double ec_startTime;
double ec_endTime;
double ec_strategy_startTime;
double ec_strategy_endTime;
/*****************experiment timers***************/

/*****************thread variables****************/
pthread_t *threads;
int *tranferDone;
int *dataLeft;
volatile int numTransferred;
pthread_mutex_t lockOfNT = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lockOfDecodeSum = PTHREAD_MUTEX_INITIALIZER;
int decode_sum = 0;

/*********for collecitve thread*********/
#define STRATEGY_STAGE 1
#define NORMAL_STAGE 0
volatile int current_stage;
pthread_mutex_t lock_hasStraggler = PTHREAD_MUTEX_INITIALIZER;
volatile int hasStraggler;
volatile int currentStragger;
pthread_mutex_t lock_firstStraggler = PTHREAD_MUTEX_INITIALIZER;
volatile int firstStraggler = 1;
pthread_mutex_t lock_rcListHead = PTHREAD_MUTEX_INITIALIZER;
LIST_HEAD(rcListHead);
pthread_mutex_t lock_strategyIsReady  =PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t strategyReady = PTHREAD_COND_INITIALIZER;
pthread_cond_t cond_hasStraggler = PTHREAD_COND_INITIALIZER;
pthread_cond_t active_parity = PTHREAD_COND_INITIALIZER;
int strategyIsReady;
IOR_offset_t *currentPosOfThread;
int* strategyReceived;
IOR_offset_t next_pairCnt;
int leftThreads;
int start_flag;

    volatile int start_strategy;
pthread_mutex_t lock_start_strategy = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond_start_strategy = PTHREAD_COND_INITIALIZER;
pthread_mutex_t lock_response_number = PTHREAD_MUTEX_INITIALIZER;
int response_number;
int *parity_hangup;
int RC_flag;
int RC_id = -1;
IOR_offset_t RC_count;

char **sample_data = NULL;
char **sample_coding = NULL;
int *sample_matrix = NULL;
int *sample_bitmatrix = NULL;
int **sample_schedule = NULL;
int sample_erasures[] = {0,-1};

IOR_offset_t min_offset_of_stripes(int n, int m){
    //min offset of stripes from n to m
    int i;
    IOR_offset_t res = currentPosOfThread[n];
    for(i=n+1; i<=m;i++){
        res = MIN(res,currentPosOfThread[i]);
    }
    return res;
}

IOR_offset_t max_offset_of_stripes(int n, int m){
    int i;
    IOR_offset_t res = currentPosOfThread[n];
    for(i=n+1; i<=m;i++){
        res = MAX(res,currentPosOfThread[n]);
    }
    return res;
}

IOR_offset_t kth_large_offset_of_stripes(int n, int k){
    //n个stripe中，第k大的数
    IOR_offset_t *A = (IOR_offset_t*)malloc(sizeof(IOR_offset_t)*n);
    memcpy(A, currentPosOfThread,sizeof(IOR_offset_t)*n);
    
    int i,j,temp;
    for(i=0;i<(n-1);i++){
        for(j=i;j<n;j++){
            if(A[i]<A[j]){
                temp = A[i];
                A[i] = A[j];
                A[j] = temp;
            }
        }
        if(i==(k-1)){
            return A[i];
        }
    }
}
/*********for collecitve thread*********/

int hitStonewall;
char ***omp_data;
char ***omp_coding;
int omp_thread_num = 8;
ec_decode_thread_args ec_decode_arg;
int decode_res = 0;
int finished_thread = 0;


typedef struct rcBlock{
    struct list_head list; //linux list head structure
    IOR_offset_t start;
    IOR_offset_t end;
}rcBlock;

/*****************thread parameters****************/
void *
ec_read_thread(ec_read_thread_args *arg)
{
    int id = arg->id;
    int k = arg->test->ec_k;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;
    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    int pairCnt = 0;
    double startTime = 0;
    double endTime = 0;
    double average_xfer_time = 0;
    startTime = GetTimeStamp();

    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {
        offset = offsetArray[pairCnt];
        offset = offset / arg->test->ec_k;
        //XferStartTime = GetTimeStamp() - startTime;
        if (id < k)
        {
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_data)[id], arg->test->ec_stripe_size, arg->test, offset);
        }
        else
        {
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_coding)[id - k], arg->test->ec_stripe_size, arg->test, offset);
        }
        //XferEndTime = GetTimeStamp() - startTime;
        //fprintf(stdout, "#Xferid=%d,startTime=%0.2lf,endTIme=%0.2lf,duration=%2lf\n",id, XferStartTime,XferEndTime,XferEndTime-XferStartTime);
        pairCnt++;
        dataLeft[id]--;
    }
    //fprintf(stdout, "reading file%d complete, size: %lld\n", id, transferred_size);

    pthread_mutex_lock(&lockOfNT);
    tranferDone[id] = 1;
    numTransferred += 1;
    pthread_mutex_unlock(&lockOfNT);
    endTime = GetTimeStamp();
    average_xfer_time = (endTime - startTime) / pairCnt;
    fprintf(stdout, "average_xfer_time = %lf\n", average_xfer_time);
    ec_timers[id] += endTime - startTime;
}

double parity_time[2];
int parity_target[2] = {0,1};
int parity_start[2];
int parity_number[2];
int slow_start;
int slow_num;
double slow_time;
int slow_target;
double decode_time0;
double decode_time1;
double decode_num;

int batch_size = 32;
int small_batch = 16;
char *batch_buffer0;
char *batch_buffer1;

pthread_mutex_t buffernum0 = PTHREAD_MUTEX_INITIALIZER;
int bnum0;
pthread_mutex_t buffernum1 = PTHREAD_MUTEX_INITIALIZER;
int bnum1;

void*
ec_slowread_thread(ec_read_thread_args *arg){
    double startTime = GetTimeStamp();
    int i;
    for(i=0;i<slow_num;i++){
        IOR_Xfer_ec(arg->access, (arg->fds)[slow_target], (arg->ec_data)[slow_target], arg->test->ec_stripe_size, arg->test, arg->offSetArray[slow_start+i]);
    }
    double endTime = GetTimeStamp();
    slow_time = endTime - startTime;
}

void*
ec_parity_thread0(ec_read_thread_args *arg){
    double startTime = GetTimeStamp();
    int i;
    int nouse;
    int it = parity_number[0]/small_batch;
    int it2 = parity_number[0]%small_batch;
    // fprintf(stdout,"it1 = %d, it2 = %d\n", it,it2);
    for(i=0;i<it;i++){

        IOR_Xfer_ec(arg->access, (arg->fds)[6+parity_target[0]], batch_buffer0, small_batch*arg->test->ec_stripe_size, arg->test, arg->offSetArray[parity_start[0]+i*small_batch]);
        //nouse = 1400000; //parity = 4
        nouse = 500000; //parity <= 3
        while(nouse--){
            ;
        }
        pthread_mutex_lock(&buffernum0);
        bnum0+=small_batch;
        // fprintf(stdout,"parity 0 add small batch\n");
        pthread_mutex_unlock(&buffernum0);
    }

    for(i=0;i<it2;i++){

        IOR_Xfer_ec(arg->access, (arg->fds)[6+parity_target[0]], batch_buffer0, arg->test->ec_stripe_size, arg->test, arg->offSetArray[parity_start[0]+it*small_batch+i]);
        //nouse = 1400000; //parity = 4
        nouse = 500000; //parity <= 3
        while(nouse--){
            ;
        }
        pthread_mutex_lock(&buffernum0);
        bnum0+=1;
        // fprintf(stdout,"parity 0 add 1\n");
        pthread_mutex_unlock(&buffernum0);
    }

    double endTime = GetTimeStamp();
    parity_time[parity_target[0]] = endTime - startTime;
}

void *
ec_parity_thread1(ec_read_thread_args *arg)
{
    double startTime = GetTimeStamp();
    int i;
    int nouse;
    // fprintf(stdout, "parity_number[1] = %d\n", parity_number[1]);
    int it = parity_number[1] / small_batch;
    int it2 = parity_number[1] % small_batch;
    // fprintf(stdout, "it1 = %d, it2 = %d\n", it, it2);
    for (i = 0; i < it; i++)
    {

        IOR_Xfer_ec(arg->access, (arg->fds)[6 + parity_target[1]], batch_buffer1, small_batch * arg->test->ec_stripe_size, arg->test, arg->offSetArray[parity_start[1] + i * small_batch]);
        //nouse = 1400000; //parity = 4
        nouse = 500000; //parity <= 3
        while (nouse--)
        {
            ;
        }                   
        pthread_mutex_lock(&buffernum1);
        bnum1 += small_batch;
        // fprintf(stdout,"parity 1 add small batch\n");
        pthread_mutex_unlock(&buffernum1);
    }

    for (i = 0; i < it2; i++)
    {

        IOR_Xfer_ec(arg->access, (arg->fds)[6 + parity_target[1]], batch_buffer1, arg->test->ec_stripe_size, arg->test, arg->offSetArray[parity_start[1]+it*small_batch + i]);
        //nouse = 1400000; //parity = 4
        nouse = 500000; //parity <= 3
        while (nouse--)
        {
            ;
        }
        pthread_mutex_lock(&buffernum1);
        // fprintf(stdout,"parity 1 add 1\n");
        bnum1 += 1;
        pthread_mutex_unlock(&buffernum1);
    }
    double endTime = GetTimeStamp();
    parity_time[parity_target[1]] = endTime - startTime;
}


void *
ec_adaptive_decode0()
{
    double startTime = GetTimeStamp();
    int i;
    for (i = 0; i < parity_number[0]; i++)
    {
        
        while (1)
        {
            pthread_mutex_lock(&buffernum0);
            if(bnum0>0){
                bnum0--;
                // fprintf(stdout, "decode0 consume 1\n");
                pthread_mutex_unlock(&buffernum0);
                break;
            }
            pthread_mutex_unlock(&buffernum0);            
        }
        
        jerasure_schedule_decode_lazy(6, 2, 8, sample_matrix, sample_erasures, sample_data, sample_coding, 524288, 8, 1);
        //sleep(0.1);
    }
    double endTime = GetTimeStamp();
    decode_time0 = endTime - startTime;
}

void *
ec_adaptive_decode1()
{
    double startTime = GetTimeStamp();
    int i;
    for (i = 0; i < parity_number[1]; i++)
    {       
        while (1)
        {
            pthread_mutex_lock(&buffernum1);
            if(bnum1>0){
                bnum1--;
                // fprintf(stdout, "decode1 consume 1\n");
                pthread_mutex_unlock(&buffernum1);
                break;
            }
            pthread_mutex_unlock(&buffernum1);            
        }
        
        jerasure_schedule_decode_lazy(6, 2, 8, sample_matrix, sample_erasures, sample_data, sample_coding, 524288, 8, 1);
        //sleep(0.1);
    }
    double endTime = GetTimeStamp();
    decode_time1 = endTime - startTime;
}

void *
ec_request_test(ec_read_thread_args *arg)
{
    int id = arg->id;
    int k = arg->test->ec_k;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;
    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    int pairCnt = 0;
    double startTime = 0;
    double endTime = 0;
    double average_xfer_time = 0;
    char* buffer = (char*)malloc(sizeof(char)*16*arg->test->transferSize);
    startTime = GetTimeStamp();

    /*adaptive parameter*/
    int isStraggler = 0;
    double times_over_threshold = 0;
    double times_below_threshold = 0;
    double xfer_startTime;
    double xfer_endTime;
    double tempTime;
    double upper_threshold = 0.02;
    double lower_threshold = 0.008;
    double duration = 0.00;
    int C = 3;
    int S = 1;
    int num_0 = 3;
    int num_1 = 1;
    int should_decode = 0;
    int should_read = 0;
    int should_readfrom0 = 0;
    int should_readfrom1 = 0;
    IOR_offset_t window_size = 8;
    IOR_offset_t window_threshold = 1024;
    IOR_offset_t total_restruction = 0;

    int i;

    /********************request test*******************/

        double encode_test_start = GetTimeStamp();
        for(i=0;i<100;i++){
            jerasure_schedule_decode_lazy(6, 2, 8, sample_matrix, sample_erasures, sample_data, sample_coding, arg->test->ec_stripe_size, 8, 1);
        }       
        double encode_test_end = GetTimeStamp();
        fprintf(stdout, "##encode latency of %ld: %lf\n", arg->test->ec_stripe_size,(encode_test_end- encode_test_start)/100);
    
        // offset = offsetArray[0];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: single 512K: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[1];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, arg->test->ec_stripe_size, arg->test, offset);
        // offset = offsetArray[2];
        // tempTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: two 512K: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[3];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 2 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: single 1MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[5];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 2 * arg->test->ec_stripe_size, arg->test, offset);
        // offset = offsetArray[7];
        // tempTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 2 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: two 1MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[9];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 4 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: single 2MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[13];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 4 * arg->test->ec_stripe_size, arg->test, offset);
        // offset = offsetArray[17];
        // tempTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 4 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: two 2MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[21];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 8 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: single 4MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[29];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 8 * arg->test->ec_stripe_size, arg->test, offset);
        // offset = offsetArray[37];
        // tempTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 8 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: two 4MB: %lf\n", xfer_endTime - xfer_startTime);
        // offset = offsetArray[45];
        // xfer_startTime = GetTimeStamp() - startTime;
        // transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], buffer, 16 * arg->test->ec_stripe_size, arg->test, offset);
        // xfer_endTime = GetTimeStamp() - startTime;
        // fprintf(stdout, "##request test: single 8MB: %lf\n", xfer_endTime - xfer_startTime);

        // double encode_test_start = GetTimeStamp();
        // jerasure_schedule_decode_lazy(6, 2, 8, sample_matrix, sample_erasures, sample_data, sample_coding, 524288, 8, 1);
        // double encode_test_end = GetTimeStamp();
        // fprintf(stdout, "##encode latency: %lf\n", encode_test_end- encode_test_start);
}

void *
ec_adaptive_thread(ec_read_thread_args *arg)
{
    int id = arg->id;
    int k = arg->test->ec_k;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;
    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    int pairCnt = 0;
    double startTime = 0;
    double endTime = 0;
    double average_xfer_time = 0;
    startTime = GetTimeStamp();

    /*adaptive parameter*/
    batch_buffer0 = (char*)malloc(sizeof(char)*batch_size*524288);
    batch_buffer1 = (char*)malloc(sizeof(char)*batch_size*524288);

    int isStraggler = 0;
    double times_over_threshold = 0;
    double times_below_threshold = 0;
    double xfer_startTime;
    double xfer_endTime;
    double upper_threshold = 0.008;
    double lower_threshold = 0.005;
    double duration = 0.00;
    double C_latency;
    double S_latency;
    double P_latency;
    
    int should_decode = 0;
    int should_read = 0;
    int should_readfrom0 =0;
    int should_readfrom1 = 0;
    IOR_offset_t window_size = 8;
    IOR_offset_t window_threshold = 1024;
    IOR_offset_t total_restruction = 0;

    /********************request test*******************/
    
    /********************request test*******************/
    offset = 0;
    /*adaptive parameter*/
    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {
        /****************straggler_detector******************/
        if (duration > upper_threshold && !isStraggler)
        { //不是straggler，等待进入
            times_over_threshold++;
        }
        else if (duration < lower_threshold && isStraggler)
        { //是straggler，等待退出
            times_below_threshold++;
        }
        else if (duration <= lower_threshold && !isStraggler)
        { //不是straggler，进入状态清零
            times_over_threshold = 0;
        }
        else if (duration >= upper_threshold && isStraggler)
        { //是straggler，退出状态清零
            times_below_threshold = 0;
        }

        if (times_over_threshold >= 8)
        {
            isStraggler = 1;
            RC_id = id;
            times_over_threshold = 0;
            fprintf(stdout, "process %d: thread %d into straggler, offset: %lld\n", rank, id, pairCnt);
        }
        if (times_below_threshold >= 8)
        {
            isStraggler = 0;
            times_below_threshold = 0;
            if (currentStragger == id)
            {
                pthread_mutex_lock(&lock_hasStraggler);
                hasStraggler = 0;
                currentStragger = -1;
                pthread_mutex_unlock(&lock_hasStraggler);
            }
            fprintf(stdout, "process %d:thread %d quit straggler, offset: %lld\n", rank, id, pairCnt);
        }
        /****************straggler_detector******************/
        int temp_pairCnt;
        
        pthread_t parity_threads[2];
        pthread_t slow_read;
        pthread_t decode_thread0;
        pthread_t decode_thread1;
        /****************window_size******************/
        while(isStraggler){
            //window_size = 2*window_size;
            //window_size = 16;
            // if(window_size>=window_threshold){
            //     window_size/=2;
            // }
            // if(window_size >= (num_reconstruct-pairCnt)){
            //     window_size = num_reconstruct-pairCnt;
            // }
            // if(window_size >= (num_reconstruct-pairCnt)){
            //     window_size = num_reconstruct-pairCnt;
            // }
            // 
            // should_read = batch_size;
            double C = 8;
            double S = 1;
            double num_0 = 1;
            double num_1 = 1;
            should_read = 0;
            should_decode = batch_size/S*C;
            should_readfrom0 = should_decode/(num_0+num_1)*num_0;
            should_readfrom1 = should_decode - should_readfrom0;
            // should_decode = 1;
            // should_read = 1;
            // should_readfrom0 = 1;
            // should_readfrom1 = 0;
            window_size = should_decode+should_read;
            if(window_size >= (num_reconstruct-pairCnt)){
                //break;
                isStraggler = 0;
                break;
            }
            fprintf(stdout, "windowsize: %d, decode: %d, read: %d, read from 0: %d, read from 1: %d\n",window_size, should_decode,should_read,should_readfrom0,should_readfrom1);
            temp_pairCnt = pairCnt;
            parity_start[0] = temp_pairCnt;
            parity_start[1] = temp_pairCnt+should_readfrom0;
            parity_number[0] = should_readfrom0;
            parity_number[1] = should_readfrom1;
            slow_target = id;
            slow_num = should_read;
            slow_start = temp_pairCnt+should_readfrom0+should_readfrom1;
            decode_num = should_decode;
            pthread_create(&parity_threads[0], NULL, ec_parity_thread0, arg);
            pthread_create(&parity_threads[1], NULL, ec_parity_thread1, arg);
            pthread_create(&slow_read, NULL, ec_slowread_thread, arg);
            //pthread_join(parity_threads[0], NULL);
            pthread_create(&decode_thread0, NULL ,ec_adaptive_decode0, NULL);
            pthread_create(&decode_thread1, NULL ,ec_adaptive_decode1, NULL);
            //pthread_join(parity_threads[0], NULL);
            pthread_join(parity_threads[1], NULL);
            pthread_join(parity_threads[0], NULL);
            pthread_join(slow_read, NULL);
            pthread_join(decode_thread0, NULL);
            pthread_join(decode_thread1, NULL);
            fprintf(stdout, "parity time0: %lf, parity time1: %lf\n", parity_time[0],parity_time[1]);
            fprintf(stdout, "slow read time: %lf\n", slow_time);
            fprintf(stdout, "decode time0: %lf\n", decode_time0);
            pairCnt += window_size;
            isStraggler = 1;
            times_over_threshold = 9;
            total_restruction += window_size;
            fprintf(stdout, "exit straggler in while\n");
            break;
        }

        /****************window_size******************/
        if(!isStraggler){
            offset = offsetArray[pairCnt];
            offset = offset / arg->test->ec_k;
            //XferStartTime = GetTimeStamp() - startTime;
            
            xfer_startTime = GetTimeStamp() -startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_data)[id], arg->test->ec_stripe_size, arg->test, offset);
        
            
            xfer_endTime = GetTimeStamp() - startTime;
            duration = xfer_endTime - xfer_startTime;
            //fprintf(stdout, "#Xferid=%d,startTime=%0.2lf,endTIme=%0.2lf,duration=%2lf\n",id, XferStartTime,XferEndTime,XferEndTime-XferStartTime);
            pairCnt++;
            dataLeft[id]--;
        }

        // if(pairCnt>=num_reconstruct)
        //     pairCnt = num_reconstruct;
    }
    //fprintf(stdout, "reading file%d complete, size: %lld\n", id, transferred_size);

    pthread_mutex_lock(&lockOfNT);
    tranferDone[id] = 1;
    numTransferred += 1;
    pthread_mutex_unlock(&lockOfNT);
    endTime = GetTimeStamp();
    average_xfer_time = (endTime - startTime) / pairCnt;
    fprintf(stdout, "#thread %d: total reconstruction: %ld\n", id, total_restruction);
    fprintf(stdout, "#thread %d: total paircnt: %ld\n", id, num_reconstruct);
    fprintf(stdout, "#thread %d: average_xfer_time = %lf\n", id, average_xfer_time);
    ec_timers[id] += endTime - startTime;
}

void *
ec_collective_thread2(ec_read_thread_args *arg)
{
    int id = arg->id;
    int k = arg->test->ec_k;
    int m = arg->test->ec_m;
    int w = arg->test->ec_w;
    enum Coding_Technique method = arg->test->ec_method;
    IOR_offset_t ec_blocksize = arg->test->ec_stripe_size;
    int ec_packetsize = arg->test->ec_packetsize;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;

    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    IOR_offset_t pairCnt = 0;
    double startTime = 0;
    double endTime = 0;

    int isOriginReader = id < k ? 1 : 0;
    int isStraggler = 0;
    int transferTime = 0;
    double xfer_startTime = 0;
    double xfer_endTime = 0;
    double duration = 0;
    double times_over_threshold = 0;
    double times_below_threshold = 0;
    //double threshold = 0.02;
    double upper_threshold = 0.02;
    double lower_threshold = 0.01;
    startTime = GetTimeStamp();

    fprintf(stdout,"break at 3524\n");
    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {
        offset = offsetArray[pairCnt];
        offset = offset / arg->test->ec_k;

        /*****************collective oeration*********************/
        /****************is_straggler******************/
        if (duration > upper_threshold && !isStraggler)
        { //不是straggler，等待进入
            times_over_threshold++;
        }
        else if (duration < lower_threshold && isStraggler)
        { //是straggler，等待退出
            times_below_threshold++;
        }
        else if (duration <= lower_threshold && !isStraggler)
        { //不是straggler，进入状态清零
            times_over_threshold = 0;
        }
        else if (duration >= upper_threshold && isStraggler)
        { //是straggler，退出状态清零
            times_below_threshold = 0;
        }

        if (times_over_threshold >= 10)
        {
            isStraggler = 1;
            times_over_threshold = 0;
            fprintf(stdout, "thread %d into straggler, offset: %lld\n", id, pairCnt);
        }
        if (times_below_threshold >= 5)
        {
            isStraggler = 0;
            times_below_threshold = 0;
            if (currentStragger == id)
            {
                pthread_mutex_lock(&lock_hasStraggler);
                hasStraggler = 0;
                currentStragger = -1;
                pthread_mutex_unlock(&lock_hasStraggler);
            }
            fprintf(stdout, "thread %d quit straggler, offset: %lld\n", id, pairCnt);
        }
        /****************is_straggler******************/

        /*****************collective oeration*********************/
        pthread_mutex_lock(&lock_hasStraggler);
        if(isStraggler && !hasStraggler){
            fprintf(stdout, "thread %d meets stress and there is no straggler currently.\n", id);
            
            hasStraggler = 1;
            currentStragger = id;
            fprintf(stdout, "thread %d will be the stragger.\n", id);
            current_stage = STRATEGY_STAGE;
            fprintf(stdout, "enter into strategy stage.\n", id);
            pthread_cond_broadcast(&cond_hasStraggler);
            fprintf(stdout, "wake up threads waiting for straggler\n", id);
            pthread_mutex_unlock(&lock_hasStraggler);
            int z;
            for (z = 0; z < (k + m); z++)
            {
                fprintf(stdout, "thread %d's current position: %lld\n", z, currentPosOfThread[z]);
            }
        }else{
            pthread_mutex_unlock(&lock_hasStraggler);
        }
        
        if(current_stage == STRATEGY_STAGE){
            if(id == currentStragger){
                fprintf(stdout, "thread %d is strategy maker, making strategy...\n", id);
                IOR_offset_t min_stripe_offset = min_offset_of_stripes(0,k-1);
                fprintf(stdout, "min_stripe_offset = %lld\n", min_stripe_offset);
                IOR_offset_t kth_large_stripe_offset;
                if(min_stripe_offset >= max_offset_of_stripes(k,k+m-1)){
                    fprintf(stdout, "stripe reader is faster\n");
                    next_pairCnt = kth_large_offset_of_stripes(k+m,k);
                    fprintf(stdout, "all readers should read at least from %lld\n", next_pairCnt);
                    if(id>=k){
                        fprintf(stdout,"hang up parity reader since no compute & slower\n");
                        parity_hangup[id-k] = 1;
                    }
                    pairCnt = MAX(pairCnt, next_pairCnt); //set pairCnt to larger one
                    currentPosOfThread[id] = pairCnt;
                    //fprintf(stdout, "thread %d: try to make parity reader work.\n", id);
                    fprintf(stdout, "thread %d: strategy is ready, broadcasting...\n", id);
                    pthread_mutex_lock(&lock_strategyIsReady);
                    strategyIsReady = 1;
                    pthread_cond_broadcast(&strategyReady);
                    pthread_mutex_unlock(&lock_strategyIsReady);
                    fprintf(stdout, "thread %d: waiting for response...\n", id);
                    while (1)
                    {
                        pthread_mutex_lock(&lock_response_number);
                        if(response_number == (leftThreads-1)){
                            fprintf(stdout, "thread %d: broadcasting success!\n", id);
                            pthread_mutex_unlock(&lock_response_number);
                            break;
                        }
                        pthread_mutex_unlock(&lock_response_number);
                    }
                    pthread_mutex_lock(&start_strategy);
                    strategyIsReady = 0;
                    current_stage = NORMAL_STAGE;
                    response_number = 0;
                    start_flag = 1;
                    pthread_cond_broadcast(&cond_start_strategy);
                    pthread_mutex_unlock(&start_strategy);
                    continue;
                }else{
                    kth_large_stripe_offset = kth_large_offset_of_stripes(k+m,k);
                    RC_count = kth_large_stripe_offset - min_stripe_offset;
                    next_pairCnt = kth_large_stripe_offset;
                    fprintf(stdout, "thread %d: try to recompute from %lld to %lld.\n", id, min_stripe_offset,kth_large_stripe_offset);
                    RC_flag = 1;
                    RC_id = id;
                    fprintf(stdout, "all readers should read at least from %lld\n", next_pairCnt);
                    pairCnt = MAX(pairCnt, next_pairCnt);
                    currentPosOfThread[id] = pairCnt;
                    fprintf(stdout, "thread %d: strategy is ready, broadcasting...\n", id);
                    pthread_mutex_lock(&lock_strategyIsReady);
                    strategyIsReady = 1;
                    pthread_cond_broadcast(&strategyReady);
                    pthread_mutex_unlock(&lock_strategyIsReady);
                    fprintf(stdout, "thread %d: waiting for response...\n", id);
                    while (1)
                    {
                        pthread_mutex_lock(&lock_response_number);
                        if(response_number == (leftThreads-1)){
                            fprintf(stdout, "thread %d: broadcasting success!\n", id);
                            pthread_mutex_unlock(&lock_response_number);
                            break;
                        }
                        pthread_mutex_unlock(&lock_response_number);
                    }
                    pthread_mutex_lock(&start_strategy);
                    strategyIsReady = 0;
                    current_stage = NORMAL_STAGE;
                    response_number = 0;
                    start_flag = 1;
                    pthread_cond_broadcast(&cond_start_strategy);
                    pthread_mutex_unlock(&start_strategy);
                    continue;
                }

            }else{//strategy receiver
                fprintf(stdout, "thread %d:waiting for strategy...\n", id);
                pthread_mutex_lock(&lock_strategyIsReady);
                while (!strategyIsReady)
                {
                    pthread_cond_wait(&strategyReady, &lock_strategyIsReady);
                }
                fprintf(stdout, "thread %d:strategy received!\n", id);
                pairCnt = MAX(pairCnt, next_pairCnt);
                currentPosOfThread[id] = pairCnt;
                fprintf(stdout, "thread %d:current position: %lld\n", id, currentPosOfThread[id]);
                pthread_mutex_lock(&lock_response_number);
                response_number++;
                pthread_mutex_unlock(&lock_response_number);
                pthread_mutex_unlock(&lock_strategyIsReady);

                fprintf(stdout, "thread %d:waiting for start signal...\n", id);
                pthread_mutex_lock(&lock_start_strategy);
                while(!start_flag){
                    pthread_cond_wait(&cond_start_strategy, &lock_start_strategy);
                }
                fprintf(stdout, "thread %d:start strategy!\n", id);
                pthread_mutex_unlock(&lock_start_strategy);
                continue;
            }

        }else{//status = normal
            if(id==RC_id && RC_flag == 1){
                fprintf(stdout,"start recomputing...\n");
                //clear status
                isStraggler = 0;
                RC_flag = 0;
                RC_id = -1;
                times_over_threshold = 0;
                times_below_threshold = 0;
                pthread_mutex_lock(&lock_hasStraggler);
                hasStraggler = 0;
                currentStragger = -1;
                pthread_mutex_unlock(&lock_hasStraggler);
                /***decode***/
                int j;
                if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                {
                    for (j = 0; j < RC_count; j++)
                    {
                        //decode_startTime = GetTimeStamp();
                        jerasure_matrix_decode(k, m, w, sample_matrix, 1, sample_erasures, sample_data, sample_coding, ec_blocksize);
                        //decode_endTime = GetTimeStamp();
                        //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                    }
                }
                else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                {
                    for (j = 0; j < RC_count; j++)
                    {
                        //decode_startTime = GetTimeStamp();
                        jerasure_schedule_decode_lazy(k, m, w, sample_bitmatrix, sample_erasures, sample_data, sample_coding, ec_blocksize, ec_packetsize, 1);
                        //decode_endTime = GetTimeStamp();
                        //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                    }
                }
                fprintf(stdout, "thread %d decode complete, current ready offset: %lld\n", id, pairCnt);
                /***decode***/
                pairCnt = next_pairCnt; //keep up with the latest
                continue;
            }else if(id<k){
                ;
            }else if(id>=k){
                if(parity_hangup[id-k]==0){
                    ;
                }else{
                    parity_hangup[id-k] = 0;
                    pthread_mutex_lock(&lock_hasStraggler);
                    hasStraggler = 0;
                    fprintf(stdout,"thread %d hang up until next straggler...\n", id);
                    while(!hasStraggler){
                        pthread_cond_wait(&cond_hasStraggler, &lock_hasStraggler);
                    }
                    pthread_mutex_unlock(&lock_hasStraggler);
                    fprintf(stdout, "thread %d: wake up!\n", id);
                    continue;
                }
            }
        }
        /*****************collective oeration*********************/

        if (id < k)
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_data)[id], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        else
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_coding)[id - k], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        xfer_endTime = GetTimeStamp() - startTime;
        duration = xfer_endTime - xfer_startTime;
        if (pairCnt == 8192)
        {
            leftThreads--;
            fprintf(stdout, "thread %d duration: %0.4lf pairCnt = %lld\n", id, duration, pairCnt);
        }

        // dataLeft[id]--;
    }
    //fprintf(stdout, "reading file%d complete, size: %lld\n", id, transferred_size);

    // pthread_mutex_lock(&lockOfNT);
    // tranferDone[id] = 1;
    // numTransferred += 1;
    // pthread_mutex_unlock(&lockOfNT);
    endTime = GetTimeStamp();
    transferTime = endTime - startTime;
    ec_timers[id] += transferTime;
}

void *
ec_collective_thread3(ec_read_thread_args *arg) //violate
{
    int id = arg->id;
    int k = arg->test->ec_k;
    int m = arg->test->ec_m;
    int w = arg->test->ec_w;
    enum Coding_Technique method = arg->test->ec_method;
    IOR_offset_t ec_blocksize = arg->test->ec_stripe_size;
    int ec_packetsize = arg->test->ec_packetsize;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;

    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    IOR_offset_t pairCnt = 0;
    double startTime = 0;
    double endTime = 0;

    int isOriginReader = id < k ? 1 : 0;
    int isStraggler = 0;
    int transferTime = 0;
    double xfer_startTime = 0;
    double xfer_endTime = 0;
    double duration = 0;
    double times_over_threshold = 0;
    double times_below_threshold = 0;
    //double threshold = 0.02;
    double upper_threshold = 0.02;
    double lower_threshold = 0.01;
    startTime = GetTimeStamp();

    //fprintf(stdout, "break at 3524\n");
    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {
        /*****************collective oeration*********************/
        /****************is_straggler******************/
        if (duration > upper_threshold && !isStraggler)
        { //不是straggler，等待进入
            times_over_threshold++;
        }
        else if (duration < lower_threshold && isStraggler)
        { //是straggler，等待退出
            times_below_threshold++;
        }
        else if (duration <= lower_threshold && !isStraggler)
        { //不是straggler，进入状态清零
            times_over_threshold = 0;
        }
        else if (duration >= upper_threshold && isStraggler)
        { //是straggler，退出状态清零
            times_below_threshold = 0;
        }

        if (times_over_threshold >= 10)
        {
            isStraggler = 1;
            RC_id = id;
            times_over_threshold = 0;
            fprintf(stdout, "process %d: thread %d into straggler, offset: %lld\n", rank,id, pairCnt);
        }
        if (times_below_threshold >= 5)
        {
            isStraggler = 0;
            times_below_threshold = 0;
            if (currentStragger == id)
            {
                pthread_mutex_lock(&lock_hasStraggler);
                hasStraggler = 0;
                currentStragger = -1;
                pthread_mutex_unlock(&lock_hasStraggler);
            }
            fprintf(stdout, "process %d:thread %d quit straggler, offset: %lld\n", rank,id, pairCnt);
        }
        /****************is_straggler******************/

        /*****************collective oeration*********************/
        if(id==RC_id)
            fprintf(stdout,"process %d:thread %d start compute from %lld\n", rank,id, pairCnt);

        while(id == RC_id){
        
            /***decode***/
            if(pairCnt == num_reconstruct){
                fprintf(stdout,"process %d:thread %d end with compute %lld\n", rank,id, pairCnt);
                break;
            }

            // MPI_Reduce(&local_finished,&global_finished, 1, MPI_INT, MPI_SUM, rank, testComm);
            // fprintf(stdout,"global_finished = %d\n", global_finished);

            // if(global_finished == numTasksWorld-1){
            //     fprintf(stdout,"process %d's file can be immediatly reconstructed\n", rank);
            // }
            //fprintf(stdout,"process %d:thread %d current pairCnt =  %lld\n", rank,id, pairCnt);
            
            if(pairCnt < kth_large_offset_of_stripes(k+m,k)){
                if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                {
                
                    //decode_startTime = GetTimeStamp();
                    jerasure_matrix_decode(k, m, w, sample_matrix, 1, sample_erasures, sample_data, sample_coding, ec_blocksize);
                    //decode_endTime = GetTimeStamp();
                    //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                
                }
                else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                {
                
                    //decode_startTime = GetTimeStamp();
                    jerasure_schedule_decode_lazy(k, m, w, sample_bitmatrix, sample_erasures, sample_data, sample_coding, ec_blocksize, ec_packetsize, 1);
                    //decode_endTime = GetTimeStamp();
                    //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                
                }
                pairCnt++;
                currentPosOfThread[id]++;
                //fprintf(stdout, "thread %d decode complete, current ready offset: %lld\n", id, pairCnt);
                /***decode***/
            }
            
            
        }

        if(pairCnt == num_reconstruct){
            fprintf(stdout, "process %d:thread %d end with compute2 %lld\n", rank, id, pairCnt);
            break;
        }

        offset = offsetArray[pairCnt];
        offset = offset / arg->test->ec_k;

        //fprintf(stdout, "thread %d stop compute, current offset %lld\n", id, pairCnt);
        /*****************collective oeration*********************/
        
        
        if (id < k)
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_data)[id], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        else
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_coding)[id - k], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        xfer_endTime = GetTimeStamp() - startTime;
        duration = xfer_endTime - xfer_startTime;
        if (pairCnt == num_reconstruct)
        {
            //leftThreads--;
            fprintf(stdout, "process %d:thread %d end with transfer\n", rank, id);
            //fprintf(stdout, "process %d: thread %d duration: %0.4lf pairCnt = %lld\n", rank, id, duration, pairCnt);
        }

        // dataLeft[id]--;
    }
    //fprintf(stdout, "reading file%d complete, size: %lld\n", id, transferred_size);

    // pthread_mutex_lock(&lockOfNT);
    // tranferDone[id] = 1;
    // numTransferred += 1;
    // pthread_mutex_unlock(&lockOfNT);
    //fprintf(stdout, "process %d: thread %d kth =  %lld\n", rank, id, kth_large_offset_of_stripes(k + m, k));
    //fprintf(stdout, "process %d: thread %d  pairCnt = %lld\n", rank, id, pairCnt);
    //fprintf(stdout, "process %d: num_reconstruct = %lld\n", rank, num_reconstruct);
    finished_thread++;
    if(id==0 && finished_thread==(k+m)){
        *(arg->local_finished) = 1;
        fprintf(stdout, "process %d: local_finished = %d\n", rank, *(arg->local_finished));
    }
    endTime = GetTimeStamp();
    transferTime = endTime - startTime;
    ec_timers[id] += transferTime;
}

void *
ec_collective_thread4(ec_read_thread_args *arg)
{
    int id = arg->id;
    int k = arg->test->ec_k;
    int m = arg->test->ec_m;
    int w = arg->test->ec_w;
    enum Coding_Technique method = arg->test->ec_method;
    IOR_offset_t ec_blocksize = arg->test->ec_stripe_size;
    int ec_packetsize = arg->test->ec_packetsize;
    //void *fd = arg->fds->fd[id];
    //fprintf(stdout,"reading file%...\n",id);

    IOR_offset_t transferred_size = 0;
    IOR_offset_t *offsetArray = arg->offSetArray;
    IOR_offset_t offset;

    int num_reconstruct = (arg->test->blockSize / arg->test->transferSize) * arg->test->segmentCount;
    dataLeft[id] = num_reconstruct;
    IOR_offset_t pairCnt = 0;
    double startTime = 0;
    double endTime = 0;

    int isOriginReader = id < k ? 1 : 0;
    int isStraggler = 0;
    int transferTime = 0;
    double xfer_startTime = 0;
    double xfer_endTime = 0;
    double duration = 0;
    double times_over_threshold = 0;
    double times_below_threshold = 0;
    //double threshold = 0.02;
    double upper_threshold = 0.02;
    double lower_threshold = 0.01;
    startTime = GetTimeStamp();

    //fprintf(stdout, "break at 3524\n");
    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {

        /*****************collective oeration*********************/
        /****************is_straggler******************/
        if (duration > upper_threshold && !isStraggler)
        { //不是straggler，等待进入
            times_over_threshold++;
        }
        else if (duration < lower_threshold && isStraggler)
        { //是straggler，等待退出
            times_below_threshold++;
        }
        else if (duration <= lower_threshold && !isStraggler)
        { //不是straggler，进入状态清零
            times_over_threshold = 0;
        }
        else if (duration >= upper_threshold && isStraggler)
        { //是straggler，退出状态清零
            times_below_threshold = 0;
        }

        if (times_over_threshold >= 10)
        {
            isStraggler = 1;
            RC_id = id;
            times_over_threshold = 0;
            fprintf(stdout, "process %d: thread %d into straggler, offset: %lld\n", rank, id, pairCnt);
        }
        if (times_below_threshold >= 5)
        {
            isStraggler = 0;
            times_below_threshold = 0;
            if (currentStragger == id)
            {
                pthread_mutex_lock(&lock_hasStraggler);
                hasStraggler = 0;
                currentStragger = -1;
                pthread_mutex_unlock(&lock_hasStraggler);
            }
            fprintf(stdout, "process %d:thread %d quit straggler, offset: %lld\n", rank, id, pairCnt);
        }
        /****************is_straggler******************/

        /*****************collective oeration*********************/
        if (id == RC_id)
            fprintf(stdout, "process %d:thread %d start compute from %lld\n", rank, id, pairCnt);

        while (id == RC_id)
        {

            /***decode***/
            if (pairCnt == num_reconstruct)
            {
                fprintf(stdout, "process %d:thread %d end with compute %lld\n", rank, id, pairCnt);
                break;
            }

            //fprintf(stdout,"process %d:thread %d current pairCnt =  %lld\n", rank,id, pairCnt);

            if (pairCnt < kth_large_offset_of_stripes(k + m, k))
            {
                if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                {

                    //decode_startTime = GetTimeStamp();
                    jerasure_matrix_decode(k, m, w, sample_matrix, 1, sample_erasures, sample_data, sample_coding, ec_blocksize);
                    //decode_endTime = GetTimeStamp();
                    //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                }
                else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                {

                    //decode_startTime = GetTimeStamp();
                    jerasure_schedule_decode_lazy(k, m, w, sample_bitmatrix, sample_erasures, sample_data, sample_coding, ec_blocksize, ec_packetsize, 1);
                    //decode_endTime = GetTimeStamp();
                    //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                }
                pairCnt++;
                currentPosOfThread[id]++;
                //fprintf(stdout, "thread %d decode complete, current ready offset: %lld\n", id, pairCnt);
                /***decode***/
            }
        }

        if (pairCnt == num_reconstruct)
        {
            fprintf(stdout, "process %d:thread %d end with compute2 %lld\n", rank, id, pairCnt);
            break;
        }

        offset = offsetArray[pairCnt];
        offset = offset / arg->test->ec_k;

        //fprintf(stdout, "thread %d stop compute, current offset %lld\n", id, pairCnt);
        /*****************collective oeration*********************/

        if (id < k)
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_data)[id], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        else
        {
            xfer_startTime = GetTimeStamp() - startTime;
            transferred_size = IOR_Xfer_ec(arg->access, (arg->fds)[id], (arg->ec_coding)[id - k], arg->test->ec_stripe_size, arg->test, offset);
            pairCnt++;
            currentPosOfThread[id]++;
        }
        xfer_endTime = GetTimeStamp() - startTime;
        duration = xfer_endTime - xfer_startTime;
        if (pairCnt == num_reconstruct)
        {
            //leftThreads--;
            fprintf(stdout, "process %d:thread %d end with transfer\n", rank, id);
            //fprintf(stdout, "process %d: thread %d duration: %0.4lf pairCnt = %lld\n", rank, id, duration, pairCnt);
        }

        // dataLeft[id]--;
    }
    //fprintf(stdout, "reading file%d complete, size: %lld\n", id, transferred_size);

    // pthread_mutex_lock(&lockOfNT);
    // tranferDone[id] = 1;
    // numTransferred += 1;
    // pthread_mutex_unlock(&lockOfNT);
    fprintf(stdout, "process %d: thread %d kth =  %lld\n", rank, id, kth_large_offset_of_stripes(k + m, k));
    fprintf(stdout, "process %d: thread %d  pairCnt = %lld\n", rank, id, pairCnt);
    fprintf(stdout, "process %d: num_reconstruct = %lld\n", rank, num_reconstruct);
    endTime = GetTimeStamp();
    transferTime = endTime - startTime;
    ec_timers[id] += transferTime;
}

void *
ec_decode_thread(int *target){
    
    //double decode_startTime;
    //double decode_endTime;
    enum Coding_Technique method = ec_decode_arg.method;
    //int decode_res_local = 0;
    //decode_startTime = GetTimeStamp();
    if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
    {
        decode_res = jerasure_matrix_decode(ec_decode_arg.k, ec_decode_arg.m, ec_decode_arg.w, ec_decode_arg.ec_matrix, 1, ec_decode_arg.erasures, omp_data[*target], omp_coding[*target], ec_decode_arg.ec_blocksize);
    }
    else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
    {
        decode_res = jerasure_schedule_decode_lazy(ec_decode_arg.k, ec_decode_arg.m, ec_decode_arg.w, ec_decode_arg.ec_bitmatrix, ec_decode_arg.erasures, omp_data[*target], omp_coding[*target], ec_decode_arg.ec_blocksize, ec_decode_arg.ec_packetsize, 1);
    }
    //decode_endTime = GetTimeStamp();
    //fprintf(stdout,"in decode thread, time is %lf\n", decode_endTime-decode_startTime);
    pthread_mutex_lock(&lockOfDecodeSum);
    decode_sum++;
    pthread_mutex_unlock(&lockOfDecodeSum);
}

IOR_offset_t
WriteOrRead_CL(IOR_param_t *test,
               void **ec_fds,
               int access)
{
    if (test->ec_verbose >= VERBOSE_3)
    {
        fprintf(stdout, "process %d, in WriteOrRead_CL\n", rank);
    }
    int errors = 0;
    IOR_offset_t amtXferred,
        transfer,
        transferCount = 0,
        pairCnt = 0,
        *offsetArray;
    int pretendRank;
    void *buffer = NULL;
    void *checkBuffer = NULL;
    void *readCheckBuffer = NULL;
    IOR_offset_t dataMoved = 0; /* for data rate calculation */
    double startForStonewall;
    //int hitStonewall;

    /********************************ec_parameters******************************/
    //int readins = 1;
    int total_stripe_num = test->ec_k + test->ec_m;

    IOR_offset_t *ec_amtXferred = (IOR_offset_t *)malloc(sizeof(IOR_offset_t) * total_stripe_num);
    IOR_offset_t ec_blocksize = test->ec_stripe_size;
    ec_read_thread_args *ec_read_args;

    Coding_Technique method = test->ec_method; // coding technique (parameter)
    int k = test->ec_k,
        m = test->ec_m,
        w = test->ec_w,
        ec_packetsize = test->ec_packetsize; // parameters
    
    int k1 = numTasksWorld-1;
    int m1 = 1;
    MPI_Status status;

    char **ec_data = NULL;
    char **ec_coding = NULL;
    int *ec_matrix = NULL;
    int *ec_bitmatrix = NULL;
    int **ec_schedule = NULL;

    int i;

    int local_finished = 0;
    int global_finished = 0;
    /********************************ec_parameters******************************/

    /* initialize values */
    pretendRank = (rank + rankOffset) % test->numTasks;

    if (test->randomOffset)
    {
        offsetArray = GetOffsetArrayRandom(test, pretendRank, access);
    }
    else
    {
        offsetArray = GetOffsetArraySequential(test, pretendRank);
    }

    SetupXferBuffers(&buffer, &checkBuffer, &readCheckBuffer,
                     test, pretendRank, access);

    //fprintf(stdout, "break at 2907\n");
    /* check for stonewall */
    startForStonewall = GetTimeStamp();
    hitStonewall = ((test->deadlineForStonewalling != 0) && ((GetTimeStamp() - startForStonewall) > test->deadlineForStonewalling));

    /* loop over offsets to access */

    /*****************************init ec********************************/
    if (access == WRITE)
    {
        //prepare ec_data 1 - k  RS(k1,1) 
        char **ec_data0 = (char **)malloc(sizeof(char*) * k1);
        for(i = 0; i<k1; i++){
            ec_data0[i] = (char *)malloc(sizeof(char) * k* ec_blocksize);
            fprintf(stdout, "waiting for MPI_Recv\n");
            MPI_Recv(ec_data0[i], k* ec_blocksize, MPI_CHAR, i, i, testComm, &status);
        }

        fprintf(stdout, "break after MPI_Recv\n");

        int j;
        char ***ec_data1= (char***)malloc(sizeof(char**) * k);
        for(i=0;i<k;i++){
            ec_data1[i] = (char**)malloc(sizeof(char*)*k1);
            for(j=0;j<k1;j++){              
                ec_data1[i][j] = ec_data0[j] + i * ec_blocksize;
                fprintf(stdout, "break after ec_data1\n");
            }
        } 
        char ***ec_coding1 = (char***)malloc(sizeof(char**) * k);
        for(i=0;i<k;i++){
            ec_coding1[i] = (char**)malloc(sizeof(char*) * m1);
            for(j=0;j<m1;j++){
                ec_coding1[i][j] = (char *)malloc(sizeof(char) * ec_blocksize);
                fprintf(stdout, "break after ec_coding1\n");
            }
        } 

        fprintf(stdout, "break after ec_data1 and ec_coding1\n");

        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            ec_matrix = reed_sol_vandermonde_coding_matrix(k1, m1, w);
            break;
        case Reed_Sol_R6_Op:
            break;
        case Cauchy_Orig:
            ec_matrix = cauchy_original_coding_matrix(k1, m1, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k1, m1, w, ec_matrix);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k1, m1, w, ec_bitmatrix);
            break;
        case Cauchy_Good:
            ec_matrix = cauchy_good_general_coding_matrix(k1, m1, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k1, m1, w, ec_matrix);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k1, m1, w, ec_bitmatrix);
            break;
        case Liberation:
            ec_bitmatrix = liberation_coding_bitmatrix(k1, w);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k1, m1, w, ec_bitmatrix);
            break;
        case Blaum_Roth:
            ec_bitmatrix = blaum_roth_coding_bitmatrix(k1, w);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k1, m1, w, ec_bitmatrix);
            break;
        case Liber8tion:
            ec_bitmatrix = liber8tion_coding_bitmatrix(k1);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k1, m1, w, ec_bitmatrix);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }

        for(i=0;i<k;i++){
            switch (method)
            {
            case No_Coding:
                break;
            case Reed_Sol_Van:
                jerasure_matrix_encode(k1, m1, w, ec_matrix, ec_data1[i], ec_coding1[i], ec_blocksize);
                break;
            case Reed_Sol_R6_Op:
                reed_sol_r6_encode(k1, w, ec_data1[i], ec_coding1[i], ec_blocksize);
                break;
            case Cauchy_Orig:
                jerasure_schedule_encode(k1, m1, w, ec_schedule, ec_data1[i], ec_coding1[i], ec_blocksize, ec_packetsize);
                break;
            case Cauchy_Good:
                jerasure_schedule_encode(k1, m1, w, ec_schedule, ec_data1[i], ec_coding1[i], ec_blocksize, ec_packetsize);
                break;
            case Liberation:
                jerasure_schedule_encode(k1, m1, w, ec_schedule, ec_data1[i], ec_coding1[i], ec_blocksize, ec_packetsize);
                break;
            case Blaum_Roth:
                jerasure_schedule_encode(k1, m1, w, ec_schedule, ec_data1[i], ec_coding1[i], ec_blocksize, ec_packetsize);
                break;
            case Liber8tion:
                jerasure_schedule_encode(k1, m1, w, ec_schedule, ec_data1[i], ec_coding1[i], ec_blocksize, ec_packetsize);
                break;
            case RDP:
            case EVENODD:
                assert(0);
            }
        }
       //end prepare ec_data 1 - k  RS(k1,1) 

        fprintf(stdout, "bread after encode\n");
        
        ec_data = (char **)malloc(sizeof(char *) * k);
        ec_coding = (char **)malloc(sizeof(char *) * m);
        for (i = 0; i < m; i++)
        {
            ec_coding[i] = (char *)malloc(sizeof(char) * test->ec_stripe_size);
            if (ec_coding[i] == NULL)
            {
                ERR("malloc coding fail");
            }
        }
        //init matrix
        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            ec_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
            break;
        case Reed_Sol_R6_Op:
            break;
        case Cauchy_Orig:
            ec_matrix = cauchy_original_coding_matrix(k, m, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
            break;
        case Cauchy_Good:
            ec_matrix = cauchy_good_general_coding_matrix(k, m, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
            break;
        case Liberation:
            ec_bitmatrix = liberation_coding_bitmatrix(k, w);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
            break;
        case Blaum_Roth:
            ec_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
            break;
        case Liber8tion:
            ec_bitmatrix = liber8tion_coding_bitmatrix(k);
            ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        //fprintf(stdout, "break at 2963\n");
        for (i = 0; i < k; i++)
        {
            ec_data[i] = ec_coding1[i][0];
        }

        fprintf(stdout, "break after ec_data[i] = ec_coding1[i][0]\n");
        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            jerasure_matrix_encode(k, m, w, ec_matrix, ec_data, ec_coding, ec_blocksize);
            break;
        case Reed_Sol_R6_Op:
            reed_sol_r6_encode(k, w, ec_data, ec_coding, ec_blocksize);
            break;
        case Cauchy_Orig:
            jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
            break;
        case Cauchy_Good:
            jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
            break;
        case Liberation:
            jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
            break;
        case Blaum_Roth:
            jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
            break;
        case Liber8tion:
            jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        fprintf(stdout, "break after ec_data,ec_coding\n");
    }
    else if (access == READ)
    {

        /******prepare sample_data & sample_coding********/
        sample_data = (char **)malloc(sizeof(char *) * k);
        sample_coding = (char **)malloc(sizeof(char *) * m);
        for (i = 0; i < m; i++)
        {
            sample_coding[i] = (char *)malloc(sizeof(char) * test->ec_stripe_size);
            if (sample_coding[i] == NULL)
            {
                ERR("malloc sample_coding fail");
            }
        }
        //init matrix
        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            sample_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
            break;
        case Reed_Sol_R6_Op:
            break;
        case Cauchy_Orig:
            sample_matrix = cauchy_original_coding_matrix(k, m, w);
            sample_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, sample_matrix);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Cauchy_Good:
            sample_matrix = cauchy_good_general_coding_matrix(k, m, w);
            sample_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, sample_matrix);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Liberation:
            sample_bitmatrix = liberation_coding_bitmatrix(k, w);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Blaum_Roth:
            sample_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Liber8tion:
            sample_bitmatrix = liber8tion_coding_bitmatrix(k);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        //fprintf(stdout, "break at 2963\n");
        for (i = 0; i < k; i++)
        {
            sample_data[i] = buffer + (i * ec_blocksize);
        }

        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            jerasure_matrix_encode(k, m, w, sample_matrix, sample_data, sample_coding, ec_blocksize);
            break;
        case Reed_Sol_R6_Op:
            reed_sol_r6_encode(k, w, sample_data, sample_coding, ec_blocksize);
            break;
        case Cauchy_Orig:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Cauchy_Good:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Liberation:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Blaum_Roth:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Liber8tion:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        /******prepare sample ec-data********/

        ec_startTime = GetTimeStamp();

        ec_timers = (double *)malloc(sizeof(double) * total_stripe_num);
        threads = (pthread_t *)malloc(sizeof(pthread_t) * total_stripe_num);
        ec_read_args = (ec_read_thread_args *)malloc(sizeof(ec_read_thread_args) * total_stripe_num);
        tranferDone = (int *)malloc(sizeof(int) * total_stripe_num);
        memset(tranferDone, 0, sizeof(int) * total_stripe_num);
        dataLeft = (int *)malloc(sizeof(int) * total_stripe_num);
        memset(tranferDone, 0, sizeof(int) * total_stripe_num);

        currentPosOfThread = (IOR_offset_t *)malloc(sizeof(IOR_offset_t) * total_stripe_num);
        memset(currentPosOfThread, 0, sizeof(IOR_offset_t) * total_stripe_num);
        strategyReceived = (int *)malloc(sizeof(int) * total_stripe_num);
        memset(strategyReceived, 0, sizeof(int) * total_stripe_num);
        parity_hangup = (int *)malloc(sizeof(int) * m);
        memset(parity_hangup, 0, sizeof(int) * m);
        leftThreads = total_stripe_num;

        ec_data = (char **)malloc(sizeof(char *) * k);
        ec_coding = (char **)malloc(sizeof(char *) * m);

        for (i = 0; i < k; i++)
        {
            ec_data[i] = (char *)malloc(sizeof(char) * ec_blocksize);
            if (ec_data[i] == NULL)
            {
                ERR("malloc data fail");
            }
        }

        for (i = 0; i < m; i++)
        {
            ec_coding[i] = (char *)malloc(sizeof(char) * ec_blocksize);
            if (ec_coding[i] == NULL)
            {
                ERR("malloc coding fail");
            }
        }

        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            ec_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
            break;
        case Reed_Sol_R6_Op:
            ec_matrix = reed_sol_r6_coding_matrix(k, w);
            break;
        case Cauchy_Orig:
            ec_matrix = cauchy_original_coding_matrix(k, m, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
            break;
        case Cauchy_Good:
            ec_matrix = cauchy_good_general_coding_matrix(k, m, w);
            ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
            break;
        case Liberation:
            ec_bitmatrix = liberation_coding_bitmatrix(k, w);
            break;
        case Blaum_Roth:
            ec_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
            break;
        case Liber8tion:
            ec_bitmatrix = liber8tion_coding_bitmatrix(k);
        }

        for (i = 0; i < total_stripe_num; i++)
        {
            ec_read_args[i].fds = ec_fds;
            ec_read_args[i].id = i;
            ec_read_args[i].test = test;
            ec_read_args[i].access = access;
            ec_read_args[i].ec_data = ec_data;
            ec_read_args[i].ec_coding = ec_coding;
            ec_read_args[i].method = method;
            ec_read_args[i].ec_matrix = ec_matrix;
            ec_read_args[i].ec_bitmatrix = ec_bitmatrix;
            ec_read_args[i].offSetArray = offsetArray;
            ec_read_args[i].local_finished = &local_finished;
        }
    }

    //fprintf(stdout, "break at 30708\n");
    /*****************************init ec********************************/

    /*****************************ec_read********************************/
    if (access == READ)
    {
        numTransferred = 0;

        if (test->ec_strategy == COLLECTIVE_THREAD)
        {
            for (i = 0; i < total_stripe_num; i++)
            {
                pthread_create(&threads[i], NULL, ec_collective_thread3, &ec_read_args[i]);
            }
        }
        else
        {
            for (i = 0; i < total_stripe_num; i++)
            {
                pthread_create(&threads[i], NULL, ec_read_thread, &ec_read_args[i]);
            }
        }

        /*ec when k stripes arrive*/
        fprintf(stdout, "ec_strategy = %d\n", test->ec_strategy);
        if (test->ec_strategy == IMMEDIATE_EC)
        {
            if (test->ec_verbose >= VERBOSE_0)
            {
                fprintf(stdout, "using immediate EC strategy!\n");
            }
            while (1)
            {
                pthread_mutex_lock(&lockOfNT);
                if (numTransferred == k)
                {
                    ec_strategy_startTime = GetTimeStamp();
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "when k stripes arrive:\n");
                        for (i = 0; i < (k + m); i++)
                        {
                            fprintf(stdout, "stripe %d: %d\n", i, tranferDone[i]);
                        }
                    }
                    /*end unfinished jobs*/
                    int canceled;
                    /*end unfinished jobs*/
                    int originReceived = 1;
                    int *erasures = (int *)malloc(sizeof(int) * total_stripe_num);
                    for (i = 0; i < k; i++)
                    { //if original stripes received or not
                        if (tranferDone[i] == 0)
                        {
                            originReceived = 0;
                            break;
                        }
                    }
                    if (originReceived)
                    {
                        for (i = 0; i < m; i++)
                        {
                            canceled = pthread_cancel(threads[k + i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to oringin stripes received\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (origin received)\n", threads[k + i]);
                            }
                        }
                        pthread_mutex_unlock(&lockOfNT);
                        break;
                    }

                    /*recompute*/
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "in recompute phase:\n");
                    }
                    int numerased = 0;
                    int num_iteration = 0;

                    for (i = 0; i < total_stripe_num; i++)
                    {
                        if (tranferDone[i] == 0)
                        {
                            canceled = pthread_cancel(threads[i]);
                            num_iteration = MAX(num_iteration, dataLeft[i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to latency\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (latency)\n", threads[k + i]);
                            }
                            erasures[numerased] = i;
                            numerased++;
                        }
                    }
                    erasures[numerased] = -1;

                    //reconstruct for njum_reconstruct times

                    int j;
                    double decode_startTime;
                    double decode_endTime;
                    double total_decode_time = 0;
                    if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                    {
                        for (j = 0; j < num_iteration; j++)
                        {
                            decode_startTime = GetTimeStamp();
                            i = jerasure_matrix_decode(k, m, w, ec_matrix, 1, erasures, ec_data, ec_coding, ec_blocksize);
                            decode_endTime = GetTimeStamp();
                            //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                            total_decode_time += decode_endTime - decode_startTime;
                        }
                    }
                    else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                    {
                        for (j = 0; j < num_iteration; j++)
                        {
                            decode_startTime = GetTimeStamp();
                            i = jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, ec_data, ec_coding, ec_blocksize, ec_packetsize, 1);
                            decode_endTime = GetTimeStamp();
                            //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                            total_decode_time += decode_endTime - decode_startTime;
                        }
                    }

                    double average_decode_time = total_decode_time / j;
                    if (i == -1)
                    {
                        pthread_mutex_unlock(&lockOfNT);
                        ERR("decode failed!");
                    }
                    else
                    {
                        if (test->ec_verbose >= VERBOSE_0)
                        {
                            fprintf(stdout, "decode success for %d times!\n", j);
                            fprintf(stdout, "average_decode_time = %lf\n", average_decode_time);
                        }
                    }
                    pthread_mutex_unlock(&lockOfNT);
                    break;
                }
                pthread_mutex_unlock(&lockOfNT);
            }
        }
        else if (test->ec_strategy == POST_PARALLELISM)
        {
            if (test->ec_verbose >= VERBOSE_0)
            {
                fprintf(stdout, "using Post Parallelism Processing!\n");
            }
            while (1)
            {
                pthread_mutex_lock(&lockOfNT);
                if (numTransferred == k)
                {
                    ec_strategy_startTime = GetTimeStamp();
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "when k stripes arrive:\n");
                        for (i = 0; i < (k + m); i++)
                        {
                            fprintf(stdout, "stripe %d: %d\n", i, tranferDone[i]);
                        }
                    }
                    /*end unfinished jobs*/
                    int canceled;
                    /*end unfinished jobs*/
                    int originReceived = 1;
                    int *erasures = (int *)malloc(sizeof(int) * total_stripe_num);
                    for (i = 0; i < k; i++)
                    { //if original stripes received or not
                        if (tranferDone[i] == 0)
                        {
                            originReceived = 0;
                            break;
                        }
                    }
                    if (originReceived)
                    {
                        for (i = 0; i < m; i++)
                        {
                            canceled = pthread_cancel(threads[k + i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to oringin stripes received\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (origin received)\n", threads[k + i]);
                            }
                        }
                        pthread_mutex_unlock(&lockOfNT);
                        break;
                    }

                    /*recompute*/
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "in recompute phase:\n");
                    }
                    int numerased = 0;
                    int num_iteration = 0;

                    for (i = 0; i < total_stripe_num; i++)
                    {
                        if (tranferDone[i] == 0)
                        {
                            canceled = pthread_cancel(threads[i]);
                            num_iteration = MAX(num_iteration, dataLeft[i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to latency\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (latency)\n", threads[k + i]);
                            }
                            erasures[numerased] = i;
                            numerased++;
                        }
                    }
                    erasures[numerased] = -1;

                    //reconstruct for njum_reconstruct times
                    int j;

                    omp_data = (char ***)malloc(sizeof(char **) * omp_thread_num);
                    omp_coding = (char ***)malloc(sizeof(char **) * omp_thread_num);

                    for (j = 0; j < omp_thread_num; j++)
                    {
                        omp_data[j] = (char **)malloc(sizeof(char *) * k);
                        omp_coding[j] = (char **)malloc(sizeof(char *) * m);

                        for (i = 0; i < k; i++)
                        {
                            omp_data[j][i] = (char *)malloc(sizeof(char) * ec_blocksize);
                            memcpy(omp_data[j][i], ec_data[i], sizeof(char) * ec_blocksize);
                        }

                        for (i = 0; i < m; i++)
                        {
                            omp_coding[j][i] = (char *)malloc(sizeof(char) * ec_blocksize);
                            memcpy(omp_coding[j][i], ec_coding[i], sizeof(char) * ec_blocksize);
                        }
                    }

                    /****************************thread_pool*******************************/
                    ec_decode_arg.method = method;
                    ec_decode_arg.k = k;
                    ec_decode_arg.m = m;
                    ec_decode_arg.w = w;
                    ec_decode_arg.ec_packetsize = ec_packetsize;
                    ec_decode_arg.ec_blocksize = ec_blocksize;
                    ec_decode_arg.ec_matrix = ec_matrix;
                    ec_decode_arg.ec_bitmatrix = ec_bitmatrix;
                    ec_decode_arg.erasures = erasures;

                    pool_init(omp_thread_num);

                    for (j = 0; j < num_iteration; j++)
                    {
                        int target = j % omp_thread_num;
                        //fprintf(stdout, "add worker\n");
                        pool_add_worker(ec_decode_thread, &target);
                    }

                    while (1)
                    {
                        pthread_mutex_lock(&lockOfDecodeSum);
                        if (decode_sum == num_iteration)
                        {
                            pthread_mutex_unlock(&lockOfDecodeSum);
                            break;
                        }
                        pthread_mutex_unlock(&lockOfDecodeSum);
                    }
                    pool_destroy();
                    /****************************thread_pool*******************************/

                    /****************************omp*******************************/
                    // int omp_id=0;
                    // i=0;
                    // double omp_startTime;
                    // double omp_endTime;
                    // int num_total=0;
                    // //i = jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize, ec_packetsize, 1);
                    // fprintf(stdout,"ec_method=%d\n", method);
                    // if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                    // {
                    //     #pragma omp parallel for reduction(+:i) num_threads(omp_thread_num) private(omp_id)
                    //     for (j = 0; j < num_iteration; j++){
                    //         omp_id = omp_get_thread_num();
                    //         //fprintf(stdout,"ompid = %d\n", omp_id);
                    //         i = jerasure_matrix_decode(k, m, w, ec_matrix, 1, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize);
                    //     }

                    // }
                    // else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                    // {
                    //     #pragma omp parallel for num_threads(omp_thread_num) //private(omp_id,omp_startTime,omp_endTime)
                    //     for (j = 0; j < num_iteration; j++){
                    //         //omp_id = omp_get_thread_num();
                    //         //fprintf(stdout,"ompid = %d\n", omp_id);
                    //         //omp_startTime = GetTimeStamp();
                    //         jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize, ec_packetsize, 1);
                    //         //omp_endTime = GetTimeStamp();
                    //         num_total++;
                    //         //fprintf(stdout, "time = %lf\n", omp_endTime-omp_startTime);
                    //     }
                    // }
                    /****************************omp*******************************/

                    if (decode_res != 0)
                    {
                        pthread_mutex_unlock(&lockOfNT);
                        ERR("decode failed!");
                    }
                    else
                    {
                        if (test->ec_verbose >= VERBOSE_0)
                        {
                            fprintf(stdout, "decode success for %d times!\n", decode_sum);
                        }
                    }
                    pthread_mutex_unlock(&lockOfNT);
                    break;
                }
                pthread_mutex_unlock(&lockOfNT);
            }
        }

        ec_strategy_endTime = GetTimeStamp();
        
        while(global_finished!=(numTasksWorld-1)){
            MPI_Allreduce(&local_finished,&global_finished,1,MPI_INT,MPI_SUM,testComm);
            fprintf("process %d, current global finished: %d\n", rank, global_finished);
            sleep(1);
        }

        if(rank == 0 && global_finished == (numTasksWorld-1)){
            fprintf("process %d, only lack one file", rank, global_finished);
        }

        for (i = 0; i < total_stripe_num; i++)
        {
            pthread_join(threads[i], NULL);
        }

        /*ec when k stripes arrive*/
        amtXferred = ec_blocksize * k;
        /*************************ec multi-thread read*************************/

        // amtXferred = IOR_Xfer(access, fd, buffer, transfer, test);
        // if (amtXferred != transfer)
        //     ERR("cannot read from file");
        ec_endTime = GetTimeStamp();
    }

    while ((offsetArray[pairCnt] != -1) && !hitStonewall)
    {

        test->offset = offsetArray[pairCnt];
        test->offset = test->offset / test->ec_k;
        //ec_offset
        /*
         * fills each transfer with a unique pattern
         * containing the offset into the file
         */
        if (test->storeFileOffset == TRUE)
        {
            FillBuffer(buffer, test, test->offset, pretendRank);
        }
        transfer = test->transferSize;
        if (access == WRITE)
        {
            /*ec transfer*/
            // for(i=0;i<total_stripe_num;i++){
            //     ec_amtXferred[i] = IOR_Xfer(access, ec_fds[i], ec_data[0], ec_blocksize, test);
            // }
            for (i = 0; i < k; i++)
            {
                ec_amtXferred[i] = IOR_Xfer(access, ec_fds[i], ec_data[i], ec_blocksize, test);
            }
            for (i = 0; i < m; i++)
            {
                ec_amtXferred[k + i] = IOR_Xfer(access, ec_fds[k + i], ec_coding[i], ec_blocksize, test);
            }
            /*ec transfer*/
            //amtXferred = IOR_Xfer(access, fd, buffer, transfer, test); //origin_mark
            amtXferred = ec_blocksize * k;
            for (i = 0; i < total_stripe_num; i++)
            {
                if (ec_amtXferred[i] != ec_blocksize)
                    ERR("write to file failed!");
            }

            //fprintf(stdout, "break at 3102\n");
        }
        else if (access == READ)
        {
            amtXferred = test->ec_stripe_size * test->ec_k;
        }
        else if (access == WRITECHECK)
        {
            // memset(checkBuffer, 'a', transfer);
            // amtXferred = IOR_Xfer(access, fd, checkBuffer, transfer, test);
            // if (amtXferred != transfer)
            //     ERR("cannot read from file write check");
            // transferCount++;
            // errors += CompareBuffers(buffer, checkBuffer, transfer,
            //                          transferCount, test, WRITECHECK);
        }
        else if (access == READCHECK)
        {
            fprintf(stdout, "in READCHECK");
            // ReadCheck(fd, buffer, checkBuffer, readCheckBuffer, test,
            //           transfer, test->blockSize, &amtXferred,
            //           &transferCount, access, &errors);
        }
        dataMoved += amtXferred;
        pairCnt++;

        hitStonewall = ((test->deadlineForStonewalling != 0) && ((GetTimeStamp() - startForStonewall) > test->deadlineForStonewalling));
    }

print_info:
    /*print ec time info*/

    if (access == READ)
    {
        for (i = 0; i < total_stripe_num; i++)
        {
            fprintf(stdout, "#read time of stripe %d : %lf\n", i, ec_timers[i]);
        }
        fprintf(stdout, "#Total read time of %d stripes: %lf\n", total_stripe_num, ec_endTime - ec_startTime);
        fprintf(stdout, "#EC strategy time of %d : %lf\n", test->ec_strategy, ec_strategy_endTime - ec_strategy_startTime);
    }

    totalErrorCount += CountErrors(test, access, errors);

    FreeBuffers(access, checkBuffer, readCheckBuffer, buffer, offsetArray);

    if (access == WRITE && test->fsync == TRUE)
    {
        //IOR_Fsync(fd, test); /*fsync after all accesses*/
        for (i = 0; i < total_stripe_num; i++)
        {
            IOR_Fsync(ec_fds[i], test);
        }
    }
    free(ec_data);
    free(ec_coding);

    return (dataMoved);
} /* WriteOrRead_CL() */

IOR_offset_t
WriteOrRead_ec(IOR_param_t *test,
            void **ec_fds,
            int access)
{
    if(test->ec_verbose >= VERBOSE_3){
        fprintf(stdout, "in WriteOrRead_ec\n");
    }
    int errors = 0;
    IOR_offset_t amtXferred,
        transfer,
        transferCount = 0,
        pairCnt = 0,
        *offsetArray;
    int pretendRank;
    void *buffer = NULL;
    void *checkBuffer = NULL;
    void *readCheckBuffer = NULL;
    IOR_offset_t dataMoved = 0; /* for data rate calculation */
    double startForStonewall;
    //int hitStonewall;

    /********************************ec_parameters******************************/
    //int readins = 1;
    int total_stripe_num = test->ec_k + test->ec_m;

    IOR_offset_t *ec_amtXferred = (IOR_offset_t *)malloc(sizeof(IOR_offset_t) * total_stripe_num);
    IOR_offset_t ec_blocksize = test->ec_stripe_size;
    fprintf(stdout, "ec_blocksize = %ld\n",ec_blocksize);
    ec_read_thread_args *ec_read_args;

    Coding_Technique method = test->ec_method;      // coding technique (parameter)
    int k = test->ec_k, 
        m = test->ec_m, 
        w = test->ec_w, 
        ec_packetsize = test->ec_packetsize; // parameters
    
    char **ec_data = NULL;
    char **ec_coding = NULL;
    int *ec_matrix = NULL;
    int *ec_bitmatrix = NULL;
    int **ec_schedule = NULL;

    int i;
    
    int local_finished = 0;
    int global_finished = 0;
    /********************************ec_parameters******************************/

    /* initialize values */
    pretendRank = (rank + rankOffset) % test->numTasks;

    if (test->randomOffset)
    {
        offsetArray = GetOffsetArrayRandom(test, pretendRank, access);
    }
    else
    {
        offsetArray = GetOffsetArraySequential(test, pretendRank);
    }

    SetupXferBuffers(&buffer, &checkBuffer, &readCheckBuffer,
                     test, pretendRank, access);

    //fprintf(stdout, "break at 2907\n");
    /* check for stonewall */
    startForStonewall = GetTimeStamp();
    hitStonewall = ((test->deadlineForStonewalling != 0) && ((GetTimeStamp() - startForStonewall) > test->deadlineForStonewalling));

    /* loop over offsets to access */
    
    /*****************************init ec********************************/
    if(access == WRITE){
        //allocate space
        ec_data = (char **)malloc(sizeof(char *) * k);
        ec_coding = (char **)malloc(sizeof(char *) * m);
        for (i = 0; i < m; i++)
        {
            ec_coding[i] = (char *)malloc(sizeof(char) * test->ec_stripe_size);
            if (ec_coding[i] == NULL)
            {
                ERR("malloc coding fail");
            }
        }
        //init matrix
        switch (method)
        {
            case No_Coding:
                break;
            case Reed_Sol_Van:
                ec_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
                break;
            case Reed_Sol_R6_Op:
                break;
            case Cauchy_Orig:
                ec_matrix = cauchy_original_coding_matrix(k, m, w);
                ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
                ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
                break;
            case Cauchy_Good:
                ec_matrix = cauchy_good_general_coding_matrix(k, m, w);
                ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
                ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
                break;
            case Liberation:
                ec_bitmatrix = liberation_coding_bitmatrix(k, w);
                ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
                break;
            case Blaum_Roth:
                ec_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
                ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
                break;
            case Liber8tion:
                ec_bitmatrix = liber8tion_coding_bitmatrix(k);
                ec_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, ec_bitmatrix);
                break;
            case RDP:
            case EVENODD:
                assert(0);
        }
        //fprintf(stdout, "break at 2963\n");
        for (i = 0; i < k; i++)
        {
            ec_data[i] = buffer + (i * ec_blocksize);
        }

        // char *sendBuffer = (char*)malloc(sizeof(char)*k*ec_blocksize);
        // for(i=0;i<k;i++){
        //     memcpy(sendBuffer+i*ec_blocksize, ec_data[i], ec_blocksize);
        // }
        // fprintf(stdout,"process %d starts to sending...\n", rank);
        // MPI_Send(sendBuffer, ec_blocksize*k, MPI_CHAR, numTasksWorld-1, rank, testComm);
        
        switch (method)
        {
            case No_Coding:
                break;
            case Reed_Sol_Van:
                jerasure_matrix_encode(k, m, w, ec_matrix, ec_data, ec_coding, ec_blocksize);
                break;
            case Reed_Sol_R6_Op:
                reed_sol_r6_encode(k, w, ec_data, ec_coding, ec_blocksize);
                break;
            case Cauchy_Orig:
                jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
                break;
            case Cauchy_Good:
                jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
                break;
            case Liberation:
                jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
                break;
            case Blaum_Roth:
                jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
                break;
            case Liber8tion:
                jerasure_schedule_encode(k, m, w, ec_schedule, ec_data, ec_coding, ec_blocksize, ec_packetsize);
                break;
            case RDP:
            case EVENODD:
                assert(0);
        }
        fprintf(stdout, "break 2998\n"); 
    }else if(access == READ){

        /******prepare sample_data & sample_coding********/
        sample_data = (char **)malloc(sizeof(char *) * k);
        sample_coding = (char **)malloc(sizeof(char *) * m);
        for (i = 0; i < m; i++)
        {
            sample_coding[i] = (char *)malloc(sizeof(char) * test->ec_stripe_size);
            if (sample_coding[i] == NULL)
            {
                ERR("malloc sample_coding fail");
            }
        }
        //init matrix
        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            sample_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
            break;
        case Reed_Sol_R6_Op:
            break;
        case Cauchy_Orig:
            sample_matrix = cauchy_original_coding_matrix(k, m, w);
            sample_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, sample_matrix);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Cauchy_Good:
            sample_matrix = cauchy_good_general_coding_matrix(k, m, w);
            sample_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, sample_matrix);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Liberation:
            sample_bitmatrix = liberation_coding_bitmatrix(k, w);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Blaum_Roth:
            sample_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case Liber8tion:
            sample_bitmatrix = liber8tion_coding_bitmatrix(k);
            sample_schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, sample_bitmatrix);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        //fprintf(stdout, "break at 2963\n");
        for (i = 0; i < k; i++)
        {
            sample_data[i] = buffer + (i * ec_blocksize);
        }

        switch (method)
        {
        case No_Coding:
            break;
        case Reed_Sol_Van:
            jerasure_matrix_encode(k, m, w, sample_matrix, sample_data, sample_coding, ec_blocksize);
            break;
        case Reed_Sol_R6_Op:
            reed_sol_r6_encode(k, w, sample_data, sample_coding, ec_blocksize);
            break;
        case Cauchy_Orig:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Cauchy_Good:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Liberation:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Blaum_Roth:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case Liber8tion:
            jerasure_schedule_encode(k, m, w, sample_schedule, sample_data, sample_coding, ec_blocksize, ec_packetsize);
            break;
        case RDP:
        case EVENODD:
            assert(0);
        }
        /******prepare sample ec-data********/

        ec_startTime = GetTimeStamp();

        ec_timers = (double *)malloc(sizeof(double) * total_stripe_num);
        threads = (pthread_t *)malloc(sizeof(pthread_t) * total_stripe_num);
        ec_read_args = (ec_read_thread_args *)malloc(sizeof(ec_read_thread_args) * total_stripe_num);
        tranferDone = (int *)malloc(sizeof(int) * total_stripe_num);
        memset(tranferDone,0,sizeof(int) * total_stripe_num);
        dataLeft = (int *)malloc(sizeof(int) * total_stripe_num);
        memset(tranferDone, 0, sizeof(int) * total_stripe_num);

        currentPosOfThread = (IOR_offset_t*)malloc(sizeof(IOR_offset_t)*total_stripe_num);
        memset(currentPosOfThread,0,sizeof(IOR_offset_t)*total_stripe_num);
        strategyReceived = (int *)malloc(sizeof(int)*total_stripe_num);
        memset(strategyReceived,0,sizeof(int)*total_stripe_num);
        parity_hangup = (int*)malloc(sizeof(int) * m);
        memset(parity_hangup,0,sizeof(int)*m);
        leftThreads = total_stripe_num;

        ec_data = (char **)malloc(sizeof(char *) * k);
        ec_coding = (char **)malloc(sizeof(char *) * m);

        for (i = 0; i < k; i++)
        {
            ec_data[i] = (char *)malloc(sizeof(char) * ec_blocksize);
            if (ec_data[i] == NULL)
            {
                ERR("malloc data fail");
            }
        }

        for (i = 0; i < m; i++)
        {
            ec_coding[i] = (char *)malloc(sizeof(char) * ec_blocksize);
            if (ec_coding[i] == NULL)
            {
                ERR("malloc coding fail");
            }
        }

        switch (method)
        {
            case No_Coding:
                break;
            case Reed_Sol_Van:
                ec_matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
                break;
            case Reed_Sol_R6_Op:
                ec_matrix = reed_sol_r6_coding_matrix(k, w);
                break;
            case Cauchy_Orig:
                ec_matrix = cauchy_original_coding_matrix(k, m, w);
                ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
                break;
            case Cauchy_Good:
                ec_matrix = cauchy_good_general_coding_matrix(k, m, w);
                ec_bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, ec_matrix);
                break;
            case Liberation:
                ec_bitmatrix = liberation_coding_bitmatrix(k, w);
                break;
            case Blaum_Roth:
                ec_bitmatrix = blaum_roth_coding_bitmatrix(k, w);
                break;
            case Liber8tion:
                ec_bitmatrix = liber8tion_coding_bitmatrix(k);
        }


        for (i = 0; i < total_stripe_num; i++)
        {
            ec_read_args[i].fds = ec_fds;
            ec_read_args[i].id = i;
            ec_read_args[i].test = test;
            ec_read_args[i].access = access;
            ec_read_args[i].ec_data = ec_data;
            ec_read_args[i].ec_coding = ec_coding;
            ec_read_args[i].method = method;
            ec_read_args[i].ec_matrix = ec_matrix;
            ec_read_args[i].ec_bitmatrix = ec_bitmatrix;
            ec_read_args[i].offSetArray = offsetArray;
            ec_read_args[i].local_finished = &local_finished;
        }
    }

    
    //fprintf(stdout, "break at 30708\n");
    /*****************************init ec********************************/
    

    /*****************************ec_read********************************/
    if(access == READ){
        numTransferred = 0;
        pthread_t testThread;
        pthread_create(&testThread, NULL, ec_request_test, &ec_read_args[0]);
        pthread_join(testThread, NULL);

        if(test->ec_strategy == COLLECTIVE_THREAD){
            for (i = 0; i < total_stripe_num; i++)
            {
                pthread_create(&threads[i], NULL, ec_collective_thread4, &ec_read_args[i]);
            }
        }else if(test->ec_strategy == ADAPTIVE_EC){
            for (i = 0; i < 1; i++)
            {
                pthread_create(&threads[i], NULL, ec_adaptive_thread, &ec_read_args[i]);
            }
        }else{
            for (i = 0; i < total_stripe_num; i++)
            {
                pthread_create(&threads[i], NULL, ec_read_thread, &ec_read_args[i]);
            }
        }
        

        /*ec when k stripes arrive*/
        fprintf(stdout, "ec_strategy = %d\n", test->ec_strategy);
        if (test->ec_strategy == IMMEDIATE_EC)
        {
            if(test->ec_verbose >= VERBOSE_0){
                fprintf(stdout, "using immediate EC strategy!\n");
            }
            while (1)
            {
                pthread_mutex_lock(&lockOfNT);
                if (numTransferred == k)
                {
                    ec_strategy_startTime = GetTimeStamp();
                    if(test->ec_verbose >= VERBOSE_2){
                        fprintf(stdout, "when k stripes arrive:\n");
                        for(i=0;i<(k+m);i++){
                            fprintf(stdout, "stripe %d: %d\n", i, tranferDone[i]);
                        }
                    }
                    /*end unfinished jobs*/
                    int canceled;
                    /*end unfinished jobs*/
                    int originReceived = 1;
                    int *erasures = (int *)malloc(sizeof(int) * total_stripe_num);
                    for (i = 0; i < k; i++)
                    { //if original stripes received or not
                        if (tranferDone[i] == 0)
                        {
                            originReceived = 0;
                            break;
                        }
                    }
                    if (originReceived)
                    {
                        for (i = 0; i < m; i++)
                        {
                            canceled = pthread_cancel(threads[k + i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to oringin stripes received\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (origin received)\n", threads[k + i]);
                            }
                        }
                        pthread_mutex_unlock(&lockOfNT);
                        break;
                    }

                    /*recompute*/
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "in recompute phase:\n");
                    }
                    int numerased = 0;
                    int num_iteration = 0;

                    for (i = 0; i < total_stripe_num; i++)
                    {
                        if (tranferDone[i] == 0)
                        {
                            canceled = pthread_cancel(threads[i]);
                            num_iteration = MAX(num_iteration, dataLeft[i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to latency\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (latency)\n", threads[k + i]);
                            }
                            erasures[numerased] = i;
                            numerased++;
                        }
                    }
                    erasures[numerased] = -1;

                    //reconstruct for njum_reconstruct times
                    
                    int j;
                    double decode_startTime;
                    double decode_endTime;
                    double total_decode_time = 0;
                    if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                    {
                        for(j = 0; j < num_iteration; j++){
                            decode_startTime = GetTimeStamp();
                            i = jerasure_matrix_decode(k, m, w, ec_matrix, 1, erasures, ec_data, ec_coding, ec_blocksize);
                            decode_endTime = GetTimeStamp();
                            //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                            total_decode_time += decode_endTime-decode_startTime;
                        }    
                            
                    }
                    else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                    {
                        for (j = 0; j < num_iteration; j++){
                            decode_startTime = GetTimeStamp();
                            i = jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, ec_data, ec_coding, ec_blocksize, ec_packetsize, 1);
                            decode_endTime = GetTimeStamp();
                            //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
                            total_decode_time += decode_endTime-decode_startTime;
                        }
                    }
                    
                    double average_decode_time = total_decode_time/j;
                    if (i == -1)
                    {
                        pthread_mutex_unlock(&lockOfNT);
                        ERR("decode failed!");
                    }
                    else
                    {
                        if (test->ec_verbose >= VERBOSE_0)
                        {
                            fprintf(stdout, "decode success for %d times!\n", j);
                            fprintf(stdout, "average_decode_time = %lf\n", average_decode_time);
                        }
                    }
                    pthread_mutex_unlock(&lockOfNT);
                    break;
                }
                pthread_mutex_unlock(&lockOfNT);
            }
        }else if (test->ec_strategy == POST_PARALLELISM)
        {
            if (test->ec_verbose >= VERBOSE_0)
            {
                fprintf(stdout, "using Post Parallelism Processing!\n");
            }
            while (1)
            {
                pthread_mutex_lock(&lockOfNT);
                if (numTransferred == k)
                {
                    ec_strategy_startTime = GetTimeStamp();
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "when k stripes arrive:\n");
                        for (i = 0; i < (k + m); i++)
                        {
                            fprintf(stdout, "stripe %d: %d\n", i, tranferDone[i]);
                        }
                    }
                    /*end unfinished jobs*/
                    int canceled;
                    /*end unfinished jobs*/
                    int originReceived = 1;
                    int *erasures = (int *)malloc(sizeof(int) * total_stripe_num);
                    for (i = 0; i < k; i++)
                    { //if original stripes received or not
                        if (tranferDone[i] == 0)
                        {
                            originReceived = 0;
                            break;
                        }
                    }
                    if (originReceived)
                    {
                        for (i = 0; i < m; i++)
                        {
                            canceled = pthread_cancel(threads[k + i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to oringin stripes received\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (origin received)\n", threads[k + i]);
                            }
                        }
                        pthread_mutex_unlock(&lockOfNT);
                        break;
                    }

                    /*recompute*/
                    if (test->ec_verbose >= VERBOSE_2)
                    {
                        fprintf(stdout, "in recompute phase:\n");
                    }
                    int numerased = 0;
                    int num_iteration = 0;

                    for (i = 0; i < total_stripe_num; i++)
                    {
                        if (tranferDone[i] == 0)
                        {
                            canceled = pthread_cancel(threads[i]);
                            num_iteration = MAX(num_iteration, dataLeft[i]);
                            if (canceled == 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d due to latency\n", threads[k + i]);
                            }
                            else if (canceled != 0 && test->ec_verbose >= VERBOSE_2)
                            {
                                fprintf(stdout, "cenceled thread %d failed (latency)\n", threads[k + i]);
                            }
                            erasures[numerased] = i;
                            numerased++;
                        }
                    }
                    erasures[numerased] = -1;

                    //reconstruct for njum_reconstruct times
                    int j;

                    omp_data = (char ***)malloc(sizeof(char**)*omp_thread_num);
                    omp_coding = (char ***)malloc(sizeof(char**)*omp_thread_num);

                    for(j=0;j<omp_thread_num;j++){
                        omp_data[j] = (char**)malloc(sizeof(char*)*k);
                        omp_coding[j] = (char**)malloc(sizeof(char*)*m);

                        for(i=0;i<k;i++){
                            omp_data[j][i] = (char*)malloc(sizeof(char) * ec_blocksize);
                            memcpy(omp_data[j][i], ec_data[i], sizeof(char) * ec_blocksize);
                        }

                        for(i=0;i<m;i++){
                            omp_coding[j][i] = (char*)malloc(sizeof(char) * ec_blocksize);
                            memcpy(omp_coding[j][i],ec_coding[i],sizeof(char)*ec_blocksize);
                        }
                    }

                    /****************************thread_pool*******************************/
                    ec_decode_arg.method = method;
                    ec_decode_arg.k = k;
                    ec_decode_arg.m = m;
                    ec_decode_arg.w = w;
                    ec_decode_arg.ec_packetsize = ec_packetsize;
                    ec_decode_arg.ec_blocksize = ec_blocksize;
                    ec_decode_arg.ec_matrix = ec_matrix;
                    ec_decode_arg.ec_bitmatrix = ec_bitmatrix;
                    ec_decode_arg.erasures = erasures;

                    pool_init(omp_thread_num);
                       
                    for (j = 0; j < num_iteration; j++)
                    {
                        int target = j%omp_thread_num;
                        //fprintf(stdout, "add worker\n");
                        pool_add_worker(ec_decode_thread, &target);
                    }

                    while(1){
                        pthread_mutex_lock(&lockOfDecodeSum);
                        if(decode_sum == num_iteration){
                            pthread_mutex_unlock(&lockOfDecodeSum);
                            break;
                        }
                        pthread_mutex_unlock(&lockOfDecodeSum);
                    }
                    pool_destroy();
                    /****************************thread_pool*******************************/

                    /****************************omp*******************************/
                    // int omp_id=0;
                    // i=0;
                    // double omp_startTime;
                    // double omp_endTime;
                    // int num_total=0;
                    // //i = jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize, ec_packetsize, 1);
                    // fprintf(stdout,"ec_method=%d\n", method);
                    // if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
                    // {
                    //     #pragma omp parallel for reduction(+:i) num_threads(omp_thread_num) private(omp_id) 
                    //     for (j = 0; j < num_iteration; j++){
                    //         omp_id = omp_get_thread_num();
                    //         //fprintf(stdout,"ompid = %d\n", omp_id);
                    //         i = jerasure_matrix_decode(k, m, w, ec_matrix, 1, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize);
                    //     }
                            
                    // }
                    // else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
                    // {
                    //     #pragma omp parallel for num_threads(omp_thread_num) //private(omp_id,omp_startTime,omp_endTime)
                    //     for (j = 0; j < num_iteration; j++){
                    //         //omp_id = omp_get_thread_num();
                    //         //fprintf(stdout,"ompid = %d\n", omp_id);
                    //         //omp_startTime = GetTimeStamp();
                    //         jerasure_schedule_decode_lazy(k, m, w, ec_bitmatrix, erasures, omp_data[omp_id], omp_coding[omp_id], ec_blocksize, ec_packetsize, 1);
                    //         //omp_endTime = GetTimeStamp();
                    //         num_total++;
                    //         //fprintf(stdout, "time = %lf\n", omp_endTime-omp_startTime);
                    //     }          
                    // }
                    /****************************omp*******************************/

                    if (decode_res != 0)
                    {
                        pthread_mutex_unlock(&lockOfNT);
                        ERR("decode failed!");
                    }
                    else
                    {
                        if (test->ec_verbose >= VERBOSE_0)
                        {
                            fprintf(stdout, "decode success for %d times!\n", decode_sum);
                        }
                    }
                    pthread_mutex_unlock(&lockOfNT);
                    break;
                }
                pthread_mutex_unlock(&lockOfNT);
            }
        }

        ec_strategy_endTime = GetTimeStamp();

        /*************************inter-process repair******************************/
        // if(test->ec_strategy == COLLECTIVE_THREAD){
        //     while(1){
        //         //fprintf(stdout, "in while1\n");
        //         MPI_Allreduce(&local_finished,&global_finished,1,MPI_INT,MPI_SUM,testComm);
        //         fprintf(stdout, "process %d, current global finished: %d\n", rank, global_finished);
        //         if(global_finished == (numTasksWorld-1)){
        //             break;
        //         }               
        //         sleep(1);
        //     }

        //     if(rank == 0 && global_finished == (numTasksWorld-1)){
        //         fprintf(stdout,"process %d, only lack one file!!!!!!!!!!!!!!!!!1\n", rank, global_finished);
        //     }

        //     int cancel;
        //     if(local_finished == 0){
        //         for(i=0;i<total_stripe_num;i++){
        //             cancel = pthread_cancel(threads[i]);
        //             if(cancel == 0){
        //                 fprintf(stdout,"process %d: cancel thread %d success!\n", rank, i);
            
        //             }
        //         }
        //     }else if(rank == (numTasksWorld-1)){
        //         //prepare matrix
        //         int k1 = numTasksWorld-1;
        //         int m1 = 1;
        //         int erasures[2] = {0,-1};
        //         switch (method)
        //         {
        //             case No_Coding:
        //                 break;
        //             case Reed_Sol_Van:
        //                 ec_matrix = reed_sol_vandermonde_coding_matrix(k1, m1, w);
        //                 break;
        //             case Reed_Sol_R6_Op:
        //                 ec_matrix = reed_sol_r6_coding_matrix(k1, w);
        //                 break;
        //             case Cauchy_Orig:
        //                 ec_matrix = cauchy_original_coding_matrix(k1, m1, w);
        //                 ec_bitmatrix = jerasure_matrix_to_bitmatrix(k1, m1, w, ec_matrix);
        //                 break;
        //             case Cauchy_Good:
        //                 ec_matrix = cauchy_good_general_coding_matrix(k1, m1, w);
        //                 ec_bitmatrix = jerasure_matrix_to_bitmatrix(k1, m1, w, ec_matrix);
        //                 break;
        //             case Liberation:
        //                 ec_bitmatrix = liberation_coding_bitmatrix(k1, w);
        //                 break;
        //             case Blaum_Roth:
        //                 ec_bitmatrix = blaum_roth_coding_bitmatrix(k1, w);
        //                 break;
        //             case Liber8tion:
        //                 ec_bitmatrix = liber8tion_coding_bitmatrix(k1);
        //         }
                
        //         //reveive other data
        //         int j;
        //         MPI_Status status;
        //         char **mpi_data = (char**)malloc(sizeof(char*)*k1);
        //         for(i=0;i<k1;i++){
        //             mpi_data[i] = (char*)malloc(sizeof(char)*ec_blocksize);
        //             if(mpi_data[i]!=NULL){
        //                 fprintf(stdout,"allocate mpi_data success!\n");
        //             }
        //         }
        //         char **mpi_coding = (char**)malloc(sizeof(char*)*m1);
        //         for(i=0;i<k;i++){
        //             for(j=0;j<m1;j++){
        //                 mpi_coding[j] = ec_data[i]; //coding是自己的data
        //             }
        //             for(j=0;j<(k1-1);j++){ //data从别的进程来,先忽略0号进程
        //                 MPI_Recv(mpi_data[j+1],ec_blocksize,MPI_CHAR,j+1,j+1,testComm, &status);
        //             }

        //             if (method == Reed_Sol_Van || method == Reed_Sol_R6_Op)
        //             {
        //                 for (j = 0; j < 1500; j++)
        //                 {
        //                     //decode_startTime = GetTimeStamp();
        //                     jerasure_matrix_decode(k1, m1, w, ec_matrix, 1, erasures, mpi_data, mpi_coding, ec_blocksize);
        //                     //decode_endTime = GetTimeStamp();
        //                     //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
        //                 }
        //             }
        //             else if (method == Cauchy_Orig || method == Cauchy_Good || method == Liberation || method == Blaum_Roth || method == Liber8tion)
        //             {
        //                 for (j = 0; j < 1500; j++)
        //                 {
        //                     //decode_startTime = GetTimeStamp();
        //                     jerasure_schedule_decode_lazy(k1, m1, w, ec_bitmatrix, erasures, mpi_data, mpi_coding, ec_blocksize, ec_packetsize, 1);
        //                     //decode_endTime = GetTimeStamp();
        //                     //fprintf(stdout,"time = %lf\n",decode_endTime-decode_startTime);
        //                 }
        //             }
        //             fprintf(stdout, "process %d: decode %d times!!\n", rank,j);

        //         }

        //         fprintf(stdout, "process %d: I/O complete!!\n", rank);
                
        //     }else{
        //         for(i=0;i<k;i++){
        //             MPI_Send(ec_data[i],ec_blocksize,MPI_CHAR,numTasksWorld-1,rank,testComm);
        //         }
        //         fprintf(stdout, "process %d: finish sending!!\n", rank);
        //     }

        //     /*************************inter-process repair******************************/
        //     if(local_finished!=0){
        //         for (i = 0; i < total_stripe_num; i++)
        //         {
        //             pthread_join(threads[i], NULL);
        //         }
        //     }
            
        //     /*ec when k stripes arrive*/
        //     fprintf(stdout, "process %d: break after join..\n", rank);
        // }else{
        //     for (i = 0; i < total_stripe_num; i++)
        //     {
        //         pthread_join(threads[i], NULL);
        //     }
        // }

        for (i = 0; i < 1; i++)
        {
            pthread_join(threads[i], NULL);
        }
        
        amtXferred = ec_blocksize * k;
        /*************************ec multi-thread read*************************/

        // amtXferred = IOR_Xfer(access, fd, buffer, transfer, test);
        // if (amtXferred != transfer)
        //     ERR("cannot read from file");
        ec_endTime = GetTimeStamp();
    }


    while ((offsetArray[pairCnt] != -1) && !hitStonewall) {
    
        test->offset = offsetArray[pairCnt];
        test->offset = test->offset/test->ec_k;
        //ec_offset
        /*
         * fills each transfer with a unique pattern
         * containing the offset into the file
         */
        if (test->storeFileOffset == TRUE)
        {
            FillBuffer(buffer, test, test->offset, pretendRank);
        }
        transfer = test->transferSize;
        if (access == WRITE)
        {
            /*ec transfer*/
            // for(i=0;i<total_stripe_num;i++){
            //     ec_amtXferred[i] = IOR_Xfer(access, ec_fds[i], ec_data[0], ec_blocksize, test);
            // }
            for(i=0;i<k;i++){
                ec_amtXferred[i] = IOR_Xfer(access, ec_fds[i], ec_data[i], ec_blocksize, test);
            }
            for(i=0;i<m;i++){
                ec_amtXferred[k+i] = IOR_Xfer(access, ec_fds[k+i], ec_coding[i], ec_blocksize, test);
            }
            /*ec transfer*/
            //amtXferred = IOR_Xfer(access, fd, buffer, transfer, test); //origin_mark
            amtXferred = ec_blocksize*k;
            for (i = 0; i < total_stripe_num; i++)
            {
                if (ec_amtXferred[i] != ec_blocksize)
                    ERR("write to file failed!");
            }
            
            //fprintf(stdout, "break at 3102\n");
        }
        else if (access == READ)
        {
            amtXferred = test->ec_stripe_size * test->ec_k;
        }
        else if (access == WRITECHECK)
        {
            // memset(checkBuffer, 'a', transfer);
            // amtXferred = IOR_Xfer(access, fd, checkBuffer, transfer, test);
            // if (amtXferred != transfer)
            //     ERR("cannot read from file write check");
            // transferCount++;
            // errors += CompareBuffers(buffer, checkBuffer, transfer,
            //                          transferCount, test, WRITECHECK);
        }
        else if (access == READCHECK)
        {
            fprintf(stdout, "in READCHECK");
            // ReadCheck(fd, buffer, checkBuffer, readCheckBuffer, test,
            //           transfer, test->blockSize, &amtXferred,
            //           &transferCount, access, &errors);
        }
        dataMoved += amtXferred;
        pairCnt++;

        hitStonewall = ((test->deadlineForStonewalling != 0) && ((GetTimeStamp() - startForStonewall) > test->deadlineForStonewalling));
    
    }

print_info:
    /*print ec time info*/

    if(access == READ){
        for (i = 0; i < total_stripe_num; i++)
        {
            fprintf(stdout, "#read time of stripe %d : %lf\n", i, ec_timers[i]);
        }
        fprintf(stdout, "#Total read time of %d stripes: %lf\n", total_stripe_num, ec_endTime - ec_startTime);
        fprintf(stdout, "#EC strategy time of %d : %lf\n", test->ec_strategy, ec_strategy_endTime - ec_strategy_startTime);
    }
    


    totalErrorCount += CountErrors(test, access, errors);

    FreeBuffers(access, checkBuffer, readCheckBuffer, buffer, offsetArray);

    if (access == WRITE && test->fsync == TRUE)
    {
        //IOR_Fsync(fd, test); /*fsync after all accesses*/
        for(i=0;i<total_stripe_num;i++){
            IOR_Fsync(ec_fds[i], test);
        }
    }
    free(ec_data);
    free(ec_coding);
    
    return (dataMoved);
} /* WriteOrRead_ec() */

/****************************ec functions*******************************************/

/******************************************************************************/
/*
 * Write times taken during each iteration of the test.
 */
     
void
WriteTimes(IOR_param_t  * test,
           double      ** timer,
           int            iteration,
           int            writeOrRead)
{
    char accessType[MAX_STR],
         timerName[MAX_STR];
    int  i,
         start,
         stop;

    if (writeOrRead == WRITE) {
        start = 0;
        stop = 6;
        strcpy(accessType, "WRITE");
    } else if (writeOrRead == READ) {
        start = 6;
        stop = 12;
        strcpy(accessType, "READ");
    } else {
        ERR("incorrect WRITE/READ option");
    }

    for (i = start; i < stop; i++) {
        switch(i) {
        case 0: strcpy(timerName, "write open start"); break;
        case 1: strcpy(timerName, "write open stop"); break;
        case 2: strcpy(timerName, "write start"); break;
        case 3: strcpy(timerName, "write stop"); break;
        case 4: strcpy(timerName, "write close start"); break;
        case 5: strcpy(timerName, "write close stop"); break;
        case 6: strcpy(timerName, "read open start"); break;
        case 7: strcpy(timerName, "read open stop"); break;
        case 8: strcpy(timerName, "read start"); break;
        case 9: strcpy(timerName, "read stop"); break;
        case 10: strcpy(timerName, "read close start"); break;
        case 11: strcpy(timerName, "read close stop"); break;
        default: strcpy(timerName, "invalid timer"); break;
        }
        fprintf(stdout, "Test %d: Iter=%d, Task=%d, Time=%f, %s\n",
                test->id, iteration, (int)rank, timer[i][iteration], timerName);
    }
} /* WriteTimes() */
