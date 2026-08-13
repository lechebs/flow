#ifndef TIMEIT_H
#define TIMEIT_H

#include <stdio.h>
#include <string.h>
#include <time.h>

#define TIMEIT(func_call) TIMEITN(func_call, 10)

#ifdef TIMER

#define TIMEITN(func_call, avg_iter)                              \
do {                                                              \
    long elapsed_ns_avg = 0;                                      \
    /* func_call;*/ /* Warmup run. */                             \
    struct timespec start, stop;                                  \
    clock_gettime(CLOCK_MONOTONIC, &start);                       \
    for (int i = 0; i < avg_iter; ++i) {                          \
        func_call;                                                \
    }                                                             \
    clock_gettime(CLOCK_MONOTONIC, &stop);                        \
    elapsed_ns_avg = (stop.tv_sec - start.tv_sec) * 1e9 +         \
                     (stop.tv_nsec - start.tv_nsec);              \
    char *func_name = #func_call;                                 \
    char func_name_buff[32] = {'\0'};                             \
    strncpy(func_name_buff, func_name, 31);                       \
    char *ptr = strstr(func_name_buff, "(");                      \
    memset(ptr, ' ', 32 - (ptr - func_name_buff));                \
    func_name_buff[31] = '\0';                                    \
    printf("%s%8.2f ms [%d runs avg]\n",                          \
           func_name_buff,                                        \
           elapsed_ns_avg / (1e6 * avg_iter), avg_iter);          \
} while (0)

#define TIMER_CREATE(name)                      \
    struct timespec _timer_##name##_start,      \
                    _timer_##name##_curr;       \
    static long _timer_##name##_tot_elapsed_ns; \
    static long _timer_##name##_elapsed_count;  \

#define TIMER_RESTART(name)                                  \
do {                                                         \
    clock_gettime(CLOCK_MONOTONIC, &_timer_##name##_start);  \
} while (0)

#define TIMER_ELAPSED(name, cond)                                  \
do {                                                               \
    if (cond) {                                                    \
        clock_gettime(CLOCK_MONOTONIC, &_timer_##name##_curr);     \
        long elapsed_ns = (_timer_##name##_curr.tv_sec -           \
                           _timer_##name##_start.tv_sec) * 1e9 +   \
                          (_timer_##name##_curr.tv_nsec -          \
                           _timer_##name##_start.tv_nsec);         \
        _timer_##name##_tot_elapsed_ns += elapsed_ns;              \
        _timer_##name##_elapsed_count++;                           \
        if (!(_timer_##name##_elapsed_count % TIMER_LOG_FREQ)) {   \
            printf("%-40s%8.2f ms %8.2f ms (avg %ld runs)\n",      \
                   #name, elapsed_ns / 1e6,                        \
                   _timer_##name##_tot_elapsed_ns / 1e6 /          \
                   _timer_##name##_elapsed_count,                  \
                   _timer_##name##_elapsed_count);                 \
        }                                                          \
    }                                                              \
} while (0)

#define TIMER_NEWLINE(cond)     \
do {                            \
    if (cond) { printf("\n"); } \
} while (0)

#else

#define TIMEITN(func_call, ...) func_call;
#define TIMER_CREATE(...) ;
#define TIMER_RESTART(...) ;
#define TIMER_ELAPSED(...) ;
#define TIMER_NEWLINE(...) ;

#endif

#endif
