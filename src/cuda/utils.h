#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <stdio.h>

#define CUDA_CHECK_LAUNCH(...)                 \
do {                                           \
    __VA_ARGS__;                               \
    CUDA_CHECK(cudaGetLastError());            \
} while (0)

#define CUDA_CHECK(api_call)                   \
do {                                           \
    cudaError_t err = api_call;                \
    if (err != cudaSuccess) {                  \
        printf("%s:%d %s: %s\n",               \
               __FILE__,                       \
               __LINE__,                       \
               cudaGetErrorName(err),          \
               cudaGetErrorString(err));       \
    }                                          \
} while (0)

#endif
