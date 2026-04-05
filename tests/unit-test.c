#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "test.h"
#include "ftype.h"
#include "field.h"
#include "consts.h"
#include "ddecomp.h"

DEFINE_NU(1.00)
DEFINE_DT(0.01)
DEFINE_DX(0.01)

DEF_TEST(test_vtranspose)
{
#ifdef FLOAT
    float __attribute__((aligned(32))) m[64];
    for (int i = 0; i < 64; ++i) {
        m[i] = ((float) rand()) / RAND_MAX;
    }

    __m256 r0 = _mm256_load_ps(m);
    __m256 r1 = _mm256_load_ps(m + 8);
    __m256 r2 = _mm256_load_ps(m + 16);
    __m256 r3 = _mm256_load_ps(m + 24);
    __m256 r4 = _mm256_load_ps(m + 32);
    __m256 r5 = _mm256_load_ps(m + 40);
    __m256 r6 = _mm256_load_ps(m + 48);
    __m256 r7 = _mm256_load_ps(m + 56);

    vtranspose(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);

    float __attribute__((aligned(32))) t[64];

    _mm256_store_ps(t, r0);
    _mm256_store_ps(t + 8, r1);
    _mm256_store_ps(t + 16, r2);
    _mm256_store_ps(t + 24, r3);
    _mm256_store_ps(t + 32, r4);
    _mm256_store_ps(t + 40, r5);
    _mm256_store_ps(t + 48, r6);
    _mm256_store_ps(t + 56, r7);

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            EXPECT_EQUAL(m[i * 8 + j], t[j * 8 + i]);
        }
    }

#else
    double __attribute__((aligned(32))) m[16];
    for (int i = 0; i < 16; ++i) {
        m[i] = ((double) rand()) / RAND_MAX;
    }

    __m256d r0 = _mm256_load_pd(m);
    __m256d r1 = _mm256_load_pd(m + 4);
    __m256d r2 = _mm256_load_pd(m + 8);
    __m256d r3 = _mm256_load_pd(m + 12);

    vtranspose(&r0, &r1, &r2, &r3);

    double __attribute__((aligned(32))) t[16];

    _mm256_store_pd(t, r0);
    _mm256_store_pd(t + 4, r1);
    _mm256_store_pd(t + 8, r2);
    _mm256_store_pd(t + 12, r3);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_EQUAL(m[i * 4 + j], t[j * 4 + i]);
        }
    }
#endif
}

DEF_TEST(test_tdma_mod, ArenaAllocator *arena);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    srand(42 + get_proc_rank(MPI_COMM_WORLD));

    ArenaAllocator arena;
    arena_init(&arena, 1ul << 20);

    if (get_proc_rank(MPI_COMM_WORLD) == 0) {
        RUN_TEST(test_vtranspose);
    }

    RUN_TEST(test_tdma_mod, &arena);

    MPI_Finalize();

    return 0;
}
