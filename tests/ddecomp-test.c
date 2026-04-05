#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "test.h"
#include "alloc.h"
#include "ftype.h"

#include "ddecomp.c"

#ifdef FLOAT
#define MPI_FTYPE MPI_FLOAT
#else
#define MPI_FTYPE MPI_DOUBLE
#endif

#define TOL 1e-8

static ftype rand_uniform()
{
    return ((ftype) rand()) / RAND_MAX;
}

static void allgather_diag_rd(ftype *dst,
                              ftype d0,
                              ftype dm)
{
    ftype buff[2] = { d0, dm };
    MPI_Allgather(buff, 2, MPI_FTYPE, dst, 2, MPI_FTYPE, MPI_COMM_WORLD);
}

DEF_TEST(test_tdma_mod, ArenaAllocator *arena)
{
    arena_enter(arena);

    const int num_procs = get_num_procs(MPI_COMM_WORLD);
    const int proc_rank = get_proc_rank(MPI_COMM_WORLD);

    const uint32_t n = 64;
    const uint32_t n_loc = n / num_procs;
    /* Allocating space for the local system. */
    ftype *a = arena_push_count(arena, ftype, n_loc);
    ftype *b = arena_push_count(arena, ftype, n_loc);
    ftype *c = arena_push_count(arena, ftype, n_loc);
    ftype *f = arena_push_count(arena, ftype, n_loc);
    /* Randomly initializing local system coefficients. */
    for (uint32_t i = 0; i < n_loc; ++i) {
        a[i] = rand_uniform();
        b[i] = rand_uniform();
        c[i] = rand_uniform();
        f[i] = rand_uniform();
    }

    /* Copying initial coefficients to allow validation. */
    ftype *a_cp = arena_push_count(arena, ftype, n_loc);
    ftype *c_cp = arena_push_count(arena, ftype, n_loc);
    ftype *f_cp = arena_push_count(arena, ftype, n_loc);
    memcpy(a_cp, a, n_loc * sizeof(ftype));
    memcpy(c_cp, c, n_loc * sizeof(ftype));
    memcpy(f_cp, f, n_loc * sizeof(ftype));

    uint32_t n_rd = num_procs * 2;
    /* Allocating space for the reduced system. */
    ftype *a_rd = arena_push_count(arena, ftype, n_rd);
    ftype *b_rd = arena_push_count(arena, ftype, n_rd);
    ftype *c_rd = arena_push_count(arena, ftype, n_rd);
    ftype *f_rd = arena_push_count(arena, ftype, n_rd);
    /* Reduced system already has only 1s on the diagonal. */
    for (uint32_t i = 0; i < n_rd; ++i) {
        b_rd[i] = 1.0;
    }

    /* Computing reduced system coefficients. */
    tdma_mod_reduce(b, a, c, f, n_loc);

    /* Assembling reduced system. For a non-trivial
     * implementation, consider packing into a single
     * message all of the coefficients. */
    allgather_diag_rd(a_rd, a[0], a[n_loc - 1]);
    allgather_diag_rd(c_rd, c[0], c[n_loc - 1]);
    allgather_diag_rd(f_rd, f[0], f[n_loc - 1]);

    /* Solving reduced system. */
    tdma_solve(b_rd, a_rd, c_rd, f_rd, n_rd);

    ftype u0 = f_rd[2 * proc_rank];
    ftype um = f_rd[2 * proc_rank + 1];
    /* Substituting edge unknowns into local system. */
    tdma_mod_substitute(a, c, f, u0, um, n_loc);

    /* Getting halo unknowns. */
    ftype ulh, urh;
    if (proc_rank > 0) {
        ulh = f_rd[2 * proc_rank - 1];
    } else {
        ulh = 0;
    }
    if (proc_rank < num_procs - 1) {
        urh = f_rd[2 * proc_rank + 2];
    } else {
        urh = 0;
    }

    /* Verifying solution. */
    ASSERT_EQUALF(a_cp[0] * ulh + b[0] * f[0] + c_cp[0] * f[1], f_cp[0], TOL);
    for (uint32_t i = 1; i < n_loc - 1; ++i) {
        ASSERT_EQUALF(a_cp[i] * f[i - 1] + b[i] * f[i] +
                      c_cp[i] * f[i + 1], f_cp[i], TOL);
    }
    ASSERT_EQUALF(a_cp[n_loc - 1] * f[n_loc - 2] +
                  b[n_loc - 1] * f[n_loc - 1] +
                  c_cp[n_loc - 1] * urh, f_cp[n_loc - 1], TOL);

    arena_exit(arena);
}
