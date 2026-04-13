#ifndef DDECOMP_H
#define DDECOMP_H

#include <mpi.h>

#include "alloc.h"
#include "field.h"

#ifdef FLOAT
#define MPI_FTYPE MPI_FLOAT
#else
#define MPI_FTYPE MPI_DOUBLE
#endif

enum DDecompType { DDECOMP_SLAB_Z, DDECOMP_PENCIL_XY };

extern const enum DDecompType _ddecomp_type;
/* Could be useful to fold branches for boundary code. */
#define DEFINE_DDECOMP_TYPE(type) \
    const enum DDecompType _ddecomp_type = type;

typedef struct {
    MPI_Comm comm_z;

    int num_procs;
    int proc_rank;

    field_size local_size;
    field_size global_size;

} DDecomp;

static inline int get_num_procs(MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

static inline int get_proc_rank(MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

void allgather_diag_rd(ftype *dst,
                       ftype d0,
                       ftype dm);

void tdma_solve(const ftype *restrict b,
                ftype *restrict a,
                ftype *restrict c,
                ftype *restrict f, /* f gets overwritten by u */
                uint32_t n);

void tdma_mod_reduce(const ftype *restrict b,
                     ftype *restrict a,
                     ftype *restrict c,
                     ftype *restrict f,
                     uint32_t n);

void tdma_mod_substitute(const ftype *restrict a,
                         const ftype *restrict c,
                         ftype *restrict f,
                         ftype u0,
                         ftype um, /* m=n-1 */
                         uint32_t n);

DDecomp *ddecomp_create(uint32_t global_depth,
                        uint32_t global_height,
                        uint32_t global_width,
                        ArenaAllocator *arena);

#endif
