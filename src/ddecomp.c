#include "ddecomp.h"

#include <mpi.h>

#include "alloc.h"
#include "ftype.h"
#include "field.h"

void allgather_diag_rd(ftype *dst,
                       ftype d0,
                       ftype dm)
{
    ftype buff[2] = { d0, dm };
    MPI_Allgather(buff, 2, MPI_FTYPE, dst, 2, MPI_FTYPE, MPI_COMM_WORLD);
}

void tdma_solve(const ftype *restrict b,
                ftype *restrict a,
                ftype *restrict c,
                ftype *restrict f, /* f gets overwritten by u */
                uint32_t n)
{
    f[0] /= b[0];
    c[0] /= b[0];
    /* Forward sweep. */
    for (uint32_t i = 1; i < n; ++i) {
        ftype r = 1.0 / (b[i] - c[i - 1] * a[i]);
        f[i] = r * (f[i] - f[i - 1] * a[i]);
        c[i] *= r;
    }

    /* Backward sweep. */
    for (uint32_t i = 1; i < n; ++i) {
        f[n - 1 - i] -= c[n - 1 - i] * f[n - i];
    }
}

void tdma_mod_reduce(const ftype *restrict b,
                     ftype *restrict a,
                     ftype *restrict c,
                     ftype *restrict f,
                     uint32_t n)
{
    #pragma GCC unroll(2)
    for (uint32_t i = 0; i < 2; ++i) {
        f[i] /= b[i];
        a[i] /= b[i];
        c[i] /= b[i];
    }
    /* Forward sweep. */
    for (uint32_t i = 2; i < n; ++i) {
        ftype r = 1.0 / (b[i] - c[i - 1] * a[i]);
        f[i] = r * (f[i] - f[i - 1] * a[i]);
        a[i] *= r * -a[i - 1];
        c[i] *= r;
    }

    /* Backward sweep. */
    for (uint32_t i = 2; i < n - 1; ++i) {
        f[n - 1 - i] -= c[n - 1 - i] * f[n - i];
        a[n - 1 - i] -= c[n - 1 - i] * a[n - i];
        c[n - 1 - i] *= -c[n - i];
    }
    ftype r = 1.0 / (1.0 - c[0] * a[1]);
    f[0] = r * (f[0] - c[0] * f[1]);
    c[0] *= r * -c[1];
    a[0] *= r;
}

void tdma_mod_substitute(const ftype *restrict a,
                         const ftype *restrict c,
                         ftype *restrict f,
                         ftype u0,
                         ftype um, /* m=n-1 */
                         uint32_t n)
{
    f[0] = u0;
    for (uint32_t i = 1; i < n - 1; ++i) {
        f[i] -= (u0 * a[i] + um * c[i]);
    }
    f[n - 1] = um;
}

DDecomp *ddecomp_create(uint32_t global_depth,
                        uint32_t global_height,
                        uint32_t global_width,
                        ArenaAllocator *arena)
{
    DDecomp *ddecomp = arena_push_noalign(arena, sizeof(DDecomp));
    ddecomp->comm_z = MPI_COMM_WORLD;

    field_size size = { global_width, global_height, global_depth };
    ddecomp->global_size = size;
    /* WARNING: Assumes global_depth is multiple of num_procs. */
    size.depth = size.depth / get_num_procs(ddecomp->comm_z);
    ddecomp->local_size = size;

    return ddecomp;
}
