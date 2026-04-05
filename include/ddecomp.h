#ifndef DDECOMP_H
#define DDECOMP_H

typedef struct {
    MPI_Comm comm;
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

#endif
