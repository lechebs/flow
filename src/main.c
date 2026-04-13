#include <stdio.h>
#include <mpi.h>

#include "alloc.h"
#include "boundary.h"
#include "consts.h"
#include "solver.h"
#include "timeit.h"
#include "field.h"
#include "output.h"
#include "thread-array.h"

#define DEPTH 64
#define HEIGHT 64
#define WIDTH 64

#define NUM_TIMESTEPS 10

#define NUM_THREADS 1

DEFINE_NU(1.0)
DEFINE_DT(0.025)
DEFINE_DX(M_PI / WIDTH)

DEFINE_CONSTANT_FORCING(0, 0, 0)

DEFINE_CONSTANT_BC_U(0, 0, 0, BC_LEFT)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_RIGHT)
DEFINE_CONSTANT_BC_U(0, 0, 0.1, BC_TOP)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BOTTOM)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_FRONT)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BACK)

typedef struct {
    Solver *solver;
    OutputVTK *output;
} SimulationData;

static void *run_simulation(void *t_data)
{
    ArenaAllocator *arena = thread_get_arena(t_data);
    SimulationData *sim_data = thread_get_shared_data(t_data);

    Solver *solver = sim_data->solver;
    OutputVTK *output = sim_data->output;

    uint32_t t_id = ((Thread *) t_data)->t_id;
    int proc_rank = get_proc_rank(MPI_COMM_WORLD);

    char output_file_name[64];
    sprintf(output_file_name,
            "output/solution-cavity-%d-0.vtk", proc_rank);
    output_vtk_write(output, output_file_name, t_data);
    thread_wait_on_barrier(t_data);

    TIMER_CREATE(solver_step_aggregate);
    TIMER_CREATE(output_vtk_write);

    for (int t = 1; t < NUM_TIMESTEPS + 1; ++t) {
        TIMER_RESTART(solver_step_aggregate);
        solver_step(solver, t, t_data);

        TIMER_RESTART(output_vtk_write);
        sprintf(output_file_name,
                "output/solution-cavity-%d-%d.vtk", proc_rank, t);
        output_vtk_write(output, output_file_name, t_data);
        TIMER_ELAPSED(output_vtk_write, t_id == 0 && proc_rank == 0);

        TIMER_ELAPSED(solver_step_aggregate, t_id == 0 && proc_rank == 0);
        if (t_id == 0 && proc_rank == 0) { printf("\n"); }
    }

    return 0;
}

int main(int argc, char *argv[])
{
    /* TODO: Use MPI_init_thread when using multiple threads. */
    MPI_Init(&argc, &argv);

    ArenaAllocator arena;
    arena_init(&arena, 1lu << 34);

    DDecomp *ddecomp = ddecomp_create(DEPTH, HEIGHT, WIDTH, &arena);

    Solver *solver = solver_alloc(ddecomp, &arena);
    solver_init(solver, &arena);

    OutputVTK *output = output_vtk_create(ddecomp->local_size, _DX, &arena);

    //output_vtk_attach_field(output, solver_get_porosity(solver),
    //                        "porosity", &arena);
    output_vtk_attach_field(output, solver_get_pressure(solver),
                            "pressure", &arena);
    output_vtk_attach_field3(output, solver_get_velocity(solver),
                             "velocity", &arena);

    ThreadArray *t_array = thread_array_create(NUM_THREADS, &arena);
    SimulationData data = { solver, output };
    thread_array_set_shared_data(t_array, &data);

    thread_array_run(t_array, run_simulation, &arena);
    thread_array_destroy(t_array);

    arena_destroy(&arena);

    MPI_Finalize();

    return 0;
}
