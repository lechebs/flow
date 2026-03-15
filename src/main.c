#include <stdio.h>

#include "alloc.h"
#include "boundary.h"
#include "consts.h"
#include "solver.h"
#include "timeit.h"
#include "field.h"
#include "output.h"
#include "thread-array.h"

#define DEPTH 512
#define HEIGHT 512
#define WIDTH 512

#define NUM_TIMESTEPS 20

#define NUM_THREADS 32

#define OUTPUT_PATH "/g100_work/pMI25_DAER/cavity/solution.vtk"

DEFINE_NU(1.0)
DEFINE_DT(0.025)
DEFINE_DX(1.0 / WIDTH)

DEFINE_CONSTANT_FORCING(0, 0, 0)

DEFINE_CONSTANT_BC_U(0, 0, 0, BC_LEFT)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_RIGHT)
DEFINE_CONSTANT_BC_U(0, 0.1, 0, BC_TOP)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BOTTOM)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_FRONT)
DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BACK)

typedef struct {
    Solver *solver;
    OutputVTK *output;
} SimulationData;

static void *run_simulation(void *t_data)
{
    SimulationData *sim_data = thread_get_shared_data(t_data);

    Solver *solver = sim_data->solver;
    OutputVTK *output = sim_data->output;

    solver_init(solver, t_data);
    thread_wait_on_barrier(t_data);

    uint32_t t_id = ((Thread *) t_data)->t_id;
    //output_vtk_write(output, "output/solution-cavity-0.vtk", t_data);
    //thread_wait_on_barrier(t_data);

    TIMER_CREATE(solver_total);
    TIMER_RESTART(solver_total);

    TIMER_CREATE(solver_step_aggregate);
    for (int t = 1; t < NUM_TIMESTEPS + 1; ++t) {
        TIMER_RESTART(solver_step_aggregate);
        solver_step(solver, t, t_data);

        //char output_file_name[64];
        //sprintf(output_file_name, "output/solution-cavity-%d.vtk", t);
        //TIMEITN(output_vtk_write(output, output_file_name, t_data), 1);

        TIMER_ELAPSED(solver_step_aggregate, t_id == 0);
        if (t_id == 0) { printf("\n"); }
    }

    //TIMEITN(output_vtk_write(output, OUTPUT_PATH, t_data), 1);

    TIMER_ELAPSED(solver_total, t_id == 0);
    /* Just in case the compiler optimizes away the whole program. */
    if (t_id == 0 && solver_get_velocity(solver).x[t_id] +
		     solver_get_velocity(solver).y[t_id] +
		     solver_get_velocity(solver).z[t_id] != 1e15) {

	printf("done\n");
    }

    return 0;
}

int main(void)
{
    ArenaAllocator arena;
    arena_init(&arena, 1lu << 35);

    Solver *solver = solver_alloc(DEPTH, HEIGHT, WIDTH, &arena);

    field_size size = { WIDTH, HEIGHT, DEPTH };
    OutputVTK *output = output_vtk_create(size, _DX, &arena);

    output_vtk_attach_field(output, solver_get_porosity(solver),
                            "porosity", &arena);
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

    return 0;
}
