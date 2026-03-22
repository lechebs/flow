#include <stdio.h>
#include <math.h>

#include "alloc.h"
#include "boundary.h"
#include "consts.h"
#include "solver.h"
#include "timeit.h"
#include "field.h"
#include "output.h"
#include "thread-array.h"

#define DEPTH 128
#define HEIGHT 192
#define WIDTH 256

#define FINAL_TIME 15.0

#define NUM_THREADS 8

//DEFINE_NU(0.0005)
DEFINE_NU(0.0003)
DEFINE_DT(0.001)
DEFINE_DX(0.00375)

DEFINE_CONSTANT_FORCING(0, 0, 0)

static ftype bc_zero(ftype x, ftype y, ftype z, ftype t)
{
    return 0;
}

static ftype bc_inlet(ftype x, ftype y, ftype z, ftype t)
{
    return 2 * (1.0 - exp(-t));
}

DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_LEFT)
DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_RIGHT)

DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_TOP)
DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_BOTTOM)

DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_FRONT)
DEFINE_FUNCTION_BC_U(bc_inlet, bc_zero, bc_zero, BC_BACK)

//DEFINE_CONSTANT_BC_U(0, 0, 0, BC_TOP)
//DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BOTTOM)
//DEFINE_CONSTANT_BC_U(0, 0, 0, BC_FRONT)
//DEFINE_CONSTANT_BC_U(0, 0, 0, BC_BACK)

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
    output_vtk_write(output, "output/solution-karman-0.vtk", t_data);
    thread_wait_on_barrier(t_data);

    uint32_t num_timesteps = (FINAL_TIME - _DT / 2) / _DT + 1;

    TIMER_CREATE(solver_step_aggregate);
    for (uint32_t t = 1; t < num_timesteps + 1; ++t) {
        TIMER_RESTART(solver_step_aggregate);
        solver_step(solver, t, t_data);

        if (t % 250 == 0) {
            char output_file_name[64];
            sprintf(output_file_name, "output/solution-karman-%d.vtk", t);
            output_vtk_write(output, output_file_name, t_data);
        }

        TIMER_ELAPSED(solver_step_aggregate, t_id == 0);
        TIMER_NEWLINE(t_id == 0 && t % TIMER_LOG_FREQ == 0);
    }

    return 0;
}

int main(void)
{
    ArenaAllocator arena;
    arena_init(&arena, 1lu << 34);

    Solver *solver = solver_alloc(DEPTH, HEIGHT, WIDTH, &arena);
    solver_init(solver, &arena);

    field_size size = { WIDTH, HEIGHT, DEPTH };
    OutputVTK *output = output_vtk_create(size, _DX, &arena);
    uint64_t offset = 0;//HEIGHT * WIDTH * (DEPTH / 2 - 1);

    output_vtk_attach_field(output, solver_get_porosity(solver), offset,
                            "porosity", &arena);
    output_vtk_attach_field(output, solver_get_pressure(solver), offset,
                            "pressure", &arena);
    output_vtk_attach_field3(output, solver_get_velocity(solver), offset,
                             "velocity", &arena);

    ThreadArray *t_array = thread_array_create(0, &arena);

    printf("num_threads=%d\n\n", t_array->num_threads);

    SimulationData data = { solver, output };
    thread_array_set_shared_data(t_array, &data);

    thread_array_run(t_array, run_simulation, &arena);
    thread_array_destroy(t_array);

    arena_destroy(&arena);

    return 0;
}
