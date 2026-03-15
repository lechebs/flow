#include "solver.h"
#include "momentum.h"
#include "pressure.h"
#include "ftype.h"
#include "field.h"
#include "consts.h"
#include "timeit.h"
#include "thread-array.h"

struct Solver {
    field_size domain_size;
    field porosity;
    field gamma;
    field pressure;
    field pressure_delta;
    field3 velocity_Dxx;
    field3 velocity_Dyy;
    field3 velocity_Dzz;
};

Solver *solver_alloc(uint32_t domain_depth,
                     uint32_t domain_height,
                     uint32_t domain_width,
                     ArenaAllocator *arena)
{
    Solver *solver = arena_push(arena, sizeof(Solver));
    memset(solver, 0, sizeof(Solver));

    field_size domain_size;
    domain_size.depth = domain_depth;
    domain_size.height = domain_height;
    domain_size.width = domain_width;
    solver->domain_size = domain_size;

    solver->porosity = field_alloc(domain_size, arena);
    solver->gamma = field_alloc(domain_size, arena);
    solver->pressure = field_alloc(domain_size, arena);
    solver->pressure_delta = field_alloc(domain_size, arena);

    solver->velocity_Dxx = field3_alloc_pad(domain_size, arena);
    solver->velocity_Dyy = field3_alloc_pad(domain_size, arena);
    solver->velocity_Dzz = field3_alloc_pad(domain_size, arena);

    return solver;
}

void solver_init(Solver *solver, Thread *thread)
{
    field_size domain_size = solver->domain_size;

    momentum_init(domain_size, solver->velocity_Dxx, thread);
    momentum_init(domain_size, solver->velocity_Dyy, thread);
    momentum_init(domain_size, solver->velocity_Dzz, thread);

    pressure_init(domain_size, solver->pressure, thread);
    pressure_init(domain_size, solver->pressure_delta, thread);

    ArenaAllocator *arena = thread_get_arena(thread);
    arena_enter(arena);

    field_size subdomain_size = domain_size;
    subdomain_size.depth /= thread_get_array_size(thread);
    /* Setting constant unit porosity. */
    field tmp = field_alloc(subdomain_size, arena);
    field_fill(subdomain_size, 1e20, tmp);

    /*
    for (uint32_t i = 0; i < domain_size.depth; ++i) {
        for (uint32_t j = 0; j < domain_size.height; ++j) {
            for (uint32_t k = 0; k < domain_size.width; ++k) {
                uint64_t idx = field_idx(domain_size, k, j, i);

                if (((i < 48 && i > 16) && (j > 16 && j < 26)) ||
                    ((i > 16) && (j > 46 && j < 56) && (k < 32)) ||
                    ((i < 48) && (j > 76 && j < 86) && (k > 32)) ||
                    ((i > 32) && (j > 106 && j < 116))
                ) {
                    tmp[idx] = 1e-20;
                }
            }
        }
    }
    */

    solver_set_porosity(solver, tmp, thread);

    thread_wait_on_barrier(thread);

    arena_exit(arena);
}

void solver_set_porosity(Solver *solver, const ftype *src, Thread *thread)
{
    uint32_t num_threads = thread_get_array_size(thread);
    field_copy_numalocal(solver->domain_size, src, solver->porosity,
		         thread->t_id, num_threads);
    compute_gamma(src, solver->domain_size, solver->gamma,
		  thread->t_id, num_threads);
}

void solver_step(Solver *solver, uint32_t timestep, Thread *thread)
{
    momentum_solve(solver->porosity,
                   solver->gamma,
                   solver->pressure,
                   solver->pressure_delta,
                   solver->domain_size,
                   solver->velocity_Dxx,
                   solver->velocity_Dyy,
                   solver->velocity_Dzz,
                   timestep,
                   thread);

    pressure_solve(to_const_field3(solver->velocity_Dzz),
                   solver->domain_size,
                   solver->pressure,
                   solver->pressure_delta,
                   timestep,
                   thread);
}

const_field3 solver_get_velocity(Solver *solver)
{
    return to_const_field3(solver->velocity_Dzz);
}

const_field solver_get_pressure(Solver *solver)
{
    return solver->pressure;
}

const_field solver_get_porosity(Solver *solver)
{
    return solver->porosity;
}
