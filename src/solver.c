#include <math.h>

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
    field3 velocity_pred;
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
    solver->velocity_pred = field3_alloc_pad(domain_size, arena);

    return solver;
}

void solver_init(Solver *solver, ArenaAllocator *arena)
{
    field_size domain_size = solver->domain_size;

    momentum_init(domain_size, solver->velocity_Dxx);
    momentum_init(domain_size, solver->velocity_Dyy);
    momentum_init(domain_size, solver->velocity_Dzz);
    momentum_init(domain_size, solver->velocity_pred);

    pressure_init(domain_size, solver->pressure);
    pressure_init(domain_size, solver->pressure_delta);

    arena_enter(arena);

    /* Setting constant unit porosity. */
    field tmp = field_alloc(domain_size, arena);
#ifdef FLOW_PAST_CUBE
    field_fill(domain_size, 1e8, tmp);

    for (uint32_t i = 0; i < domain_size.depth; ++i) {
        for (uint32_t j = 0; j < domain_size.height; ++j) {
            for (uint32_t k = 0; k < domain_size.width; ++k) {
                uint64_t idx = field_idx(domain_size, k, j, i);

                ftype x = k * _DX;
                ftype y = j * _DX;
                ftype z = i * _DX;

                double sdf = pow(pow(y - _DX * domain_size.height * 0.5, 8) +
                                 pow(x - _DX * domain_size.width * 0.25, 8) +
                                 pow(z - _DX * domain_size.depth * 0.5, 8), 1.0 / 8) - 0.09;

                double exp = 8 * tanh(sdf * 30);
                tmp[idx] = pow(10, exp);
            }
        }
    }
#else
    field_fill(domain_size, 1e20, tmp);
#endif

    solver_set_porosity(solver, tmp);

    arena_exit(arena);
}

void solver_set_porosity(Solver *solver, const ftype *src)
{
    field_copy(solver->domain_size, src, solver->porosity);
    compute_gamma(src, solver->domain_size, solver->gamma);
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
                   solver->velocity_pred,
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
