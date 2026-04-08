#ifndef SOLVER_H
#define SOLVER_H

#include "ftype.h"
#include "field.h"
#include "alloc.h"
#include "thread-array.h"
#include "ddecomp.h"

struct Solver;

typedef struct Solver Solver;

Solver *solver_alloc(DDecomp *ddecomp, ArenaAllocator *arena);

void solver_init(Solver *solver, ArenaAllocator *arena);

void solver_set_porosity(Solver *solver, const ftype *src);

void solver_step(Solver *solver, uint32_t timestep, Thread *thread);

const_field3 solver_get_velocity(Solver *solver);

const_field solver_get_pressure(Solver *solver);

const_field solver_get_porosity(Solver *solver);

#endif
