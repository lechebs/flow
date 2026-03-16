#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "field.h"
#include "alloc.h"
#include "thread-array.h"

void momentum_init(field_size size, field3 field);

void momentum_solve(const_field porosity,
                    const_field gamma,
                    const_field pressure,
                    const_field pressure_delta,
                    field_size size,
                    field3 velocity_Dxx,
                    field3 velocity_Dyy,
                    field3 velocity_Dzz,
                    field3 velocity_Dzz_pred,
                    uint32_t timestep,
                    Thread *thread);

#endif
