#ifndef CUDA_MOMENTUM_H
#define CUDA_MOMENTUM_H

#include "ftype.h"

#ifdef __cplusplus
extern "C" {
#endif
/* TODO: Move these somewhere else. */
void alloc_device_data(field_size domain_size);
void free_device_data();

void update_device_consts();

void launch_momentum_solve_Dyy(const ftype *__restrict__ w,
                               uint32_t depth,
                               uint32_t height,
                               uint32_t width,
                               uint32_t timestep,
                               ftype *__restrict__ tmp,
                               ftype *__restrict__ f_x,
                               ftype *__restrict__ f_y,
                               ftype *__restrict__ f_z,
                               ftype *__restrict__ u_x,
                               ftype *__restrict__ u_y,
                               ftype *__restrict__ u_z);

/* Thank god I chose to keep the interface of these function
 * as simple as possible :D */
void launch_momentum_solve_Dzz(const ftype *__restrict__ w,
                               uint32_t depth,
                               uint32_t height,
                               uint32_t width,
                               uint32_t timestep,
                               ftype *__restrict__ tmp,
                               ftype *__restrict__ f_x,
                               ftype *__restrict__ f_y,
                               ftype *__restrict__ f_z,
                               ftype *__restrict__ u_x,
                               ftype *__restrict__ u_y,
                               ftype *__restrict__ u_z);
#ifdef __cplusplus
}
#endif

#endif
