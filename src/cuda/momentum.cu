#include "ftype.h"
#include "field.h"

/* cuda/ includes */
#include "utils.h"
#include "momentum.h"

#define CUDA_H2D cudaMemcpyHostToDevice
#define CUDA_D2H cudaMemcpyDeviceToHost

struct SolveDzzParams
{
    ftype *__restrict__ w;
    ftype *__restrict__ tmp;
    ftype *__restrict__ f_x;
    ftype *__restrict__ f_y;
    ftype *__restrict__ f_z;
    ftype *__restrict__ u_x;
    ftype *__restrict__ u_y;
    ftype *__restrict__ u_z;
};

static struct SolveDzzParams d_data;

void alloc_device_data(field_size domain_size)
{
    uint64_t alloc_size = field_num_points(domain_size) * sizeof(ftype);

    /* WARNING: Some of these are padded in the cpu version. */
    CUDA_CHECK(cudaMalloc((void **) &d_data.w, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.tmp, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_z, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_z, alloc_size));
}

void free_device_data()
{
    CUDA_CHECK(cudaFree(d_data.w));
    CUDA_CHECK(cudaFree(d_data.tmp));
    CUDA_CHECK(cudaFree(d_data.f_x));
    CUDA_CHECK(cudaFree(d_data.f_y));
    CUDA_CHECK(cudaFree(d_data.f_z));
    CUDA_CHECK(cudaFree(d_data.u_x));
    CUDA_CHECK(cudaFree(d_data.u_y));
    CUDA_CHECK(cudaFree(d_data.u_z));
}

__global__ void solve_Dzz_sweep(__grid_constant__ const SolveDzzParams params)
{
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;

    params.u_x[idx] = (ftype) 1.0;
    params.u_y[idx] = (ftype) 1.0;
    params.u_z[idx] = (ftype) 1.0;
}

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
                               ftype *__restrict__ u_z)
{
    uint64_t alloc_size = depth * height * width * sizeof(ftype);

    // NOTE: The first transfer is redundant
    CUDA_CHECK(cudaMemcpy(d_data.w, w, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_x, f_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_y, f_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_z, f_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_x, u_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_y, u_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_z, u_z, alloc_size, CUDA_H2D));

    CUDA_CHECK_LAUNCH(solve_Dzz_sweep<<<1, 128>>>(d_data));

    CUDA_CHECK(cudaMemcpy(u_x, d_data.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data.u_z, alloc_size, CUDA_D2H));

    // No need to synchronize
}
