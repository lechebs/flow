#include "ftype.h"
#include "field.h"
#include "consts.h"

/* cuda/ includes */
#include "utils.h"
#include "boundary.h"
#include "momentum.h"

// TODO: Consider templating kernel on ftype

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

    uint32_t depth, height, width;
    uint32_t timestep;
};

static struct SolveDzzParams d_data;

void alloc_device_data(field_size domain_size)
{
    uint64_t alloc_size = field_num_points(domain_size) * sizeof(ftype);

    // WARNING: Some of these are padded in the cpu version.
    CUDA_CHECK(cudaMalloc((void **) &d_data.w, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.f_z, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data.u_z, alloc_size));

    field_size tmp_size = { domain_size.width,
                            domain_size.height,
                            domain_size.depth * 4 + 3 };
    uint64_t tmp_alloc_size = field_num_points(tmp_size) * sizeof(ftype);

    CUDA_CHECK(cudaMalloc((void **) &d_data.tmp, tmp_alloc_size));

    d_data.depth = domain_size.depth;
    d_data.height = domain_size.height;
    d_data.width = domain_size.width;
}

void update_device_consts()
{
    CUDA_CHECK(cudaMemcpyToSymbol(_D_DX, &_DX, sizeof(ftype), 0, CUDA_H2D));
    CUDA_CHECK(cudaMemcpyToSymbol(_D_DT, &_DT, sizeof(ftype), 0, CUDA_H2D));
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

template<typename BoundaryConditions>
__global__ void solve_Dzz_sweep(__grid_constant__ const SolveDzzParams params)
{
    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t y = blockIdx.y;

    const uint32_t depth = params.depth;
    const uint32_t height = params.height;
    const uint32_t width = params.width;

    if (x > width - 1) {
        return;
    }

    ftype *__restrict__ tmp = params.tmp;
    // WARNING: Is 256 byte alignement guaranteed?
    ftype *__restrict__ tmp_f_x = tmp + depth * height * width;
    ftype *__restrict__ tmp_f_y = tmp + depth * height * width * 2;
    ftype *__restrict__ tmp_f_z = tmp + depth * height * width * 3;

    ftype *__restrict__ tmp_u_x = tmp + depth * height * width * 4;
    ftype *__restrict__ tmp_u_y = tmp + depth * height * width * 4
                                      + height * width;
    ftype *__restrict__ tmp_u_z = tmp + depth * height * width * 4
                                      + height * width * 2;

    const uint64_t face_stride = height * width;
    const uint64_t xy_offset = width * y + x;

    const BoundaryConditions bc;

    // Apply BCs to the first face
    tmp[xy_offset] = 0; // set upper coefficient to 0
    // Enforce correct solution in rhs
    ftype3 u0 = bc.get_u_delta_front(x, y, 0, params.timestep);
    tmp_f_x[xy_offset] = u0.x;
    tmp_f_y[xy_offset] = u0.y;
    tmp_f_z[xy_offset] = u0.z;

    // Gauss reduce the remaining domain, except last face
    for (uint32_t z = 1; z < depth - 1; ++z) {
        const uint64_t xyz_offset = height * width * z + xy_offset;

        ftype w = params.w[xyz_offset];
        ftype upp_prev = tmp[xyz_offset - face_stride];
        ftype f_x_prev = tmp_f_x[xyz_offset - face_stride];
        ftype f_y_prev = tmp_f_y[xyz_offset - face_stride];
        ftype f_z_prev = tmp_f_z[xyz_offset - face_stride];

        // TODO: Operator overloading?
        ftype f_x = params.f_x[xyz_offset] - params.u_x[xyz_offset];
        ftype f_y = params.f_y[xyz_offset] - params.u_y[xyz_offset];
        ftype f_z = params.f_z[xyz_offset] - params.u_z[xyz_offset];

        ftype norm_coef = w * upp_prev + 1 + 2 * w;
        // TODO: Replace division with product
        tmp[xyz_offset] = -w / norm_coef;
        tmp_f_x[xyz_offset] = (w * f_x_prev + f_x) / norm_coef;
        tmp_f_y[xyz_offset] = (w * f_y_prev + f_y) / norm_coef;
        tmp_f_z[xyz_offset] = (w * f_z_prev + f_z) / norm_coef;
    }

    // Apply bcs to last face, solving directly

    const uint64_t xyz_offset = height * width * (depth - 1) + xy_offset;

    ftype w = params.w[xyz_offset];
    ftype upp_prev = tmp[xyz_offset - face_stride];
    ftype f_x_prev = tmp_f_x[xyz_offset - face_stride];
    ftype f_y_prev = tmp_f_y[xyz_offset - face_stride];
    ftype f_x = params.f_x[xyz_offset] - params.u_x[xyz_offset];
    ftype f_y = params.f_y[xyz_offset] - params.u_y[xyz_offset];

    ftype norm_coef = 1 + 4 * w + (ftype) (4.0 / 3.0) * w * upp_prev;

    ftype3 un = bc.get_u_delta_back(x, y, depth - 1, params.timestep);

    un.x = (w * (ftype) (8.0 / 3.0) * un.x + f_x +
            (ftype) (4.0 / 3.0) * w * f_x_prev) / norm_coef;
    un.y = (w * (ftype) (8.0 / 3.0) * un.y + f_y +
            (ftype) (4.0 / 3.0) * w * f_y_prev) / norm_coef;

    tmp_u_x[xy_offset] = un.x;
    tmp_u_y[xy_offset] = un.y;
    tmp_u_z[xy_offset] = un.z;

    params.u_x[xyz_offset] += un.x;
    params.u_y[xyz_offset] += un.y;
    params.u_z[xyz_offset] += un.z;

    // TODO: tmp_u_* can easily fit into shared memory!

    // Backward substitution
    for (uint32_t z = 1; z < depth; ++z) {
        const uint64_t xyz_offset = height * width * (depth - z - 1) +
                                    xy_offset;

        ftype f_x = tmp_f_x[xyz_offset];
        ftype f_y = tmp_f_y[xyz_offset];
        ftype f_z = tmp_f_z[xyz_offset];
        ftype upp = tmp[xyz_offset];

        ftype u_x_prev = tmp_u_x[xy_offset];
        ftype u_y_prev = tmp_u_y[xy_offset];
        ftype u_z_prev = tmp_u_z[xy_offset];

        ftype u_x = -upp * u_x_prev + f_x;
        ftype u_y = -upp * u_y_prev + f_y;
        ftype u_z = -upp * u_z_prev + f_z;

        tmp_u_x[xy_offset] = u_x;
        tmp_u_y[xy_offset] = u_y;
        tmp_u_z[xy_offset] = u_z;

        params.u_x[xyz_offset] += u_x;
        params.u_y[xyz_offset] += u_y;
        params.u_z[xyz_offset] += u_z;
    }
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

    // NOTE: The first transfer is redundant.
    CUDA_CHECK(cudaMemcpy(d_data.w, w, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_x, f_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_y, f_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.f_z, f_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_x, u_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_y, u_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data.u_z, u_z, alloc_size, CUDA_H2D));

    // Spawn one thread per point in the xy plane.

    const int num_threads = 128;
    const dim3 num_blocks((width - 1) / num_threads + 1, height);

    d_data.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    CUDA_CHECK_LAUNCH(
        solve_Dzz_sweep<BC><<<num_blocks, num_threads>>>(d_data));

    CUDA_CHECK(cudaMemcpy(u_x, d_data.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data.u_z, alloc_size, CUDA_D2H));

    /* No need to synchronize */
}
