#include "ftype.h"
#include "field.h"
#include "consts.h"
#include "timeit.h"

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

// TODO: Test 2D coarsening

#ifdef FLOAT
#define COARSE_FACTOR 4
#else
#define COARSE_FACTOR 2
#endif

#define P_128(addr) reinterpret_cast<_ftype128 *>((addr))

template<typename BoundaryConditions>
__global__ void solve_Dzz_sweep(__grid_constant__ const SolveDzzParams params)
{
    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t y = blockIdx.y;

    const uint32_t depth = params.depth;
    const uint32_t height = params.height;
    const uint32_t width = params.width;

    if (x > width / COARSE_FACTOR - 1) {
        return;
    }

    const _ftype128 *__restrict__ p_w = P_128(params.w);
    const _ftype128 *__restrict__ p_f_x = P_128(params.f_x);
    const _ftype128 *__restrict__ p_f_y = P_128(params.f_y);
    const _ftype128 *__restrict__ p_f_z = P_128(params.f_z);
    _ftype128 *__restrict__ p_u_x = P_128(params.u_x);
    _ftype128 *__restrict__ p_u_y = P_128(params.u_y);
    _ftype128 *__restrict__ p_u_z = P_128(params.u_z);

    const uint64_t num_points = depth * height * width;

    _ftype128 *__restrict__ tmp = P_128(params.tmp);
    // WARNING: Is 256 byte alignement guaranteed?
    _ftype128 *__restrict__ tmp_f_x = P_128(params.tmp + num_points);
    _ftype128 *__restrict__ tmp_f_y = P_128(params.tmp + num_points * 2);
    _ftype128 *__restrict__ tmp_f_z = P_128(params.tmp + num_points * 3);

    const uint64_t face_stride = height * width / COARSE_FACTOR;
    const uint64_t xy_offset = width / COARSE_FACTOR * y + x;

    const BoundaryConditions bc;

    _ftype128 u0_x_vec;
    _ftype128 u0_y_vec;
    _ftype128 u0_z_vec;
    _ftype128 zero_vec = { 0, 0 };

    // Apply BCs to the first face

    // TODO: I see branches in the SASS, are these being fully inlined?
    ftype3 u0_vec_0 = bc.get_u_delta_front(
        COARSE_FACTOR * x, y, 0, params.timestep);
    ftype3 u0_vec_1 = bc.get_u_delta_front(
        COARSE_FACTOR * x + 1, y, 0, params.timestep);

    u0_x_vec.x = u0_vec_0.x;
    u0_y_vec.x = u0_vec_0.y;
    u0_z_vec.x = u0_vec_0.z;

    u0_x_vec.y = u0_vec_1.x;
    u0_y_vec.y = u0_vec_1.y;
    u0_z_vec.y = u0_vec_1.z;

#ifdef FLOAT
    ftype3 u0_vec_2 = bc.get_u_delta_front(
        COARSE_FACTOR * x + 2, y, 0, params.timestep);
    ftype3 u0_vec_3 = bc.get_u_delta_front(
        COARSE_FACTOR * x + 3, y, 0, params.timestep);

    u0_x_vec.z = u0_vec_2.x;
    u0_y_vec.z = u0_vec_2.y;
    u0_z_vec.z = u0_vec_2.z;

    u0_x_vec.w = u0_vec_3.x;
    u0_y_vec.w = u0_vec_3.y;
    u0_z_vec.w = u0_vec_3.z;

    zero_vec.z = 0;
    zero_vec.w = 0;
#endif

    _ftype128 upp_redu = zero_vec;
    _ftype128 f_x_redu = u0_x_vec;
    _ftype128 f_y_redu = u0_y_vec;
    _ftype128 f_z_redu = u0_z_vec;

    tmp[xy_offset] = upp_redu; // set upper coefficient to 0
    // Enforce correct solution in rhs
    tmp_f_x[xy_offset] = f_x_redu;
    tmp_f_y[xy_offset] = f_y_redu;
    tmp_f_z[xy_offset] = f_z_redu;

    // Gauss reduce the remaining domain, except last face
    for (uint32_t z = 1; z < depth - 1; ++z) {
        const uint64_t xyz_offset =
            height * width / COARSE_FACTOR * z + xy_offset;

        _ftype128 w = p_w[xyz_offset];
        // TODO: Operator overloading?
        _ftype128 f_x = p_f_x[xyz_offset] - p_u_x[xyz_offset];
        _ftype128 f_y = p_f_y[xyz_offset] - p_u_y[xyz_offset];
        _ftype128 f_z = p_f_z[xyz_offset] - p_u_z[xyz_offset];

        _ftype128 norm_coef_inv = 1 / (w * upp_redu + 1 + 2 * w);

        upp_redu = -w * norm_coef_inv;
        f_x_redu = (w * f_x_redu + f_x) * norm_coef_inv;
        f_y_redu = (w * f_y_redu + f_y) * norm_coef_inv;
        f_z_redu = (w * f_z_redu + f_z) * norm_coef_inv;

        tmp[xyz_offset] = upp_redu;
        tmp_f_x[xyz_offset] = f_x_redu;
        tmp_f_y[xyz_offset] = f_y_redu;
        tmp_f_z[xyz_offset] = f_z_redu;
    }

    // Apply bcs to last face, solving directly

    const uint64_t xyz_offset =
        height * width / COARSE_FACTOR * (depth - 1) + xy_offset;

    _ftype128 w = p_w[xyz_offset];
    _ftype128 upp_prev = tmp[xyz_offset - face_stride];

    _ftype128 f_x = p_f_x[xyz_offset] - p_u_x[xyz_offset];
    _ftype128 f_y = p_f_y[xyz_offset] - p_u_y[xyz_offset];

    _ftype128 norm_coef_inv
        = 1 / (1 + 4 * w + (ftype) (4.0 / 3.0) * w * upp_prev);

    ftype3 un_0 = bc.get_u_delta_back(
        COARSE_FACTOR * x, y, depth - 1, params.timestep);
    ftype3 un_1 = bc.get_u_delta_back(
        COARSE_FACTOR * x + 1, y, depth - 1, params.timestep);

    _ftype128 un_x = { un_0.x, un_1.x };
    _ftype128 un_y = { un_0.y, un_1.y };
    _ftype128 un_z = { un_0.z, un_1.z };

#ifdef FLOAT
    ftype3 un_2 = bc.get_u_delta_back(
        COARSE_FACTOR * x + 2, y, depth - 1, params.timestep);
    ftype3 un_3 = bc.get_u_delta_back(
        COARSE_FACTOR * x + 3, y, depth - 1, params.timestep);

    un_x.z = un_2.x;
    un_y.z = un_2.y;
    un_z.z = un_2.z;
    un_x.w = un_3.x;
    un_y.w = un_3.y;
    un_z.w = un_3.z;
#endif

    un_x = (w * (ftype) (8.0 / 3.0) * un_x + f_x +
            (ftype) (4.0 / 3.0) * w * f_x_redu) * norm_coef_inv;
    un_y = (w * (ftype) (8.0 / 3.0) * un_y + f_y +
            (ftype) (4.0 / 3.0) * w * f_y_redu) * norm_coef_inv;

    _ftype128 u_x_prev = un_x;
    _ftype128 u_y_prev = un_y;
    _ftype128 u_z_prev = un_z;

    p_u_x[xyz_offset] += un_x;
    p_u_y[xyz_offset] += un_y;
    p_u_z[xyz_offset] += un_z;

    // Backward substitution
    for (uint32_t z = 1; z < depth; ++z) {
        const uint64_t xyz_offset =
            height * width / COARSE_FACTOR * (depth - z - 1) + xy_offset;

        _ftype128 f_x = tmp_f_x[xyz_offset];
        _ftype128 f_y = tmp_f_y[xyz_offset];
        _ftype128 f_z = tmp_f_z[xyz_offset];
        _ftype128 upp = tmp[xyz_offset];

        u_x_prev = -upp * u_x_prev + f_x;
        u_y_prev = -upp * u_y_prev + f_y;
        u_z_prev = -upp * u_z_prev + f_z;

        p_u_x[xyz_offset] += u_x_prev;
        p_u_y[xyz_offset] += u_y_prev;
        p_u_z[xyz_offset] += u_z_prev;
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

    const int num_threads = 64;
    const dim3 num_blocks((width / COARSE_FACTOR - 1) / num_threads + 1,
                          height);

    d_data.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    TIMER_CREATE(solve_Dzz_sweep);
    TIMER_RESTART(solve_Dzz_sweep);
    CUDA_CHECK_LAUNCH(
        solve_Dzz_sweep<BC><<<num_blocks, num_threads>>>(d_data));
    cudaDeviceSynchronize();
    TIMER_ELAPSED(solve_Dzz_sweep, 1);

    CUDA_CHECK(cudaMemcpy(u_x, d_data.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data.u_z, alloc_size, CUDA_D2H));

    /* No need to synchronize */
}
