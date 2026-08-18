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

#if __CUDACC_VER_MAJOR__> 12
#define GRID_CONSTANT __grid_constant__
#else
#define GRID_CONSTANT
#endif

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
                            domain_size.depth * 4 };
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
#define VEC_LENGTH 4
#else
#define VEC_LENGTH 2
#endif

// NOTE: The effect of coarsening may be visible
// on larger grids, 256^3 is not even saturating
// the 80 SMs of the v100
#define COARSE_FACTOR_Y 1

#define P_128(addr) reinterpret_cast<_ftype128 *>((addr))

struct DxxSweep {};
struct DyySweep {};
struct DzzSweep {};

template<typename BoundaryConditions, typename SweepDir>
struct SweepBoundaryConditions {};

template<typename BoundaryConditions>
struct SweepBoundaryConditions<BoundaryConditions, DyySweep>
{
    __device__ __forceinline__
    SweepBoundaryConditions(uint32_t depth_ [[maybe_unused]],
                            uint32_t height_,
                            uint32_t width_ [[maybe_unused]])
        : height(height_) {}

    __device__ __forceinline__
    ftype3 get_u_start(uint32_t x, uint32_t z, uint32_t t) const
    {
        return bc.get_u_delta_top(x, 0, z, t);
    }

    __device__ __forceinline__
    ftype3 get_u_end(uint32_t x, uint32_t z, uint32_t t) const
    {
        return bc.get_u_delta_bottom(x, height - 1, z, t);
    }

private:
    const BoundaryConditions bc;
    const uint32_t height;
};

template<typename BoundaryConditions>
struct SweepBoundaryConditions<BoundaryConditions, DzzSweep>
{
    __device__ __forceinline__
    SweepBoundaryConditions(uint32_t depth_,
                            uint32_t height_ [[maybe_unused]],
                            uint32_t width_ [[maybe_unused]])
        : depth(depth_) {}

    __device__ __forceinline__
    ftype3 get_u_start(uint32_t x, uint32_t y, uint32_t t) const
    {
        return bc.get_u_delta_front(x, y, 0, t);
    }

    __device__ __forceinline__
    ftype3 get_u_end(uint32_t x, uint32_t y, uint32_t t) const
    {
        return bc.get_u_delta_back(x, y, depth - 1, t);
    }

private:
    const BoundaryConditions bc;
    const uint32_t depth;
};

template<typename SweepDir>
struct SweepUtils {};

template<>
struct SweepUtils<DyySweep>
{
    __device__ __forceinline__
    SweepUtils(uint32_t depth_, uint32_t height_, uint32_t width_)
        : depth(depth_), height(height_), width(width_) {}

    __device__ __forceinline__ uint32_t face_height() const { return depth; }

    __device__ __forceinline__ uint32_t sweep_extent() const { return height; }

    __device__ __forceinline__ bool is_Dzz_sweep() const { return false; }

    __device__ __forceinline__ uint64_t global_idx(
        uint32_t x, uint32_t z, uint32_t y, uint32_t vec_len) const
    {
        return (height * width * z + width * y) / vec_len + x;
    }


private:
    const uint32_t depth, height, width;
};

template<>
struct SweepUtils<DzzSweep>
{
    __device__ __forceinline__
    SweepUtils(uint32_t depth_, uint32_t height_, uint32_t width_)
        : depth(depth_), height(height_), width(width_) {}

    __device__ __forceinline__ uint32_t face_height() const { return height; }

    __device__ __forceinline__ uint32_t sweep_extent() const { return depth; }

    __device__ __forceinline__ bool is_Dzz_sweep() const { return true; }

    __device__ __forceinline__ uint64_t global_idx(
        uint32_t x, uint32_t y, uint32_t z, uint32_t vec_len) const
    {
        // WARNING: x is already divided by vec_length
        return (height * width * z + width * y) / vec_len + x;
    }

private:
    const uint32_t depth, height, width;
};

// TODO: For tmp buffers, there's no strict need to proceed in the
// same direction of the f and u buffers, for the cpu version it
// could be helpful to avoid TLB misses?

// This needs to be wrapped in a struct to be partially specialized for Dzz
template<typename BoundaryConditions, typename SweepDir>
__global__ void solve_sweep(GRID_CONSTANT const SolveDzzParams params)
{
    using SweepUtils = SweepUtils<SweepDir>;
    using BC = SweepBoundaryConditions<BoundaryConditions, SweepDir>;

    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;

    const uint32_t depth = params.depth;
    const uint32_t height = params.height;
    const uint32_t width = params.width;

    if (x > width / VEC_LENGTH - 1) {
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

    const BC bc(depth, height, width);
    const SweepUtils sweep(depth, height, width);

    // NOTE: Coordinates are relative to the Dzz sweep.

    for (uint32_t y = blockIdx.y;
                  y < sweep.face_height(); y += gridDim.y) {

        _ftype128 u0_x_vec;
        _ftype128 u0_y_vec;
        _ftype128 u0_z_vec;
        _ftype128 zero_vec = { 0, 0 };

        // Apply BCs to the first face

        // TODO: I see branches in the SASS, are these being fully inlined?
        ftype3 u0_vec_0 = bc.get_u_start(
            VEC_LENGTH * x, y, params.timestep);
        ftype3 u0_vec_1 = bc.get_u_start(
            VEC_LENGTH * x + 1, y, params.timestep);

        u0_x_vec.x = u0_vec_0.x;
        u0_y_vec.x = u0_vec_0.y;
        u0_z_vec.x = u0_vec_0.z;

        u0_x_vec.y = u0_vec_1.x;
        u0_y_vec.y = u0_vec_1.y;
        u0_z_vec.y = u0_vec_1.z;

#ifdef FLOAT
        ftype3 u0_vec_2 = bc.get_u_start(
            VEC_LENGTH * x + 2, y, params.timestep);
        ftype3 u0_vec_3 = bc.get_u_start(
            VEC_LENGTH * x + 3, y, params.timestep);

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

        // TODO: Traversing front to back is slightly faster on V100 (~1%)
        const uint64_t tmp_offset = width / VEC_LENGTH * y + x;

        tmp[tmp_offset] = upp_redu; // set upper coefficient to 0
        // Enforce correct solution in rhs
        tmp_f_x[tmp_offset] = f_x_redu;
        tmp_f_y[tmp_offset] = f_y_redu;
        tmp_f_z[tmp_offset] = f_z_redu;

        // Gauss reduce the remaining domain, except last face
        for (uint32_t z = 1; z < sweep.sweep_extent() - 1; ++z) {
            const uint64_t xyz_offset =
                sweep.global_idx(x, y, z, VEC_LENGTH);

            _ftype128 w = p_w[xyz_offset];
            _ftype128 f_x = p_f_x[xyz_offset] - p_u_x[xyz_offset];
            _ftype128 f_y = p_f_y[xyz_offset] - p_u_y[xyz_offset];
            _ftype128 f_z = p_f_z[xyz_offset] - p_u_z[xyz_offset];

            _ftype128 norm_coef_inv = 1 / (w * upp_redu + 1 + 2 * w);

            upp_redu = -w * norm_coef_inv;
            f_x_redu = (w * f_x_redu + f_x) * norm_coef_inv;
            f_y_redu = (w * f_y_redu + f_y) * norm_coef_inv;
            f_z_redu = (w * f_z_redu + f_z) * norm_coef_inv;

            const uint64_t tmp_offset =
                (height * width * z + width * y) / VEC_LENGTH + x;

            tmp[tmp_offset] = upp_redu;
            tmp_f_x[tmp_offset] = f_x_redu;
            tmp_f_y[tmp_offset] = f_y_redu;
            tmp_f_z[tmp_offset] = f_z_redu;
        }

        // Apply bcs to last face, solving directly

        const uint64_t xyz_offset =
            sweep.global_idx(x, y, sweep.sweep_extent() - 1, VEC_LENGTH);

        _ftype128 w = p_w[xyz_offset];

        _ftype128 norm_coef_inv
            = 1 / (1 + 4 * w + (ftype) (4.0 / 3.0) * w * upp_redu);

        ftype3 un_0 = bc.get_u_end(VEC_LENGTH * x, y, params.timestep);
        ftype3 un_1 = bc.get_u_end(VEC_LENGTH * x + 1, y, params.timestep);

        _ftype128 un_x = { un_0.x, un_1.x };
        _ftype128 un_y = { un_0.y, un_1.y };
        _ftype128 un_z = { un_0.z, un_1.z };

#ifdef FLOAT
        ftype3 un_2 = bc.get_u_end(VEC_LENGTH * x + 2, y, params.timestep);
        ftype3 un_3 = bc.get_u_end(VEC_LENGTH * x + 3, y, params.timestep);

        un_x.z = un_2.x;
        un_y.z = un_2.y;
        un_z.z = un_2.z;
        un_x.w = un_3.x;
        un_y.w = un_3.y;
        un_z.w = un_3.z;
#endif

        _ftype128 f_x = p_f_x[xyz_offset] - p_u_x[xyz_offset];

        // TODO: Can you wrap this inside SweepBoundaryConditions?
        if (sweep.is_Dzz_sweep()) {
            _ftype128 f_y = p_f_y[xyz_offset] - p_u_y[xyz_offset];
            un_y = (w * (ftype) (8.0 / 3.0) * un_y + f_y +
                (ftype) (4.0 / 3.0) * w * f_y_redu) * norm_coef_inv;
        } else { // DyySweep
            _ftype128 f_z = p_f_z[xyz_offset] - p_u_z[xyz_offset];
            un_z = (w * (ftype) (8.0 / 3.0) * un_z + f_z +
                (ftype) (4.0 / 3.0) * w * f_z_redu) * norm_coef_inv;
        }

        un_x = (w * (ftype) (8.0 / 3.0) * un_x + f_x +
                (ftype) (4.0 / 3.0) * w * f_x_redu) * norm_coef_inv;

        _ftype128 u_x_prev = un_x;
        _ftype128 u_y_prev = un_y;
        _ftype128 u_z_prev = un_z;

        p_u_x[xyz_offset] += un_x;
        p_u_y[xyz_offset] += un_y;
        p_u_z[xyz_offset] += un_z;

        // Backward substitution
        for (uint32_t z = 1; z < sweep.sweep_extent(); ++z) {
            const uint64_t tmp_offset =
                (height * width * (sweep.sweep_extent() - z - 1) +
		 width * y) / VEC_LENGTH + x;

            _ftype128 f_x = tmp_f_x[tmp_offset];
            _ftype128 f_y = tmp_f_y[tmp_offset];
            _ftype128 f_z = tmp_f_z[tmp_offset];
            _ftype128 upp = tmp[tmp_offset];

            u_x_prev = -upp * u_x_prev + f_x;
            u_y_prev = -upp * u_y_prev + f_y;
            u_z_prev = -upp * u_z_prev + f_z;

            const uint64_t xyz_offset = sweep.global_idx(
                x, y, sweep.sweep_extent() - z - 1, VEC_LENGTH);

            p_u_x[xyz_offset] += u_x_prev;
            p_u_y[xyz_offset] += u_y_prev;
            p_u_z[xyz_offset] += u_z_prev;
        }
    }
}

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

    // Spawn one thread per point in the xz plane.

    const int num_threads = 64;
    const dim3 num_blocks((width / VEC_LENGTH - 1) / num_threads + 1,
                          depth / COARSE_FACTOR_Y);

    d_data.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    CUDA_TIMER_CREATE(solve_Dyy_sweep);
    CUDA_TIMER_RESTART(solve_Dyy_sweep);
    CUDA_CHECK_LAUNCH(
        solve_sweep<BC, DyySweep><<<num_blocks, num_threads>>>(d_data));
    cudaDeviceSynchronize();
    CUDA_TIMER_ELAPSED(solve_Dyy_sweep, 1);
    CUDA_TIMER_DESTROY(solve_Dyy_sweep);

    CUDA_CHECK(cudaMemcpy(u_x, d_data.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data.u_z, alloc_size, CUDA_D2H));

    // No need to synchronize
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
    const dim3 num_blocks((width / VEC_LENGTH - 1) / num_threads + 1,
                          height / COARSE_FACTOR_Y);

    d_data.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    CUDA_TIMER_CREATE(solve_Dzz_sweep);
    CUDA_TIMER_RESTART(solve_Dzz_sweep);
    CUDA_CHECK_LAUNCH(
        solve_sweep<BC, DzzSweep><<<num_blocks, num_threads>>>(d_data));
    cudaDeviceSynchronize();
    CUDA_TIMER_ELAPSED(solve_Dzz_sweep, 1);
    CUDA_TIMER_DESTROY(solve_Dzz_sweep);

    CUDA_CHECK(cudaMemcpy(u_x, d_data.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data.u_z, alloc_size, CUDA_D2H));

    // No need to synchronize
}
