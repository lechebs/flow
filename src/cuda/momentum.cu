#include "ftype.h"
#include "field.h"
#include "consts.h"
#include "timeit.h"

/* cuda/ includes */
#include "utils.h"
#include "forcing.h"
#include "boundary.h"
#include "momentum.h"

// TODO: Consider templating kernel on ftype

#define CUDA_H2D cudaMemcpyHostToDevice
#define CUDA_D2H cudaMemcpyDeviceToHost

#define WARP_SIZE 32

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

struct SolveDxxParams
{
    ftype *__restrict__ k;
    ftype *__restrict__ w;
    ftype *__restrict__ p;
    ftype *__restrict__ phi;
    ftype *__restrict__ tmp;
    ftype *__restrict__ eta_x;
    ftype *__restrict__ eta_y;
    ftype *__restrict__ eta_z;
    ftype *__restrict__ zeta_x;
    ftype *__restrict__ zeta_y;
    ftype *__restrict__ zeta_z;
    ftype *__restrict__ u_x;
    ftype *__restrict__ u_y;
    ftype *__restrict__ u_z;

    uint32_t depth, height, width;
    uint32_t timestep;
};

static struct SolveDzzParams d_data_dzz;
static struct SolveDxxParams d_data_dxx;

void alloc_device_data(field_size domain_size)
{
    uint64_t alloc_size = field_num_points(domain_size) * sizeof(ftype);

    // WARNING: Some of these are padded in the cpu version.
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.w, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.f_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.f_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.f_z, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.u_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.u_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.u_z, alloc_size));

    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.k, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.p, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.phi, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.eta_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.eta_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.eta_z, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.zeta_x, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.zeta_y, alloc_size));
    CUDA_CHECK(cudaMalloc((void **) &d_data_dxx.zeta_z, alloc_size));

    d_data_dxx.w = d_data_dzz.w;
    d_data_dxx.u_x = d_data_dzz.u_x;
    d_data_dxx.u_y = d_data_dzz.u_y;
    d_data_dxx.u_z = d_data_dzz.u_z;

    field_size tmp_size = { domain_size.width,
                            domain_size.height,
                            domain_size.depth * 4 };
    uint64_t tmp_alloc_size = field_num_points(tmp_size) * sizeof(ftype);

    CUDA_CHECK(cudaMalloc((void **) &d_data_dzz.tmp, tmp_alloc_size));

    d_data_dxx.tmp = d_data_dzz.tmp;

    d_data_dzz.depth = domain_size.depth;
    d_data_dzz.height = domain_size.height;
    d_data_dzz.width = domain_size.width;

    d_data_dxx.depth = domain_size.depth;
    d_data_dxx.height = domain_size.height;
    d_data_dxx.width = domain_size.width;
}

void update_device_consts()
{
    CUDA_CHECK(cudaMemcpyToSymbol(_D_DX, &_DX, sizeof(ftype), 0, CUDA_H2D));
    CUDA_CHECK(cudaMemcpyToSymbol(_D_DT, &_DT, sizeof(ftype), 0, CUDA_H2D));
    CUDA_CHECK(cudaMemcpyToSymbol(_D_NU, &_NU, sizeof(ftype), 0, CUDA_H2D));

    ftype dx_rcp = 1 / _DX;
    CUDA_CHECK(cudaMemcpyToSymbol(_D_DX_RCP, &dx_rcp, sizeof(ftype), 0, CUDA_H2D));
}

void free_device_data()
{
    CUDA_CHECK(cudaFree(d_data_dzz.w));
    CUDA_CHECK(cudaFree(d_data_dzz.tmp));
    CUDA_CHECK(cudaFree(d_data_dzz.f_x));
    CUDA_CHECK(cudaFree(d_data_dzz.f_y));
    CUDA_CHECK(cudaFree(d_data_dzz.f_z));
    CUDA_CHECK(cudaFree(d_data_dzz.u_x));
    CUDA_CHECK(cudaFree(d_data_dzz.u_y));
    CUDA_CHECK(cudaFree(d_data_dzz.u_z));

    CUDA_CHECK(cudaFree(d_data_dxx.k));
    CUDA_CHECK(cudaFree(d_data_dxx.p));
    CUDA_CHECK(cudaFree(d_data_dxx.phi));
    CUDA_CHECK(cudaFree(d_data_dxx.eta_x));
    CUDA_CHECK(cudaFree(d_data_dxx.eta_y));
    CUDA_CHECK(cudaFree(d_data_dxx.eta_z));
    CUDA_CHECK(cudaFree(d_data_dxx.zeta_x));
    CUDA_CHECK(cudaFree(d_data_dxx.zeta_y));
    CUDA_CHECK(cudaFree(d_data_dxx.zeta_z));
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
struct SweepBoundaryConditions<BoundaryConditions, DxxSweep>
{
    __device__ __forceinline__
    SweepBoundaryConditions(uint32_t depth_ [[maybe_unused]],
                            uint32_t height_ [[maybe_unused]],
                            uint32_t width_)
        : width(width_) {}

    __device__ __forceinline__
    ftype3 get_u_start(uint32_t y, uint32_t z, uint32_t t) const
    {
        return bc.get_u_delta_left(0, y, z, t);
    }

    __device__ __forceinline__
    ftype3 get_u_end(uint32_t y, uint32_t z, uint32_t t) const
    {
        return bc.get_u_delta_right(width - 1, y, z, t);
    }

private:
    const BoundaryConditions bc;
    const uint32_t width;
};

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

    // WARNING: If coarsening is used, the amount of tmp
    // storage used can be reduced! In that case indexing
    // differs by a constant offset.
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

        // Traversing front to back is slightly faster on V100 (~1%)
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

struct WarpTile
{
    __device__ __forceinline__ WarpTile() {}

    __device__ __forceinline__
    void load(const _ftype128 *__restrict__ src, uint32_t stride);

    __device__ __forceinline__
    void store_add(_ftype128 *__restrict__ dst, uint32_t stride);

    __device__ __forceinline__ void transpose();

    __device__ __forceinline__ void zero()
    {
        #pragma unroll
        for (int k = 0; k < 8; ++k) {
            frag[k].x = 0;
            frag[k].y = 0;
#ifdef FLOAT
            frag[k].z = 0;
            frag[k].w = 0;
#endif
        }
    }

    template<int comp>
    __device__ __forceinline__ ftype &col(int idx)
    {
#ifdef FLOAT
        if constexpr (comp == 0) {
            return frag[idx].x;
        } else if constexpr (comp == 1) {
            return frag[idx].y;
        } else if constexpr (comp == 2) {
            return frag[idx].z;
        } else {
            return frag[idx].w;
        }
#else
        if constexpr (comp == 0) {
            return frag[idx].x;
        } else {
            return frag[idx].y;
        }
#endif
    }

    __device__ __forceinline__ _ftype128 &cols(int idx)
    {
        return frag[idx];
    }

private:
    _ftype128 frag[8];
};

// Specialize WarpTile for shmem or regs storage

__device__ __forceinline__
void WarpTile::load(const _ftype128 *__restrict__ src,
                    uint32_t stride)
{
    const int lane_idx = threadIdx.x % WARP_SIZE;

    const _ftype128 *__restrict__ t_src =
        src + (lane_idx / 8) * 8 * stride / VEC_LENGTH + (lane_idx % 8);

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        frag[i] = *t_src;
        t_src += stride / VEC_LENGTH;
    }
}

__device__ __forceinline__
void WarpTile::store_add(_ftype128 *__restrict__ dst,
                         uint32_t stride)
{
    const int lane_idx = threadIdx.x % WARP_SIZE;

    _ftype128 *__restrict__ t_dst =
        dst + (lane_idx / 8) * 8 * stride / VEC_LENGTH + (lane_idx % 8);

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        *t_dst += frag[i];
        t_dst += stride / VEC_LENGTH;
    }
}

__device__ __forceinline__ _ftype128 shfl_xor_128(_ftype128 var,
                                                  int lane_mask)
{
    _ftype128 res;
    res.x = __shfl_xor_sync(0xffffffff, var.x, lane_mask);
    res.y = __shfl_xor_sync(0xffffffff, var.y, lane_mask);
#ifdef FLOAT
    res.z = __shfl_xor_sync(0xffffffff, var.z, lane_mask);
    res.w = __shfl_xor_sync(0xffffffff, var.w, lane_mask);
#endif
    return res;
}

__device__ __forceinline__ void WarpTile::transpose()
{
    // Keep frag in registers!

    // We are not actually tranposing the tile,
    // rather transposing each 8x8 subtile of 128 bit words

    // Each quarter warp transposes its subtile

    // TODO: You can use one shuffle by predicating
    // both the src and dst register choices

    #pragma unroll
    for (int s = 1; s <= 4; s <<= 1) {

        #pragma unroll
        for (int i = s; i < 8; i += s * 2) {
            #pragma unroll
            for (int ii = 0; ii < s; ++ii) {
                frag[i + ii] = shfl_xor_128(frag[i + ii], s);
            }
        }

        if ((threadIdx.x / s) % 2 != 0) {
            #pragma unroll
            for (int i = 0; i < 8; i += s * 2) {
                #pragma unroll
                for (int ii = 0; ii < s; ++ii) {
                    _ftype128 tmp = frag[i + ii];
                    frag[i + ii] = frag[i + ii + s];
                    frag[i + ii + s] = tmp;
                }
            }
        }

        #pragma unroll
        for (int i = s; i < 8; i += s * 2) {
            #pragma unroll
            for (int ii = 0; ii < s; ++ii) {
                frag[i + ii] = shfl_xor_128(frag[i + ii], s);
            }
        }
    }
}

__device__ __forceinline__
void compute_pp_grad(const _ftype128 *__restrict__ p,
                     const _ftype128 *__restrict__ phi,
                     uint32_t height,
                     uint32_t width,
                     bool is_last_face,
                     bool is_last_row,
                     bool is_last_col,
                     WarpTile &pp_x, // transposed on return
                     WarpTile &pp_y,
                     WarpTile &pp_z)
{
    // TODO: Receive pressure predictor directly

    // Load p tile once, aim for reuse
    pp_x.load(p, width);
    // Something like load_add could be useful, currently
    // create another WarpTile hoping that the compiler
    // will fuse the loops
    WarpTile phi_t;
    phi_t.load(phi, width);

    // Compute pressure predictor
    #pragma unroll
    for (int k = 0; k < 8; ++k) {
        pp_x.cols(k) += phi_t.cols(k);
    }

    if (is_last_face) {
        pp_z.zero();
    } else {
        // Load next tile in the z direction
        pp_z.load(p + height * width / VEC_LENGTH, width);
        phi_t.load(phi + height * width / VEC_LENGTH, width);
        // Compute pressure predictor
        #pragma unroll
        for (int k = 0; k < 8; ++k) {
            pp_z.cols(k) += phi_t.cols(k);
        }
        // Compute dp/dz
        #pragma unroll
        for (int k = 0; k < 8; ++k) {
            pp_z.cols(k) = (pp_z.cols(k) - pp_x.cols(k)) * -_D_DX_RCP;
        }
    }

    // Compute dp/dy
    #pragma unroll
    for (int k = 0; k < 7; ++k) {
        pp_y.cols(k) = (pp_x.cols(k + 1) - pp_x.cols(k)) * -_D_DX_RCP;
    }

    const int lane_idx = threadIdx.x % WARP_SIZE;
    const int quarter = lane_idx / 8;
    const int strip = lane_idx % 8;

    if (is_last_row && quarter == 3) {
        pp_y.col<0>(7) = 0;
        pp_y.col<1>(7) = 0;
#ifdef FLOAT
        pp_y.col<2>(7) = 0;
        pp_y.col<3>(7) = 0;
#endif
    } else {
        // TODO: shuffle for q < 3, gmem only for q == 3.
        const uint64_t next_row =
            (uint64_t) width / VEC_LENGTH * 8 * (quarter + 1) + strip;
        pp_y.cols(7) =
            (*(p + next_row) + *(phi + next_row) -
             pp_x.cols(7)) * -_D_DX_RCP;
    }

    pp_x.transpose();

    // Compute dp/dx
    #pragma unroll
    for (int k = 0; k < 8; ++k) {
        pp_x.col<0>(k) = (pp_x.col<1>(k) - pp_x.col<0>(k)) * -_D_DX_RCP;
#ifdef FLOAT
        pp_x.col<1>(k) = (pp_x.col<2>(k) - pp_x.col<1>(k)) * -_D_DX_RCP;
        pp_x.col<2>(k) = (pp_x.col<3>(k) - pp_x.col<2>(k)) * -_D_DX_RCP;
#endif

        if (k < 7) { // should be folded, right?
            pp_x.col<VEC_LENGTH - 1>(k) =
                (pp_x.col<0>(k + 1) -
                 pp_x.col<VEC_LENGTH - 1>(k)) * -_D_DX_RCP;
        } else {
            if (is_last_col) {
                pp_x.col<VEC_LENGTH - 1>(7) = 0;
            } else {
                const int lane_idx = threadIdx.x % WARP_SIZE;

                // WARNING: gather across 32 rows :/
                ftype pp_next = *reinterpret_cast<const ftype *>(
                    p + width / VEC_LENGTH * lane_idx + 8) +
                    *reinterpret_cast<const ftype *>(
                        phi + width / VEC_LENGTH * lane_idx + 8);

                pp_x.col<VEC_LENGTH - 1>(7) =
                    (pp_next - pp_x.col<VEC_LENGTH - 1>(7)) * -_D_DX_RCP;
            }
        }
    }

    pp_x.transpose();
}

#define COMP_X 0
#define COMP_Y 1
#define COMP_Z 2

template<typename BoundaryConditions, int comp> // x, y, z rhs components
__device__ __forceinline__
void compute_Dxx_Dyy_Dzz(const _ftype128 *__restrict__ eta,
                         const _ftype128 *__restrict__ zeta,
                         const _ftype128 *__restrict__ vel,
                         const BoundaryConditions &bc,
                         uint32_t depth,
                         uint32_t height,
                         uint32_t width,
                         // tile absolute coordinates
                         uint32_t tk,
                         uint32_t tj,
                         uint32_t ti,
                         uint32_t timestep,
                         ftype &eta_prev,
                         WarpTile &dst)
{
    const int lane_idx = threadIdx.x % WARP_SIZE;

    // Compute Dzz vel

    const bool is_first_face = ti == 0;
    const bool is_last_face = ti == depth - 1;

    WarpTile dd;
    dd.load(vel, width);

    // TODO: Can you accumulate directly in dst
    // to avoid allocating dd?

    if (is_first_face || (is_last_face && comp == COMP_Z)) {
        return;
    } else if (is_last_face) {

        WarpTile v;
        v.load(vel - height * width / VEC_LENGTH, width);

        // Something similar to load_add would be useful,
        // perhaps something more generic. Like a macro
        // or a template function that expands to a load
        // and in place add/sub to avoid allocating more
        // registers.
        // Currently let's just create another WarpTile
        // and hope that the compiler fuses the load
        // and finite difference loops.

        #pragma unroll
        for (int j = 0; j < 8; ++j) {

            dd.cols(j) = (ftype) 4.0 / 3 * v.cols(j) - 4 * dd.cols(j);

            ftype bc_vals[VEC_LENGTH];
            #pragma unroll
            for (int c = 0; c < VEC_LENGTH; ++c) {
                ftype3 val =
                    bc.get_u_back(tk + (lane_idx % 8) * VEC_LENGTH + c,
                                  tj + (lane_idx / 8) * 8 + j,
                                  depth - 1,
                                  timestep - 1);
                if constexpr (comp == COMP_X) {
                    bc_vals[c] = val.x;
                } else { // comp == COMP_Y
                    bc_vals[c] = val.y;
                }
            }

            dd.col<0>(j) += (ftype) 8.0 / 3 * bc_vals[0];
            dd.col<1>(j) += (ftype) 8.0 / 3 * bc_vals[1];
#ifdef FLOAT
            dd.col<2>(j) += (ftype) 8.0 / 3 * bc_vals[2];
            dd.col<3>(j) += (ftype) 8.0 / 3 * bc_vals[3];
#endif
        }
    } else {
        // TODO: Also here, a load_add would be useful,
        // just like the cpu version
        WarpTile v_prev, v_next;
        v_prev.load(vel - height * width / VEC_LENGTH, width);
        v_next.load(vel + height * width / VEC_LENGTH, width);

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            dd.cols(j) = v_next.cols(j) + -2 * dd.cols(j) + v_prev.cols(j);
        }
    }

    // Compute Dyy zeta

    const bool is_first_row = tj == 0;
    const bool is_last_row = tj == height - WARP_SIZE;

    WarpTile z;
    z.load(zeta, width);

    const int quarter = lane_idx / 8;
    const int strip = lane_idx % 8;

    if (!is_first_row || quarter > 0) {
        _ftype128 zprev =
            *(zeta + width / VEC_LENGTH * 8 * quarter
                   - width / VEC_LENGTH + strip);

        dd.cols(0) += zprev - 2 * z.cols(0) + z.cols(1);
    }

    #pragma unroll
    for (int j = 1; j < 7; ++j) {
        dd.cols(j) += z.cols(j - 1) - 2 * z.cols(j) + z.cols(j + 1);
    }

    if (is_last_row && comp != COMP_Y && lane_idx > 23) { // divergence :/

        ftype bc_vals[VEC_LENGTH];
        #pragma unroll
        for (int c = 0; c < VEC_LENGTH; ++c) {
            ftype3 val =
                bc.get_u_bottom(tk + (lane_idx % 8) * VEC_LENGTH + c,
                                height - 1,
                                ti,
                                timestep - 1);
            if constexpr (comp == COMP_X) {
                bc_vals[c] = val.x;
            } else { // comp == COMP_Z
                bc_vals[c] = val.z;
            }
        }

        dd.cols(7) += (ftype) 4.0 / 3 * z.cols(6) - 4 * z.cols(7);

        dd.col<0>(7) += (ftype) 8.0 / 3 * bc_vals[0];
        dd.col<1>(7) += (ftype) 8.0 / 3 * bc_vals[1];
#ifdef FLOAT
        dd.col<2>(7) += (ftype) 8.0 / 3 * bc_vals[2];
        dd.col<3>(7) += (ftype) 8.0 / 3 * bc_vals[3];
#endif

    } else if (!is_last_row || lane_idx < 24) {
        // TODO: only last subwarp has to fetch from gmem,
        // others can use shuffles
        dd.cols(7) += z.cols(6) - 2 * z.cols(7) +
                      *(zeta + width / VEC_LENGTH * (8 * (lane_idx / 8 + 1))
                        + (lane_idx % 8));
    }

    const bool is_first_col = tk == 0;
    const bool is_last_col = tk == width - 8 * VEC_LENGTH;

    // compute Dxx eta

    WarpTile e;
    e.load(eta, width);
    e.transpose();
    // cpu version can avoid allocating another
    // tile by storing eta finite difference one
    // shifted one column to the left.
    // That would be a bit cumbersome to handle
    // with register tiles here, so let's avoid it.
    WarpTile ddx;

    const ftype left_col = is_first_col ? 0 : eta_prev;

    ftype right_col = 0;
    if (is_last_col) {

        if constexpr (comp != COMP_X) {
            ftype3 u = bc.get_u_right(
                width - 1, tj + lane_idx, ti, timestep - 1);

            if constexpr (comp == COMP_Y) {
                right_col = u.y;
            } else { // comp == COMP_Z
                right_col = u.z;
            }
        }
    } else {
        right_col = (eta + lane_idx * (width / VEC_LENGTH) + 8)->x;
    }

    // k = 0
    ddx.col<0>(0) = left_col - 2 * e.col<0>(0) + e.col<1>(0);
#ifdef FLOAT
    ddx.col<1>(0) = e.col<0>(0) - 2 * e.col<1>(0) + e.col<2>(0);
    ddx.col<2>(0) = e.col<1>(0) - 2 * e.col<2>(0) + e.col<3>(0);
#endif
    ddx.col<VEC_LENGTH - 1>(0) = e.col<VEC_LENGTH - 2>(0) -
                                 2 * e.col<VEC_LENGTH - 1>(0) +
                                 e.col<0>(1);

    #pragma unroll
    for (int k = 1; k < 7; ++k) {

        ddx.col<0>(k) =
            e.col<VEC_LENGTH - 1>(k - 1) - 2 * e.col<0>(k) + e.col<1>(k);

        ddx.col<1>(k) = e.col<0>(k) - 2 * e.col<1>(k);
#ifdef FLOAT
        ddx.col<1>(k) += e.col<2>(k);
#else
        ddx.col<1>(k) += e.col<0>(k + 1);
#endif

#ifdef FLOAT
        ddx.col<2>(k) = e.col<1>(k) - 2 * e.col<2>(k) + e.col<3>(k);
        ddx.col<3>(k) = e.col<2>(k) - 2 * e.col<3>(k) + e.col<0>(k + 1);
#endif
    }

    // k = 7

    ddx.col<0>(7) = e.col<VEC_LENGTH - 1>(6) - 2 * e.col<0>(7) + e.col<1>(7);
#ifdef FLOAT
    ddx.col<1>(7) = e.col<0>(7) - 2 * e.col<1>(7) + e.col<2>(7);
    ddx.col<2>(7) = e.col<1>(7) - 2 * e.col<2>(7) + e.col<3>(7);
#endif

    if (is_last_col && comp != COMP_X) {
        ddx.col<VEC_LENGTH - 1>(7) =
            (ftype) (4.0 / 3) * e.col<VEC_LENGTH - 2>(7) -
            4 * e.col<VEC_LENGTH - 1>(7) +
            (ftype) (8.0 / 3) * right_col;
    } else {
        ddx.col<VEC_LENGTH - 1>(7) =
            e.col<VEC_LENGTH - 2>(7) - 2 * e.col<VEC_LENGTH - 1>(7) + right_col;
    }

    eta_prev = e.col<VEC_LENGTH - 1>(7);

    ddx.transpose();

    #pragma unroll
    for (int j = 0; j < 8; ++j) {
        dst.cols(j) += (dd.cols(j) + ddx.cols(j)) *
                       _D_NU * _D_DX_RCP * _D_DX_RCP;
    }
}

struct ftype128_vec
{
    __device__ __forceinline__ ftype128_vec(const ftype128_vec &v)
        : x(v.x), y(v.y), z(v.z) {}

    __device__ __forceinline__ ftype128_vec(const _ftype128 *x_,
                                            const _ftype128 *y_,
                                            const _ftype128 *z_)
        : x(x_), y(y_), z(z_) {}

    __device__ __forceinline__ ftype128_vec(const ftype *x_,
                                            const ftype *y_,
                                            const ftype *z_)
        : x(reinterpret_cast<const _ftype128 *>(x_)),
          y(reinterpret_cast<const _ftype128 *>(y_)),
          z(reinterpret_cast<const _ftype128 *>(z_)) {}

    __device__ __forceinline__ ftype128_vec operator+(uint64_t offset) const
    {
        return ftype128_vec(x + offset, y + offset, z + offset);
    }

    __device__ __forceinline__
    _ftype128 *x_cc() { return const_cast<_ftype128 *>(x); }

    __device__ __forceinline__
    _ftype128 *y_cc() { return const_cast<_ftype128 *>(y); }

    __device__ __forceinline__
    _ftype128 *z_cc() { return const_cast<_ftype128 *>(z); }

    const _ftype128 *__restrict__ x;
    const _ftype128 *__restrict__ y;
    const _ftype128 *__restrict__ z;
};

template<typename BoundaryConditions, typename ForcingTerm>
__device__ __forceinline__
void compute_Dxx_rhs_tile(const _ftype128 *__restrict__ k,
                          const _ftype128 *__restrict__ p,
                          const _ftype128 *__restrict__ phi,
                          const ftype128_vec eta,
                          const ftype128_vec zeta,
                          const ftype128_vec vel,
                          const BoundaryConditions &bc,
                          ForcingTerm f,
                          uint32_t depth,
                          uint32_t height,
                          uint32_t width,
                          // tile absolute coordinates
                          uint32_t tk,
                          uint32_t tj,
                          uint32_t ti,
                          uint32_t timestep,
                          ftype3 &eta_prev_last,
                          WarpTile &rhs_x,
                          WarpTile &rhs_y,
                          WarpTile &rhs_z)
{
    compute_pp_grad(p, phi, height, width,
                    ti == depth - 1,
                    tj == height - WARP_SIZE,
                    tk == width - 8 * VEC_LENGTH,
                    rhs_x, rhs_y, rhs_z);

    compute_Dxx_Dyy_Dzz<BoundaryConditions, COMP_X>(
        eta.x, zeta.x, vel.x, bc, depth, height, width,
        tk, tj, ti, timestep, eta_prev_last.x, rhs_x);
    compute_Dxx_Dyy_Dzz<BoundaryConditions, COMP_Y>(
        eta.y, zeta.y, vel.y, bc, depth, height, width,
        tk, tj, ti, timestep, eta_prev_last.y, rhs_y);
     compute_Dxx_Dyy_Dzz<BoundaryConditions, COMP_Z>(
        eta.z, zeta.z, vel.z, bc, depth, height, width,
        tk, tj, ti, timestep, eta_prev_last.z, rhs_z);

    const int lane_idx = threadIdx.x % WARP_SIZE;

    // I'm avoiding explicit WarpTiles here

    #pragma unroll
    for (int j = 0; j < 8; ++j) {

        const uint64_t off =
            width / VEC_LENGTH * ((lane_idx / 8) * 8 + j) + (lane_idx % 8);

        _ftype128 k_ = k[off];

        _ftype128 dt_over_beta =
            _D_DT / ((ftype) 1 + _D_DT * _D_NU / (2 * k_));

        _ftype128 coeff =
            (ftype) 1 - 2 * _D_DT * _D_NU / (2 * k_ + _D_DT * _D_NU);

        const uint32_t x = tk + VEC_LENGTH * (lane_idx % 8);
        const uint32_t y = tj + (lane_idx / 8) * 8 + j;
        const uint32_t z = ti;

        _ftype128 f_x, f_y, f_z;

        f_x.x = f.x(x + 0, y, z, timestep);
        f_y.x = f.y(x + 0, y, z, timestep);
        f_z.x = f.z(x + 0, y, z, timestep);
        f_x.y = f.x(x + 1, y, z, timestep);
        f_y.y = f.y(x + 1, y, z, timestep);
        f_z.y = f.z(x + 1, y, z, timestep);
#ifdef FLOAT
        f_x.z = f.x(x + 2, y, z, timestep);
        f_y.z = f.y(x + 2, y, z, timestep);
        f_z.z = f.z(x + 2, y, z, timestep);
        f_x.w = f.x(x + 3, y, z, timestep);
        f_y.w = f.y(x + 3, y, z, timestep);
        f_z.w = f.z(x + 3, y, z, timestep);
#endif

        rhs_x.cols(j) = coeff * vel.x[off]
                        + (rhs_x.cols(j) + f_x) * dt_over_beta
                        - eta.x[off];

        rhs_y.cols(j) = coeff * vel.y[off] +
                        (rhs_y.cols(j) + f_y) * dt_over_beta
                        - eta.y[off];

        rhs_z.cols(j) = coeff * vel.z[off] +
                        (rhs_z.cols(j) + f_z) * dt_over_beta -
                        eta.z[off];
    }

    rhs_x.transpose();
    rhs_y.transpose();
    rhs_z.transpose();
}

struct FwReduParams
{
    WarpTile w_t;
    WarpTile f_x_t;
    WarpTile f_y_t;
    WarpTile f_z_t;

    ftype *__restrict__ tmp_upp;
    ftype *__restrict__ tmp_f_x;
    ftype *__restrict__ tmp_f_y;
    ftype *__restrict__ tmp_f_z;

    ftype upp_redu;
    ftype f_x_redu;
    ftype f_y_redu;
    ftype f_z_redu;

    uint32_t tmp_stride;
};

__device__ __forceinline__ void fw_store_redu(FwReduParams &params)
{
    *params.tmp_upp = params.upp_redu;
    *params.tmp_f_x = params.f_x_redu;
    *params.tmp_f_y = params.f_y_redu;
    *params.tmp_f_z = params.f_z_redu;

    params.tmp_upp += params.tmp_stride;
    params.tmp_f_x += params.tmp_stride;
    params.tmp_f_y += params.tmp_stride;
    params.tmp_f_z += params.tmp_stride;
}

template<int comp> // float4/double2 component
__device__ __forceinline__ void fw_redu_col(FwReduParams &params, int col)
{
    ftype w = params.w_t.col<comp>(col);
    ftype f_x = params.f_x_t.col<comp>(col);
    ftype f_y = params.f_y_t.col<comp>(col);
    ftype f_z = params.f_z_t.col<comp>(col);

    ftype norm_coef_inv = 1 / (w * params.upp_redu + 1 + 2 * w);

    params.upp_redu = -w * norm_coef_inv;
    params.f_x_redu = (w * params.f_x_redu + f_x) * norm_coef_inv;
    params.f_y_redu = (w * params.f_y_redu + f_y) * norm_coef_inv;
    params.f_z_redu = (w * params.f_z_redu + f_z) * norm_coef_inv;

    fw_store_redu(params);
}

__device__ __forceinline__ void fw_redu_cols(FwReduParams &params, int col)
{
    fw_redu_col<0>(params, col);
    fw_redu_col<1>(params, col);
#ifdef FLOAT
    fw_redu_col<2>(params, col);
    fw_redu_col<3>(params, col);
#endif
}

struct BwSubParams
{
    WarpTile u_x_t;
    WarpTile u_y_t;
    WarpTile u_z_t;

    ftype *__restrict__ tmp_upp;
    ftype *__restrict__ tmp_f_x;
    ftype *__restrict__ tmp_f_y;
    ftype *__restrict__ tmp_f_z;

    ftype u_x_prev;
    ftype u_y_prev;
    ftype u_z_prev;

    uint32_t tmp_stride;
};

template<int comp>
__device__ __forceinline__ void bw_sub_col(BwSubParams &params, int col)
{
    ftype upp = *params.tmp_upp;
    ftype f_x = *params.tmp_f_x;
    ftype f_y = *params.tmp_f_y;
    ftype f_z = *params.tmp_f_z;

    // NOTE: I could avoid *_prev variables
    params.u_x_prev = -upp * params.u_x_prev + f_x;
    params.u_y_prev = -upp * params.u_y_prev + f_y;
    params.u_z_prev = -upp * params.u_z_prev + f_z;

    params.u_x_t.col<comp>(col) = params.u_x_prev;
    params.u_y_t.col<comp>(col) = params.u_y_prev;
    params.u_z_t.col<comp>(col) = params.u_z_prev;

    params.tmp_upp -= params.tmp_stride;
    params.tmp_f_x -= params.tmp_stride;
    params.tmp_f_y -= params.tmp_stride;
    params.tmp_f_z -= params.tmp_stride;
}

__device__ __forceinline__ void bw_sub_cols(BwSubParams &params, int col)
{
#ifdef FLOAT
    bw_sub_col<3>(params, col);
    bw_sub_col<2>(params, col);
#endif
    bw_sub_col<1>(params, col);
    bw_sub_col<0>(params, col);
}

template<typename BoundaryConditions, typename ForcingTerm>
__global__ void solve_sweep_Dxx(GRID_CONSTANT const SolveDxxParams params)
{
    const uint32_t y = blockDim.x * blockIdx.y + threadIdx.x;

    const uint32_t depth = params.depth;
    const uint32_t height = params.height;
    const uint32_t width = params.width;

    const uint64_t num_points = depth * height * width;

    const _ftype128 *__restrict__ p_k = P_128(params.k);
    const _ftype128 *__restrict__ p_w = P_128(params.w);
    const _ftype128 *__restrict__ p_p = P_128(params.p);
    const _ftype128 *__restrict__ p_phi = P_128(params.phi);

    ftype128_vec p_eta(params.eta_x, params.eta_y, params.eta_z);
    const ftype128_vec p_zeta(params.zeta_x, params.zeta_y, params.zeta_z);
    const ftype128_vec p_vel(params.u_x, params.u_y, params.u_z);

    ftype *__restrict__ tmp = params.tmp;
    ftype *__restrict__ tmp_f_x = params.tmp + num_points;
    ftype *__restrict__ tmp_f_y = params.tmp + num_points * 2;
    ftype *__restrict__ tmp_f_z = params.tmp + num_points * 3;

    const BoundaryConditions bc;

    FwReduParams fw_params;
    fw_params.tmp_stride = depth * height;

    for (uint32_t z = blockIdx.x; z < depth; z += gridDim.x) {

        // Traverse tmp buffers front to back, map yz to xy
        const uint64_t tmp_offset = height * z + y;

        fw_params.tmp_upp = tmp + tmp_offset;
        fw_params.tmp_f_x = tmp_f_x + tmp_offset;
        fw_params.tmp_f_y = tmp_f_y + tmp_offset;
        fw_params.tmp_f_z = tmp_f_z + tmp_offset;

        const uint64_t warp_tile_offset =
            (height * width * z + width * (y & ~31)) / VEC_LENGTH;

        ftype3 eta_prev = { 0, 0, 0 }; // Last col of prev eta tile

        // Loop peeling on first tile of the row to apply BCs

        fw_params.w_t.load(p_w + warp_tile_offset, width);
        fw_params.w_t.transpose();
        compute_Dxx_rhs_tile(p_k + warp_tile_offset,
                             p_p + warp_tile_offset,
                             p_phi + warp_tile_offset,
                             p_eta + warp_tile_offset,
                             p_zeta + warp_tile_offset,
                             p_vel + warp_tile_offset,
                             bc,
                             ForcingTerm(),
                             depth, height, width,
                             0, y & ~31, z,
                             params.timestep,
                             eta_prev,
                             fw_params.f_x_t,
                             fw_params.f_y_t,
                             fw_params.f_z_t);

        ftype3 u0 = bc.get_u_delta_left(0, y, z, params.timestep);

        fw_params.upp_redu = 0;
        fw_params.f_x_redu = u0.x;
        fw_params.f_y_redu = u0.y;
        fw_params.f_z_redu = u0.z;

        // TODO: Either a thread carries 4 rows concurrently, coarsening
        // in the y or z direction, or the stores to tmp have to
        // happen in chunks of 2/4 elements (to allow 128 bit wide
        // stores). In both cases we want vectorized stores.
        // I guess coarsening is only viable when using shared memory
        // to store tiles, since the number of alive tiles increases.

        fw_store_redu(fw_params);

        // TODO: buffer redu coefficients in _ftype128,
        // then issue vectorized stores

        fw_redu_col<1>(fw_params, 0);
#ifdef FLOAT
        fw_redu_col<2>(fw_params, 0);
        fw_redu_col<3>(fw_params, 0);
#endif
        // Consume remaining cols of tile
        #pragma unroll
        for (int tk = 1; tk < 8; ++tk) {
            fw_redu_cols(fw_params, tk);
        }

        const uint32_t tile_width = 8 * VEC_LENGTH;
        // Gauss reduce over the tiles of the row, except the last
        for (uint32_t k = tile_width;
                      k < width - tile_width; k += tile_width) {

            fw_params.w_t.load(
                p_w + warp_tile_offset + k / VEC_LENGTH, width);
            fw_params.w_t.transpose();
            compute_Dxx_rhs_tile(p_k + warp_tile_offset + k / VEC_LENGTH,
                                 p_p + warp_tile_offset + k / VEC_LENGTH,
                                 p_phi + warp_tile_offset + k / VEC_LENGTH,
                                 p_eta + warp_tile_offset + k / VEC_LENGTH,
                                 p_zeta + warp_tile_offset + k / VEC_LENGTH,
                                 p_vel + warp_tile_offset + k / VEC_LENGTH,
                                 bc,
                                 ForcingTerm(),
                                 depth, height, width,
                                 k, y & ~31, z,
                                 params.timestep,
                                 eta_prev,
                                 fw_params.f_x_t,
                                 fw_params.f_y_t,
                                 fw_params.f_z_t);
            // Consume tile
            #pragma unroll
            for (int tk = 0; tk < 8; ++tk) {
                fw_redu_cols(fw_params, tk);
            }
        }

        // Loop peeling on last tile to apply BCs and solve last col

        const uint64_t tile_offset =
            warp_tile_offset + (width - tile_width) / VEC_LENGTH;

        fw_params.w_t.load(p_w + tile_offset, width);
        fw_params.w_t.transpose();
        compute_Dxx_rhs_tile(p_k + tile_offset,
                             p_p + tile_offset,
                             p_phi + tile_offset,
                             p_eta + tile_offset,
                             p_zeta + tile_offset,
                             p_vel + tile_offset,
                             bc,
                             ForcingTerm(),
                             depth, height, width,
                             width - tile_width, y & ~31, z,
                             params.timestep,
                             eta_prev,
                             fw_params.f_x_t,
                             fw_params.f_y_t,
                             fw_params.f_z_t);

        #pragma unroll
        for (int tk = 0; tk < 7; ++tk) {
            fw_redu_cols(fw_params, tk);
        }
        fw_redu_col<0>(fw_params, 7);
#ifdef FLOAT
        fw_redu_col<1>(fw_params, 7);
        fw_redu_col<2>(fw_params, 7);
#endif

        // Directly solve last col by applying BCs

#ifdef FLOAT
        ftype w = fw_params.w_t.col<3>(7);
        ftype f_y = fw_params.f_y_t.col<3>(7);
        ftype f_z = fw_params.f_z_t.col<3>(7);
#else
        ftype w = fw_params.w_t.col<1>(7);
        ftype f_y = fw_params.f_y_t.col<1>(7);
        ftype f_z = fw_params.f_z_t.col<1>(7);
#endif

        ftype norm_coef_inv =
            1 / (1 + 4 * w + (ftype) (4.0 / 3) * w * fw_params.upp_redu);

        ftype3 un = bc.get_u_delta_right(width - 1, y, z, params.timestep);

        ftype un_y = (w * (ftype) (8.0 / 3.0) * un.y + f_y +
            (ftype) (4.0 / 3.0) * w * fw_params.f_y_redu) * norm_coef_inv;

        ftype un_z = (w * (ftype) (8.0 / 3.0) * un.z + f_z +
            (ftype) (4.0 / 3.0) * w * fw_params.f_z_redu) * norm_coef_inv;

        BwSubParams bw_params;
        bw_params.tmp_stride = depth * height;

        uint64_t bw_offset = tmp_offset + bw_params.tmp_stride * (width - 2);
        bw_params.tmp_upp = tmp + bw_offset;
        bw_params.tmp_f_x = tmp_f_x + bw_offset;
        bw_params.tmp_f_y = tmp_f_y + bw_offset;
        bw_params.tmp_f_z = tmp_f_z + bw_offset;

        // note that when storing back u, you need to store-add
        // The backward substitution starts from the last column's
        // solution, which is the BC for x and the directly solved
        // tangential components.
        bw_params.u_x_prev = un.x;
        bw_params.u_y_prev = un_y;
        bw_params.u_z_prev = un_z;

#ifdef FLOAT
        bw_params.u_x_t.col<3>(7) = un.x;
        bw_params.u_y_t.col<3>(7) = un_y;
        bw_params.u_z_t.col<3>(7) = un_z;

        bw_sub_col<2>(bw_params, 7);
        bw_sub_col<1>(bw_params, 7);
#else
        bw_params.u_x_t.col<1>(7) = un.x;
        bw_params.u_y_t.col<1>(7) = un_y;
        bw_params.u_z_t.col<1>(7) = un_z;
#endif
        bw_sub_col<0>(bw_params, 7);

        // Bw substitution over remaining tile
        #pragma unroll
        for (int tk = 1; tk < 8; ++tk) {
            bw_sub_cols(bw_params, 7 - tk);
        }

        bw_params.u_x_t.transpose();
        bw_params.u_y_t.transpose();
        bw_params.u_z_t.transpose();

        const uint64_t dst_offset =
            warp_tile_offset + width / VEC_LENGTH - 8;

        bw_params.u_x_t.store_add(p_eta.x_cc() + dst_offset, width);
        bw_params.u_y_t.store_add(p_eta.y_cc() + dst_offset, width);
        bw_params.u_z_t.store_add(p_eta.z_cc() + dst_offset, width);

        for (uint32_t k = tile_width;
                      k < width; k += tile_width) {

            #pragma unroll
            for (int tk = 0; tk < 8; ++tk) {
                bw_sub_cols(bw_params, 7 - tk);
            }

            const uint64_t tile_offset =
                warp_tile_offset + (width - k - tile_width) / VEC_LENGTH;

            bw_params.u_x_t.transpose();
            bw_params.u_y_t.transpose();
            bw_params.u_z_t.transpose();

            bw_params.u_x_t.store_add(p_eta.x_cc() + tile_offset, width);
            bw_params.u_y_t.store_add(p_eta.y_cc() + tile_offset, width);
            bw_params.u_z_t.store_add(p_eta.z_cc() + tile_offset, width);
        }
    }
}

// TODO: restrict could be avoided here and kept only in the kernel signature
void launch_momentum_solve_Dxx(const ftype *__restrict__ k,
                               const ftype *__restrict__ w,
                               const ftype *__restrict__ p,
                               const ftype *__restrict__ phi,
                               uint32_t depth,
                               uint32_t height,
                               uint32_t width,
                               uint32_t timestep,
                               ftype *__restrict__ tmp,
                               ftype *__restrict__ eta_x,
                               ftype *__restrict__ eta_y,
                               ftype *__restrict__ eta_z,
                               const ftype *__restrict__ zeta_x,
                               const ftype *__restrict__ zeta_y,
                               const ftype *__restrict__ zeta_z,
                               const ftype *__restrict__ u_x,
                               const ftype *__restrict__ u_y,
                               const ftype *__restrict__ u_z)
{
    // Perform h2d transfers
    uint64_t alloc_size = depth * height * width * sizeof(ftype);

    // NOTE: The first two transfers are redundant.
    CUDA_CHECK(cudaMemcpy(d_data_dxx.k, k, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.w, w, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.p, p, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.phi, phi, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.eta_x, eta_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.eta_y, eta_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.eta_z, eta_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.zeta_x, zeta_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.zeta_y, zeta_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.zeta_z, zeta_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.u_x, u_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.u_y, u_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dxx.u_z, u_z, alloc_size, CUDA_H2D));

    // Spawn one thread per cell in the yz plane.

    const int num_threads = 64;
    const dim3 num_blocks(depth / COARSE_FACTOR_Y, height / num_threads);

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;
    using F = ForcingTerm<ManFx, ManFy, ManFz>;

    d_data_dxx.timestep = timestep;

    CUDA_TIMER_CREATE(solve_Dxx_sweep);
    CUDA_TIMER_RESTART(solve_Dxx_sweep);
    CUDA_CHECK_LAUNCH(
        solve_sweep_Dxx<BC, F><<<num_blocks, num_threads>>>(d_data_dxx));
    cudaDeviceSynchronize();
    CUDA_TIMER_ELAPSED(solve_Dxx_sweep, 1);
    CUDA_TIMER_DESTROY(solve_Dxx_sweep);

    CUDA_CHECK(cudaMemcpy(eta_x, d_data_dxx.eta_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(eta_y, d_data_dxx.eta_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(eta_z, d_data_dxx.eta_z, alloc_size, CUDA_D2H));
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
    CUDA_CHECK(cudaMemcpy(d_data_dzz.w, w, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_x, f_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_y, f_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_z, f_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_x, u_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_y, u_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_z, u_z, alloc_size, CUDA_H2D));

    // Spawn one thread per cell in the xz plane.

    const int num_threads = 64;
    // NOTE: There's no point in rounding, the domain is expected
    // to have sizes multiple of 32
    const dim3 num_blocks((width / VEC_LENGTH - 1) / num_threads + 1,
                          depth / COARSE_FACTOR_Y);

    d_data_dzz.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    CUDA_TIMER_CREATE(solve_Dyy_sweep);
    CUDA_TIMER_RESTART(solve_Dyy_sweep);
    CUDA_CHECK_LAUNCH(
        solve_sweep<BC, DyySweep><<<num_blocks, num_threads>>>(d_data_dzz));
    cudaDeviceSynchronize();
    CUDA_TIMER_ELAPSED(solve_Dyy_sweep, 1);
    CUDA_TIMER_DESTROY(solve_Dyy_sweep);

    CUDA_CHECK(cudaMemcpy(u_x, d_data_dzz.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data_dzz.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data_dzz.u_z, alloc_size, CUDA_D2H));

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
    CUDA_CHECK(cudaMemcpy(d_data_dzz.w, w, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_x, f_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_y, f_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.f_z, f_z, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_x, u_x, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_y, u_y, alloc_size, CUDA_H2D));
    CUDA_CHECK(cudaMemcpy(d_data_dzz.u_z, u_z, alloc_size, CUDA_H2D));

    // Spawn one thread per point in the xy plane.

    const int num_threads = 64;
    const dim3 num_blocks((width / VEC_LENGTH - 1) / num_threads + 1,
                          height / COARSE_FACTOR_Y);

    d_data_dzz.timestep = timestep;

    using BC = BoundaryConditions<ManUx, ManUy, ManUz>;

    // If the domain is small enough in the depth dimension, shared
    // memory could be used to hold the reduced coefficients.

    CUDA_TIMER_CREATE(solve_Dzz_sweep);
    CUDA_TIMER_RESTART(solve_Dzz_sweep);
    CUDA_CHECK_LAUNCH(
        solve_sweep<BC, DzzSweep><<<num_blocks, num_threads>>>(d_data_dzz));
    cudaDeviceSynchronize();
    CUDA_TIMER_ELAPSED(solve_Dzz_sweep, 1);
    CUDA_TIMER_DESTROY(solve_Dzz_sweep);

    CUDA_CHECK(cudaMemcpy(u_x, d_data_dzz.u_x, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_y, d_data_dzz.u_y, alloc_size, CUDA_D2H));
    CUDA_CHECK(cudaMemcpy(u_z, d_data_dzz.u_z, alloc_size, CUDA_D2H));

    // No need to synchronize
}
