#ifndef CUDA_BOUNDARY_H
#define CUDA_BOUNDARY_H

#include "ftype.h"

// TODO: Move to separate header
__constant__ ftype _D_DX;
__constant__ ftype _D_DT;

// Dirichlet boundary conditions for velocity
template<typename FuncX, typename FuncY, typename FuncZ>
struct BoundaryConditions
{
    __device__ BoundaryConditions() {}

    __device__ __forceinline__ ftype3 get_u_front(uint32_t x,
                                                  uint32_t y,
                                                  uint32_t z,
                                                  uint32_t t) const
    {
        ftype3 u;

        u.x = u_x(x * _D_DX + _D_DX / 2, y * _D_DX, 0, t * _D_DT);
        u.y = u_y(x * _D_DX, y * _D_DX + _D_DX / 2, 0, t * _D_DT);

        ftype dux_dx = u_x(x * _D_DX + _D_DX / 2, y * _D_DX, 0, t * _D_DT) -
                       u_x(x * _D_DX - _D_DX / 2, y * _D_DX, 0, t * _D_DT);

        ftype duy_dy = u_y(x * _D_DX, y * _D_DX + _D_DX / 2, 0, t * _D_DT) -
                       u_y(x * _D_DX, y * _D_DX - _D_DX / 2, 0, t * _D_DT);

        u.z = u_z(x * _D_DX, y * _D_DX, 0, t * _D_DT) - (dux_dx + duy_dy) / 2;

        return u;
    }

    __device__ __forceinline__ ftype3 get_u_back(uint32_t x,
                                                 uint32_t y,
                                                 uint32_t z,
                                                 uint32_t t) const
    {
        ftype3 u;

        u.x = u_x(
            x * _D_DX + _D_DX / 2, y * _D_DX, z * _D_DX + _D_DX / 2, t * _D_DT);
        u.y = u_y(
            x * _D_DX, y * _D_DX + _D_DX / 2, z * _D_DX + _D_DX / 2, t * _D_DT);
        u.z = u_z(
            x * _D_DX, y * _D_DX, z * _D_DX + _D_DX / 2, t * _D_DT);

        return u;
    }

    __device__ __forceinline__ ftype3 get_u_delta_front(uint32_t x,
                                                        uint32_t y,
                                                        uint32_t z,
                                                        uint32_t t) const
    {
        ftype3 u_prev = get_u_front(x, y, z, t - 1);
        ftype3 u_curr = get_u_front(x, y, z, t);

        return { u_curr.x - u_prev.x,
                 u_curr.y - u_prev.y,
                 u_curr.z - u_prev.z };
    }

    __device__ __forceinline__ ftype3 get_u_delta_back(uint32_t x,
                                                       uint32_t y,
                                                       uint32_t z,
                                                       uint32_t t) const
    {
        ftype3 u_prev = get_u_back(x, y, z, t - 1);
        ftype3 u_curr = get_u_back(x, y, z, t);

        return { u_curr.x - u_prev.x,
                 u_curr.y - u_prev.y,
                 u_curr.z - u_prev.z };
    }

private:
    const FuncX u_x;
    const FuncY u_y;
    const FuncZ u_z;
};

// Sample manufactured solution

struct ManUx
{
    __device__ ManUx() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        return sin(x) * cos(y + t) * sin(z);
    }
};

struct ManUy
{
    __device__ ManUy() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        return cos(x) * sin(y + t) * sin(z);
    }
};

struct ManUz
{
    __device__ ManUz() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        return 2 * cos(x) * cos(y + t) * cos(z);
    }
};


#endif
