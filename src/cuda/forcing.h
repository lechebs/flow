#ifndef CUDA_FORCING_H
#define CUDA_FORCING_H

#include "ftype.h"
#include "consts.h"

#include <math.h>

template<typename FuncX, typename FuncY, typename FuncZ>
struct ForcingTerm
{
    __device__ ForcingTerm() {}

    __device__ __forceinline__ ftype x(uint32_t x,
                                       uint32_t y,
                                       uint32_t z,
                                       uint32_t t) const
    {
        return f_x(x * _D_DX + _D_DX / 2, y * _D_DX, z * _D_DX,
                   t * _D_DT - _D_DT / 2);
    }

    __device__ __forceinline__ ftype y(uint32_t x,
                                       uint32_t y,
                                       uint32_t z,
                                       uint32_t t) const
    {
        return f_y(x * _D_DX, y * _D_DX + _D_DX / 2, z * _D_DX,
                   t * _D_DT - _D_DT / 2);
    }

    __device__ __forceinline__ ftype z(uint32_t x,
                                       uint32_t y,
                                       uint32_t z,
                                       uint32_t t) const
    {
        return f_z(x * _D_DX, y * _D_DX, z * _D_DX + _D_DX / 2,
                   t * _D_DT - _D_DT / 2);
    }

private:
    const FuncX f_x;
    const FuncY f_y;
    const FuncZ f_z;
};

// Sample forcing for manufactured solution

__device__ __forceinline__
ftype man_porosity(ftype x, ftype y, ftype z, ftype t)
{
    return 2 + sin(M_PI * x) * cos(M_PI * y + t) * sin(M_PI * z);
}

struct ManFx
{
    __device__ ManFx() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        ftype k = man_porosity(x - _D_DX / 2, y, z, t);
        return sin(x) * (-sin(y + t) * sin(z) +
                         cos(y + t) * sin(z) * _D_NU * (3 + 1.0 / k) +
                         cos(y + t) * cos(z) * -3 * _D_NU);
    }
};

struct ManFy
{
    __device__ ManFy() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        ftype k = man_porosity(x, y - _D_DX / 2, z, t);
        return cos(x) * (cos(y + t) * sin(z) +
                         sin(y + t) * sin(z) * _D_NU * (3 + 1.0 / k) +
                         sin(y + t) * cos(z) * -3 * _D_NU);
    }
};

struct ManFz
{
    __device__ ManFz() {}

    __device__ __forceinline__
    ftype operator()(ftype x, ftype y, ftype z, ftype t) const
    {
        ftype k = man_porosity(x, y, z - _D_DX / 2, t);
        return cos(x) * (sin(y + t) * cos(z) * -2 +
                         cos(y + t) * cos(z) * 2 * _D_NU * (3 + 1.0 / k) +
                         cos(y + t) * sin(z) * -3 * _D_NU);
    }
};

#endif
