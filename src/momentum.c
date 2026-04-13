#include <math.h>
#include <assert.h>

#include "momentum.h"

#include "ftype.h"
#include "alloc.h"
#include "field.h"
#include "finite-diff.h"
#include "boundary.h"
#include "consts.h"
#include "timeit.h"
#include "thread-array.h"
#include "ddecomp.h"

DECLARE_FORCING()

DECLARE_BC_U(BC_LEFT)
DECLARE_BC_U(BC_TOP)
DECLARE_BC_U(BC_FRONT)
DECLARE_BC_U(BC_RIGHT)
DECLARE_BC_U(BC_BOTTOM)
DECLARE_BC_U(BC_BACK)

#ifdef AUTO_VEC
    #define vneg(v) (-(v))
#else
    #define vneg(vec) vxor(vec, SIGN_MASK)
#endif

/* TODO: Remove these. */
static vftype ZEROS;
static vftype ONES;
static vftype SIGN_MASK;

static inline __attribute__((always_inline))
void apply_start_bc(vftype u0_x,
                    vftype u0_y,
                    vftype u0_z,
                    ftype *restrict upper,
                    ftype *restrict f_x,
                    ftype *restrict f_y,
                    ftype *restrict f_z)
{
    /* u_x = u_0 + (-du_y/dy -du_z/dz) * dx/2
     * u_y = u0_y
     * u_z = u0_z */

    /* Set upper coefficient to 0 and enforce solution in rhs. */
    vstore(upper, ZEROS);
    /* du_y/dy + du_z/dz is already handled in the bc getter. */
    vstore(f_x, u0_x);
    vstore(f_y, u0_y);
    vstore(f_z, u0_z);
}

static inline __attribute__((always_inline))
vftype compute_end_bc_tang_u(vftype ws,
                             vftype ws2,
                             vftype fs_prev,
                             vftype fs,
                             vftype uns,
                             vftype norm_coeffs)
{
    //return vdiv(vadd(vfmadd(fs_prev, ws, fs),
    //                 vmul(ws2, uns)),
    //            norm_coeffs);

    return (ws * 8.0 / 3.0 * uns +
            fs + (ftype) (4.0 / 3.0) * ws * fs_prev) / norm_coeffs;
}

static inline __attribute__((always_inline))
void apply_end_bc(const ftype *restrict w,
                  const ftype *restrict upper_prev,
                  const ftype *restrict f_y_prev,
                  const ftype *restrict f_z_prev,
                  const ftype *restrict f_y,
                  const ftype *restrict f_z,
                  const ftype *restrict u_y,
                  const ftype *restrict u_z,
                  vftype un_x,
                  vftype un_y,
                  vftype un_z,
                  vftype *restrict us_x,
                  vftype *restrict us_y,
                  vftype *restrict us_z)
{
    /* u_x = un_x
     * u_y =  (2w_i un_y + f_y_i + w_i f_y_i_prev) /
     *        (1 + 3w_i + w_i upper_prev_i)
     * u_z =  (2w_i un_z + f_z_i + w_i f_z_i_prev) /
     *        (1 + 3w_i + w_i upper_prev_i) */

    vftype ws = vload(w);
    vftype upper_prevs = vload(upper_prev);
    vftype fs_y_prevs = vload(f_y_prev);
    vftype fs_z_prevs = vload(f_z_prev);
    vftype fs_y = vsub(vload(f_y), vload(u_y));
    vftype fs_z = vsub(vload(f_z), vload(u_z));
    vftype ws2 = vadd(ws, ws);
    /*
    vftype norm_coeffs = vfmadd(upper_prevs, ws,
                                vadd(ONES, vadd(ws2, ws)));
    */
    vftype norm_coeffs = 1 + 4 * ws + (ftype) (4.0 / 3.0) * ws * upper_prevs;


    *us_x = un_x;
    *us_y = compute_end_bc_tang_u(ws, ws2, fs_y_prevs,
                                  fs_y, un_y, norm_coeffs);
    *us_z = compute_end_bc_tang_u(ws, ws2, fs_z_prevs,
                                  fs_z, un_z, norm_coeffs);
}

static inline __attribute__((always_inline))
void apply_left_bc(uint32_t x,
                   uint32_t y,
                   uint32_t z,
                   uint32_t t,
                   ftype *restrict upper,
                   ftype *restrict f_x,
                   ftype *restrict f_y,
                   ftype *restrict f_z)
{
    vftype u0_x, u0_y, u0_z;
    _get_left_bc_u_delta(x, y, z, t, &u0_x, &u0_y, &u0_z);

    apply_start_bc(u0_x, u0_y, u0_z, upper, f_x, f_y, f_z);
}

static inline __attribute__((always_inline))
void apply_right_bc(const ftype *restrict w,
                    const ftype *restrict upper_prev,
                    const ftype *restrict f_y_prev,
                    const ftype *restrict f_z_prev,
                    uint32_t x,
                    uint32_t y,
                    uint32_t z,
                    uint32_t t,
                    ftype *restrict f_x,
                    ftype *restrict f_y,
                    ftype *restrict f_z,
                    vftype *restrict u_x,
                    vftype *restrict u_y,
                    vftype *restrict u_z)
{
    vftype un_x, un_y, un_z;
    _get_right_bc_u_delta(x, y, z, t, &un_x, &un_y, &un_z);

    /* Quick hack :/ */
    ftype __attribute__((aligned(32))) zeros[VLEN * 2] = {0};

    apply_end_bc(w, upper_prev, f_y_prev, f_z_prev,
                 f_y, f_z, zeros, zeros + VLEN, un_x, un_y,
                 un_z, u_x, u_y, u_z);

    vstore(f_x, *u_x);
    vstore(f_y, *u_y);
    vstore(f_z, *u_z);
}

static inline __attribute__((always_inline))
void gauss_reduce_vstrip(const ftype *restrict w,
                         ftype *restrict upper_prev,
                         const ftype *restrict f_x_src,
                         const ftype *restrict f_y_src,
                         const ftype *restrict f_z_src,
                         ftype *restrict f_x_dst,
                         ftype *restrict f_y_dst,
                         ftype *restrict f_z_dst)
{
    vftype ws = vload(w);
    vftype upper_prevs = vload(upper_prev);
    vftype f_x_prevs = vload(f_x_dst - VLEN);
    vftype f_y_prevs = vload(f_y_dst - VLEN);
    vftype f_z_prevs = vload(f_z_dst - VLEN);
    vftype fs_x = vload(f_x_src);
    vftype fs_y = vload(f_y_src);
    vftype fs_z = vload(f_z_src);
    vftype norm_coefs = vfmadd(ws, upper_prevs, vadd(ONES, vadd(ws, ws)));
    vstore(upper_prev + VLEN, vdiv(vneg(ws), norm_coefs));
    vstore(f_x_dst, vdiv(vfmadd(ws, f_x_prevs, fs_x), norm_coefs));
    vstore(f_y_dst, vdiv(vfmadd(ws, f_y_prevs, fs_y), norm_coefs));
    vstore(f_z_dst, vdiv(vfmadd(ws, f_z_prevs, fs_z), norm_coefs));
}

static inline __attribute__((always_inline))
void backward_sub_vstrip(const ftype *restrict f_x,
                         const ftype *restrict f_y,
                         const ftype *restrict f_z,
                         const ftype *restrict upper,
                         vftype *restrict u_x_prevs,
                         vftype *restrict u_y_prevs,
                         vftype *restrict u_z_prevs,
                         ftype *restrict u_x,
                         ftype *restrict u_y,
                         ftype *restrict u_z)
{
    vftype fs_x = vload(f_x);
    vftype fs_y = vload(f_y);
    vftype fs_z = vload(f_z);
    vftype uppers = vload(upper);
    *u_x_prevs = vfmadd(vneg(uppers), *u_x_prevs, fs_x);
    vstore(u_x, *u_x_prevs);
    *u_y_prevs = vfmadd(vneg(uppers), *u_y_prevs, fs_y);
    vstore(u_y, *u_y_prevs);
    *u_z_prevs = vfmadd(vneg(uppers), *u_z_prevs, fs_z);
    vstore(u_z, *u_z_prevs);
}

#define INNER_COL  0
#define INNER_ROW  0
#define INNER_FACE 0
#define LAST_COL   1
#define LAST_ROW   1
#define LAST_FACE  1

#ifdef FLOAT

#define def_vtile(name) \
    vftype name##0, name##1, name##2, name##3, \
           name##4, name##5, name##6, name##7

#define load_vtile(src, stride, tile)   \
do {                                    \
    tile##0 = vload(src + 0 * stride);  \
    tile##1 = vload(src + 1 * stride);  \
    tile##2 = vload(src + 2 * stride);  \
    tile##3 = vload(src + 3 * stride);  \
    tile##4 = vload(src + 4 * stride);  \
    tile##5 = vload(src + 5 * stride);  \
    tile##6 = vload(src + 6 * stride);  \
    tile##7 = vload(src + 7 * stride);  \
} while (0)

#define load_vtile_add(src, stride, tile)  \
do {                                       \
    tile##0 += vload(src + 0 * stride);    \
    tile##1 += vload(src + 1 * stride);    \
    tile##2 += vload(src + 2 * stride);    \
    tile##3 += vload(src + 3 * stride);    \
    tile##4 += vload(src + 4 * stride);    \
    tile##5 += vload(src + 5 * stride);    \
    tile##6 += vload(src + 6 * stride);    \
    tile##7 += vload(src + 7 * stride);    \
} while (0)

#define store_vtile(dst, stride, tile)  \
do {                                    \
    vstore(dst + 0 * stride, tile##0);  \
    vstore(dst + 1 * stride, tile##1);  \
    vstore(dst + 2 * stride, tile##2);  \
    vstore(dst + 3 * stride, tile##3);  \
    vstore(dst + 4 * stride, tile##4);  \
    vstore(dst + 5 * stride, tile##5);  \
    vstore(dst + 6 * stride, tile##6);  \
    vstore(dst + 7 * stride, tile##7);  \
} while (0)

#define fin_diff_vtiles(src_tile, src_tile_prev, dst_tile)  \
do {                                                        \
    dst_tile##0 = (src_tile##0 - src_tile_prev##0) / -_DX;  \
    dst_tile##1 = (src_tile##1 - src_tile_prev##1) / -_DX;  \
    dst_tile##2 = (src_tile##2 - src_tile_prev##2) / -_DX;  \
    dst_tile##3 = (src_tile##3 - src_tile_prev##3) / -_DX;  \
    dst_tile##4 = (src_tile##4 - src_tile_prev##4) / -_DX;  \
    dst_tile##5 = (src_tile##5 - src_tile_prev##5) / -_DX;  \
    dst_tile##6 = (src_tile##6 - src_tile_prev##6) / -_DX;  \
    dst_tile##7 = (src_tile##7 - src_tile_prev##7) / -_DX;  \
} while (0)

#define fin_diff_vtile(src_tile, dst_tile)             \
do {                                                   \
    dst_tile##0 = (src_tile##1 - src_tile##0) / -_DX;  \
    dst_tile##1 = (src_tile##2 - src_tile##1) / -_DX;  \
    dst_tile##2 = (src_tile##3 - src_tile##2) / -_DX;  \
    dst_tile##3 = (src_tile##4 - src_tile##3) / -_DX;  \
    dst_tile##4 = (src_tile##5 - src_tile##4) / -_DX;  \
    dst_tile##5 = (src_tile##6 - src_tile##5) / -_DX;  \
    dst_tile##6 = (src_tile##7 - src_tile##6) / -_DX;  \
    dst_tile##7 = (src_tile##8 - src_tile##7) / -_DX;  \
} while (0)

#define transpose_vtile_ip(tile)                        \
    vtranspose(&tile##0, &tile##1, &tile##2, &tile##3,  \
               &tile##4, &tile##5, &tile##6, &tile##7)

#else

#define def_vtile(name) \
    vftype name##0, name##1, name##2, name##3;


#define load_vtile(src, stride, tile)   \
do {                                    \
    tile##0 = vload(src + 0 * stride);  \
    tile##1 = vload(src + 1 * stride);  \
    tile##2 = vload(src + 2 * stride);  \
    tile##3 = vload(src + 3 * stride);  \
} while (0)

#define load_vtile_add(src, stride, tile)  \
do {                                       \
    tile##0 += vload(src + 0 * stride);    \
    tile##1 += vload(src + 1 * stride);    \
    tile##2 += vload(src + 2 * stride);    \
    tile##3 += vload(src + 3 * stride);    \
} while (0)

#define store_vtile(dst, stride, tile)  \
do {                                    \
    vstore(dst + 0 * stride, tile##0);  \
    vstore(dst + 1 * stride, tile##1);  \
    vstore(dst + 2 * stride, tile##2);  \
    vstore(dst + 3 * stride, tile##3);  \
} while (0)

#define fin_diff_vtiles(src_tile, src_tile_prev, dst_tile)  \
do {                                                        \
    dst_tile##0 = (src_tile##0 - src_tile_prev##0) / -_DX;  \
    dst_tile##1 = (src_tile##1 - src_tile_prev##1) / -_DX;  \
    dst_tile##2 = (src_tile##2 - src_tile_prev##2) / -_DX;  \
    dst_tile##3 = (src_tile##3 - src_tile_prev##3) / -_DX;  \
} while (0)

#define fin_diff_vtile(src_tile, dst_tile)             \
do {                                                   \
    dst_tile##0 = (src_tile##1 - src_tile##0) / -_DX;  \
    dst_tile##1 = (src_tile##2 - src_tile##1) / -_DX;  \
    dst_tile##2 = (src_tile##3 - src_tile##2) / -_DX;  \
    dst_tile##3 = (src_tile##4 - src_tile##3) / -_DX;  \
} while (0)

#define transpose_vtile_ip(tile) \
    vtranspose(&tile##0, &tile##1, &tile##2, &tile##3);

#endif

#define COMPONENT_X 0
#define COMPONENT_Y 1
#define COMPONENT_Z 2

static inline __attribute__((always_inline))
void compute_Dxx_Dyy_Dzz(const ftype *restrict eta,
                         const ftype *restrict zeta,
                         const ftype *restrict vel,
                         int is_last_face,
                         int is_last_row,
                         int is_last_col,
                         int component,
                         uint32_t i,
                         uint32_t j,
                         uint32_t k,
                         uint32_t depth,
                         uint32_t height,
                         uint32_t width,
                         uint32_t timestep,
                         ftype *restrict rhs)
{
    /* NOTE: All of the conditionals here will be
     * evaluated at compile time and optimized away. */

    def_vtile(dd);
    load_vtile(vel, width, dd);

    /* Compute Dzz vel */

    if (is_last_face && component == COMPONENT_Z) {
         /* We don't care about the value of the rhs here. */
        return;

    } else if (is_last_face) {

        def_vtile(v);
        /* WARNING: Safe but vel must be padded. */
        load_vtile(vel - height * width, width, v);

        dd0 = (ftype) 4.0 / 3 * v0 - 4 * dd0;
        dd1 = (ftype) 4.0 / 3 * v1 - 4 * dd1;
        dd2 = (ftype) 4.0 / 3 * v2 - 4 * dd2;
        dd3 = (ftype) 4.0 / 3 * v3 - 4 * dd3;
#ifdef FLOAT
        dd4 = (ftype) 4.0 / 3 * v4 - 4 * dd4;
        dd5 = (ftype) 4.0 / 3 * v5 - 4 * dd5;
        dd6 = (ftype) 4.0 / 3 * v6 - 4 * dd6;
        dd7 = (ftype) 4.0 / 3 * v7 - 4 * dd7;
#endif

        if (component == COMPONENT_X) {
            vftype _y, _z; /* Dummy variables, they will be optmized away. */
            _get_back_bc_u(k, j + 0, i, timestep - 1, &v0, &_y, &_z);
            _get_back_bc_u(k, j + 1, i, timestep - 1, &v1, &_y, &_z);
            _get_back_bc_u(k, j + 2, i, timestep - 1, &v2, &_y, &_z);
            _get_back_bc_u(k, j + 3, i, timestep - 1, &v3, &_y, &_z);
#ifdef FLOAT
            _get_back_bc_u(k, j + 0, i, timestep - 1, &v4, &_y, &_z);
            _get_back_bc_u(k, j + 1, i, timestep - 1, &v5, &_y, &_z);
            _get_back_bc_u(k, j + 2, i, timestep - 1, &v6, &_y, &_z);
            _get_back_bc_u(k, j + 3, i, timestep - 1, &v7, &_y, &_z);
#endif
        } else { /* component == COMPONENT_Y */
            vftype _x, _z; /* Dummy variables, they will be optmized away. */
            _get_back_bc_u(k, j + 0, i, timestep - 1, &_x, &v0, &_z);
            _get_back_bc_u(k, j + 1, i, timestep - 1, &_x, &v1, &_z);
            _get_back_bc_u(k, j + 2, i, timestep - 1, &_x, &v2, &_z);
            _get_back_bc_u(k, j + 3, i, timestep - 1, &_x, &v3, &_z);
#ifdef FLOAT
            _get_back_bc_u(k, j + 0, i, timestep - 1, &_x, &v4, &_z);
            _get_back_bc_u(k, j + 1, i, timestep - 1, &_x, &v5, &_z);
            _get_back_bc_u(k, j + 2, i, timestep - 1, &_x, &v6, &_z);
            _get_back_bc_u(k, j + 3, i, timestep - 1, &_x, &v7, &_z);
#endif
        }

        dd0 += (ftype) 8.0 / 3 * v0;
        dd1 += (ftype) 8.0 / 3 * v1;
        dd2 += (ftype) 8.0 / 3 * v2;
        dd3 += (ftype) 8.0 / 3 * v3;
#ifdef FLOAT
        dd4 += (ftype) 8.0 / 3 * v4;
        dd5 += (ftype) 8.0 / 3 * v5;
        dd6 += (ftype) 8.0 / 3 * v6;
        dd7 += (ftype) 8.0 / 3 * v7;
#endif

    } else {
        dd0 *= -2; dd1 *= -2; dd2 *= -2; dd3 *= -2;
#ifdef FLOAT
        dd4 *= -2; dd5 *= -2; dd6 *= -2; dd7 *= -2;
#endif
        load_vtile_add(vel - height * width, width, dd);
        load_vtile_add(vel + height * width, width, dd);
    }

    /* Compute Dyy zeta_x */
    def_vtile(zt);
    vftype ztp, ztn;

    load_vtile(zeta, width, zt);
    /* WARNING: Safe but zeta must be padded. */
    ztp = vload(zeta - width);

    dd0 += ztp - 2 * zt0 + zt1;
    dd1 += zt0 - 2 * zt1 + zt2;
    dd2 += zt1 - 2 * zt2 + zt3;
#ifdef FLOAT
    dd3 += zt2 - 2 * zt3 + zt4;
    dd4 += zt3 - 2 * zt4 + zt5;
    dd5 += zt4 - 2 * zt5 + zt6;
    dd6 += zt5 - 2 * zt6 + zt7;
#endif

    if (is_last_row && component != COMPONENT_Y) {
        if (component == COMPONENT_X) {
            vftype _y, _z;
            /* WARNING: Do not use global_j here, since it corresponds
             * to the tile global_j! Use global_j + VLEN */
            _get_bottom_bc_u(k, height - 1, i, timestep - 1, &ztn, &_y, &_z);
        } else { /* component == COMPONENT_Z */
            vftype _x, _y;
            _get_bottom_bc_u(k, height - 1, i, timestep - 1, &_x, &_y, &ztn);
        }
#ifdef FLOAT
        dd7 += (ftype) 4.0 / 3 * zt6 - 4 * zt7 + (ftype) 8.0 / 3 * ztn;
#else
        dd3 += (ftype) 4.0 / 3 * zt2 - 4 * zt3 + (ftype) 8.0 / 3 * ztn;
#endif
    } else if (!is_last_row) {
        ztn = vload(zeta + width * VLEN);
#ifdef FLOAT
        dd7 += zt6 - 2 * zt7 + ztn;
#else
        dd3 += zt2 - 2 * zt3 + ztn;
#endif
    }

     /* Compute Dxx eta_x */
    def_vtile(et);
    vftype etp, etn;
    /* WARNING: Safe but eta must be padded. */
    etp = vgather(eta - 1, width);
    load_vtile(eta, width, et);
    transpose_vtile_ip(et);

    etp = etp - 2 * et0 + et1;
    et0 = et0 - 2 * et1 + et2;
    et1 = et1 - 2 * et2 + et3;
#ifdef FLOAT
    et2 = et2 - 2 * et3 + et4;
    et3 = et3 - 2 * et4 + et5;
    et4 = et4 - 2 * et5 + et6;
    et5 = et5 - 2 * et6 + et7;
#endif

    if (is_last_col && component != COMPONENT_X) {
        if (component == COMPONENT_Y) {
            vftype _x, _z;
            _get_right_bc_u(width - 1, j, i, timestep - 1, &_x, &etn, &_z);
        } else { /* component == COMPONENT_Z */
             vftype _x, _y;
            _get_right_bc_u(width - 1, j, i, timestep - 1, &_x, &_y, &etn);
        }
#ifdef FLOAT
        et6 = (ftype) 4.0 / 3 * et6 - 4 * et7 + (ftype) 8.0 / 3 * etn;
#else
        et2 = (ftype) 4.0 / 3 * et2 - 4 * et3 + (ftype) 8.0 / 3 * etn;
#endif
    } else if (!is_last_col) {
        etn = vgather(eta + VLEN, width);
#ifdef FLOAT
        et6 = et6 - 2 * et7 + etn;
#else
        et2 = et2 - 2 * et3 + etn;
#endif
    }/* else {
        et3 = vbroadcast(0);
    }*/

#ifdef FLOAT
    vtranspose(&etp, &et0, &et1, &et2, &et3, &et4, &et5, &et6);
#else
    vtranspose(&etp, &et0, &et1, &et2);
#endif

    dd0 = (dd0 + etp) * _NU / (_DX * _DX);
    dd1 = (dd1 + et0) * _NU / (_DX * _DX);
    dd2 = (dd2 + et1) * _NU / (_DX * _DX);
    dd3 = (dd3 + et2) * _NU / (_DX * _DX);
#ifdef FLOAT
    dd4 = (dd4 + et3) * _NU / (_DX * _DX);
    dd5 = (dd5 + et4) * _NU / (_DX * _DX);
    dd6 = (dd6 + et5) * _NU / (_DX * _DX);
    dd7 = (dd7 + et6) * _NU / (_DX * _DX);
#endif

    load_vtile_add(rhs, VLEN, dd);
    store_vtile(rhs, VLEN, dd);
}

static inline __attribute__((always_inline))
void compute_rhs_vtile(const ftype *restrict porosity,
                       const ftype *restrict p,
                       const ftype *restrict phi,
                       const ftype *restrict eta_x,
                       const ftype *restrict eta_y,
                       const ftype *restrict eta_z,
                       const ftype *restrict zeta_x,
                       const ftype *restrict zeta_y,
                       const ftype *restrict zeta_z,
                       const ftype *restrict vel_x,
                       const ftype *restrict vel_y,
                       const ftype *restrict vel_z,
                       int is_last_face,
                       int is_last_row,
                       int is_last_col,
                       uint32_t i,
                       uint32_t j,
                       uint32_t k,
                       uint32_t depth,
                       uint32_t height,
                       uint32_t width,
                       uint32_t timestep,
                       ftype *restrict rhs_x_t,
                       ftype *restrict rhs_y_t,
                       ftype *restrict rhs_z_t)
{
    def_vtile(pp);
    /* Computing pressure predictor for current tile. */
    load_vtile(p, width, pp);
    load_vtile_add(phi, width, pp);

    if (!is_last_face) {
        def_vtile(ppz);
        load_vtile(p + height * width, width, ppz);
        load_vtile_add(phi + height * width, width, ppz);
        /* Computing pressure predictor z derivative. */
        fin_diff_vtiles(ppz, pp, ppz);
        store_vtile(rhs_z_t, VLEN, ppz);
    } else {
        vftype zeros = vbroadcast(0);
        vstore(rhs_z_t + 0 * VLEN, zeros);
        vstore(rhs_z_t + 1 * VLEN, zeros);
        vstore(rhs_z_t + 2 * VLEN, zeros);
        vstore(rhs_z_t + 3 * VLEN, zeros);
#ifdef FLOAT
        vstore(rhs_z_t + 4 * VLEN, zeros);
        vstore(rhs_z_t + 5 * VLEN, zeros);
        vstore(rhs_z_t + 6 * VLEN, zeros);
        vstore(rhs_z_t + 7 * VLEN, zeros);
#endif
    }

    /* Load additional row and compute dpp/dy. */
#ifdef FLOAT
    vftype pp8;
    if (is_last_row) {
        pp8 = pp7; /* Such that y derivative is zero. */
    } else {
        pp8 = vload(p + width * VLEN) + vload(phi + width * VLEN);
    }
#else
    vftype pp4;
    if (is_last_row) {
        pp4 = pp3; /* Such that y derivative is zero. */
    } else {
        pp4 = vload(p + width * VLEN) + vload(phi + width * VLEN);
    }
#endif
    def_vtile(ppy);
    fin_diff_vtile(pp, ppy);
    store_vtile(rhs_y_t, VLEN, ppy);

    /* NOTE: Why not using vloadu? */
    transpose_vtile_ip(pp);
    /* Compute dpp/dx. */
#ifdef FLOAT
    if (is_last_col) {
        pp8 = pp7; /* Such that x derivative is zero. */
    } else {
        pp8 = vgather(p + VLEN, width) + vgather(phi + VLEN, width);
    }
#else
    if (is_last_col) {
        pp4 = pp3; /* Such that x derivative is zero. */
    } else {
        pp4 = vgather(p + VLEN, width) + vgather(phi + VLEN, width);
    }
#endif
    fin_diff_vtile(pp, pp);
    transpose_vtile_ip(pp);
    store_vtile(rhs_x_t, VLEN, pp);

    compute_Dxx_Dyy_Dzz(eta_x, zeta_x, vel_x,
                        is_last_face, is_last_row, is_last_col, COMPONENT_X,
                        i, j, k, depth, height, width, timestep, rhs_x_t);

    compute_Dxx_Dyy_Dzz(eta_y, zeta_y, vel_y,
                        is_last_face, is_last_row, is_last_col, COMPONENT_Y,
                        i, j, k, depth, height, width, timestep, rhs_y_t);

    compute_Dxx_Dyy_Dzz(eta_z, zeta_z, vel_z,
                        is_last_face, is_last_row, is_last_col, COMPONENT_Z,
                        i, j, k, depth, height, width, timestep, rhs_z_t);

    /* eps_i = (1 + 2 nu dt / (2k_i + nu dt)) u_i + .. */

    for (int jj = 0; jj < VLEN; ++jj) {

        vftype k_ = vload(porosity + width * jj);
        vftype dt_over_beta = _DT / (1 + _DT * _NU / (2 * k_));
        vftype coeff = 1 - 2 * _DT * _NU / (2 * k_ + _DT * _NU);

        vftype rx = vload(rhs_x_t + VLEN * jj);
        vftype ry = vload(rhs_y_t + VLEN * jj);
        vftype rz = vload(rhs_z_t + VLEN * jj);

        vftype fx, fy, fz;
        _get_forcing(k, j + jj, i, timestep, &fx, &fy, &fz);

        rx = (fx + rx) * dt_over_beta;
        ry = (fy + ry) * dt_over_beta;
        rz = (fz + rz) * dt_over_beta;

        vftype vx = vload(vel_x + width * jj);
        vftype vy = vload(vel_y + width * jj);
        vftype vz = vload(vel_z + width * jj);

        rx += coeff * vx - vload(eta_x + width * jj);
        ry += coeff * vy - vload(eta_y + width * jj);
        rz += coeff * vz - vload(eta_z + width * jj);

        vstore(rhs_x_t + VLEN * jj, rx);
        vstore(rhs_y_t + VLEN * jj, ry);
        vstore(rhs_z_t + VLEN * jj, rz);
    }

    transpose_vtile_inplace(VLEN, rhs_x_t);
    transpose_vtile_inplace(VLEN, rhs_y_t);
    transpose_vtile_inplace(VLEN, rhs_z_t);
}

static inline __attribute__((always_inline))
void solve_vtile_row(const ftype *restrict porosity,
                     const ftype *restrict w,
                     const ftype *restrict p,
                     const ftype *restrict phi,
                     ftype *restrict eta_x,
                     ftype *restrict eta_y,
                     ftype *restrict eta_z,
                     const ftype *restrict zeta_x,
                     const ftype *restrict zeta_y,
                     const ftype *restrict zeta_z,
                     const ftype *restrict vel_x,
                     const ftype *restrict vel_y,
                     const ftype *restrict vel_z,
                     uint32_t global_i,
                     uint32_t j,
                     uint32_t depth,
                     uint32_t height,
                     uint32_t width,
                     uint32_t timestep,
                     int is_last_face,
                     int is_last_row,
                     ftype *restrict tmp)
{
    ftype *restrict tmp_upp = tmp;
    /* WARNING: Cache aliasing? */
    ftype *restrict tmp_f_x = tmp + 1 * width * VLEN;
    ftype *restrict tmp_f_y = tmp + 2 * width * VLEN;
    ftype *restrict tmp_f_z = tmp + 3 * width * VLEN;

    ftype __attribute__((aligned(32))) f_x_t[VLEN * VLEN];
    ftype __attribute__((aligned(32))) f_y_t[VLEN * VLEN];
    ftype __attribute__((aligned(32))) f_z_t[VLEN * VLEN];
    ftype __attribute__((aligned(32))) w_t[VLEN * VLEN];

    transpose_vtile(w, width, VLEN, w_t);

    compute_rhs_vtile(porosity, p, phi, eta_x, eta_y, eta_z,
                      zeta_x, zeta_y, zeta_z, vel_x, vel_y, vel_z,
                      is_last_face, is_last_row, INNER_COL,
                      global_i, j, 0, depth, height, width, timestep,
                      f_x_t, f_y_t, f_z_t);

    /* Apply BCs on the first column of the tile. */
    apply_left_bc(0, j, global_i, timestep,
                  tmp_upp, tmp_f_x, tmp_f_y, tmp_f_z);
    /* Reduce remaining columns of the tile. */
    for (int k = 1; k < VLEN; ++k) {
        gauss_reduce_vstrip(w_t + VLEN * k,
                            tmp_upp + VLEN * (k - 1),
                            f_x_t + VLEN * k,
                            f_y_t + VLEN * k,
                            f_z_t + VLEN * k,
                            tmp_f_x + VLEN * k,
                            tmp_f_y + VLEN * k,
                            tmp_f_z + VLEN * k);
    }

    /* Reduce remaining tiles except the last one. */
    for (uint32_t tk = VLEN; tk < width - VLEN; tk += VLEN) {
        /* Load and transpose next tile. */
        transpose_vtile(w + tk, width, VLEN, w_t);

        compute_rhs_vtile(porosity + tk, p + tk, phi + tk,
                          eta_x + tk, eta_y + tk, eta_z + tk,
                          zeta_x + tk, zeta_y + tk, zeta_z + tk,
                          vel_x + tk, vel_y + tk, vel_z + tk,
                          is_last_face, is_last_row, INNER_COL,
                          global_i, j, tk, depth, height, width, timestep,
                          f_x_t, f_y_t, f_z_t);
 
        for (uint32_t k = 0; k < VLEN; ++k) {
            gauss_reduce_vstrip(w_t + VLEN * k,
                                tmp_upp + VLEN * (tk + k - 1),
                                f_x_t + VLEN * k,
                                f_y_t + VLEN * k,
                                f_z_t + VLEN * k,
                                tmp_f_x + VLEN * (tk + k),
                                tmp_f_y + VLEN * (tk + k),
                                tmp_f_z + VLEN * (tk + k));
        }
    }

    transpose_vtile(w + width - VLEN, width, VLEN, w_t);

    compute_rhs_vtile(porosity + width - VLEN,
                      p + width - VLEN,
                      phi + width - VLEN,
                      eta_x + width - VLEN,
                      eta_y + width - VLEN,
                      eta_z + width - VLEN,
                      zeta_x + width - VLEN,
                      zeta_y + width - VLEN,
                      zeta_z + width - VLEN,
                      vel_x + width - VLEN,
                      vel_y + width - VLEN,
                      vel_z + width - VLEN,
                      is_last_face, is_last_row, LAST_COL,
                      global_i, j, width - VLEN,
                      depth, height, width, timestep,
                      f_x_t, f_y_t, f_z_t);

    /* Reduce last tile except last column. */
    for (int k = 0; k < VLEN - 1; ++k) {
        gauss_reduce_vstrip(w_t + VLEN * k,
                            tmp_upp + VLEN * (width - VLEN + k - 1),
                            f_x_t + VLEN * k,
                            f_y_t + VLEN * k,
                            f_z_t + VLEN * k,
                            tmp_f_x + VLEN * (width - VLEN + k),
                            tmp_f_y + VLEN * (width - VLEN + k),
                            tmp_f_z + VLEN * (width - VLEN + k));
    }

    vftype u_x_prev, u_y_prev, u_z_prev;
    /* Apply BCs on the right column. */
    apply_right_bc(w_t + VLEN * (VLEN - 1),
                   tmp_upp + VLEN * (width - 2),
                   tmp_f_y + VLEN * (width - 2),
                   tmp_f_z + VLEN * (width - 2),
                   width - 1, j, global_i,
                   timestep,
                   /* Write solutions into f_t buffers,
                    * we will reuse them for u_t buffers */
                   f_x_t + VLEN * (VLEN - 1),
                   f_y_t + VLEN * (VLEN - 1),
                   f_z_t + VLEN * (VLEN - 1),
                   &u_x_prev,
                   &u_y_prev,
                   &u_z_prev);

    /* Reuse local buffers. */
    ftype __attribute__((aligned(32))) *u_x_t = f_x_t;
    ftype __attribute__((aligned(32))) *u_y_t = f_y_t;
    ftype __attribute__((aligned(32))) *u_z_t = f_z_t;

    /* Backward substitute last tile (last col already solved). */
    for (int k = 1; k < VLEN; ++k) {
        backward_sub_vstrip(
            tmp_f_x + VLEN * (width - 1 - k),
            tmp_f_y + VLEN * (width - 1 - k),
            tmp_f_z + VLEN * (width - 1 - k),
            tmp_upp + VLEN * (width - 1 - k),
            &u_x_prev,
            &u_y_prev,
            &u_z_prev,
            u_x_t + VLEN * (VLEN - 1 - k),
            u_y_t + VLEN * (VLEN - 1 - k),
            u_z_t + VLEN * (VLEN - 1 - k));
    }
    transpose_vtile_add(u_x_t, VLEN, width, eta_x + width - VLEN);
    transpose_vtile_add(u_y_t, VLEN, width, eta_y + width - VLEN);
    transpose_vtile_add(u_z_t, VLEN, width, eta_z + width - VLEN);

    /* Backward substitute one tile at a time. */
    for (uint32_t tk = VLEN; tk < width; tk += VLEN) {
        for (int k = 0; k < VLEN; ++k) {
            backward_sub_vstrip(
                tmp_f_x + VLEN * (width - 1 - (tk + k)),
                tmp_f_y + VLEN * (width - 1 - (tk + k)),
                tmp_f_z + VLEN * (width - 1 - (tk + k)),
                tmp_upp + VLEN * (width - 1 - (tk + k)),
                &u_x_prev,
                &u_y_prev,
                &u_z_prev,
                u_x_t + VLEN * (VLEN - 1 - k),
                u_y_t + VLEN * (VLEN - 1 - k),
                u_z_t + VLEN * (VLEN - 1 - k));
        }
         /* Transpose and store. */
        transpose_vtile_add(u_x_t, VLEN, width, eta_x + width - VLEN - tk);
        transpose_vtile_add(u_y_t, VLEN, width, eta_y + width - VLEN - tk);
        transpose_vtile_add(u_z_t, VLEN, width, eta_z + width - VLEN - tk);
    }
}

static void solve_Dxx_blocks_fused_rhs(const ftype *restrict k,
                                       const ftype *restrict w,
                                       const ftype *restrict p,
                                       const ftype *restrict phi,
                                       ftype *restrict eta_x,
                                       ftype *restrict eta_y,
                                       ftype *restrict eta_z,
                                       const ftype *restrict zeta_x,
                                       const ftype *restrict zeta_y,
                                       const ftype *restrict zeta_z,
                                       const ftype *restrict vel_x,
                                       const ftype *restrict vel_y,
                                       const ftype *restrict vel_z,
                                       uint32_t depth,
                                       uint32_t height,
                                       uint32_t width,
                                       uint32_t timestep,
                                       ftype *restrict tmp,
                                       uint32_t t_id,
                                       uint32_t num_threads,
                                       int proc_rank,
                                       int num_procs)
{
    uint32_t block_depth = (depth - 1) / num_threads;

    uint32_t block_face_start = block_depth * t_id;
    uint32_t block_face_end = block_depth * (t_id + 1);

    uint32_t global_face_start = proc_rank * depth;

    /* TODO: I could try to improve cache reuse when computing the rhs
     * by tiling the sweep across the depth. */

    /* Solving each face of the domain, except the last. */
    for (uint32_t i = block_face_start; i < block_face_end; ++i) {
        for (uint32_t j = 0; j < height - VLEN; j += VLEN) {
            /* Solving each row of tiles, except the last. */
            uint64_t off = height * width * i + width * j;
            solve_vtile_row(k + off, w + off, p + off, phi + off,
                            eta_x + off, eta_y + off, eta_z + off,
                            zeta_x + off, zeta_y + off, zeta_z + off,
                            vel_x + off, vel_y + off, vel_z + off,
                            global_face_start + i, j,
                            depth, height, width, timestep,
                            INNER_FACE, INNER_ROW, tmp);
        }
        /* Solving each tile of the last row. */
        uint64_t off = height * width * i + width * (height - VLEN);
        solve_vtile_row(k + off, w + off, p + off, phi + off, 
                        eta_x + off, eta_y + off, eta_z + off,
                        zeta_x + off, zeta_y + off, zeta_z + off,
                        vel_x + off, vel_y + off, vel_z + off,
                        global_face_start + i, height - VLEN,
                        depth, height, width, timestep,
                        INNER_FACE, LAST_ROW, tmp);
    }

    uint32_t num_remainder_faces = depth - 1 -
                                   block_depth * num_threads;

    if (t_id == 0 && proc_rank == (num_procs - 1)) {
        /* Solving the last face of the domain. */
        for (uint32_t j = 0; j < height - VLEN; j += VLEN) {
            uint64_t off = height * width * (depth - 1) + width * j;
            solve_vtile_row(k + off, w + off, p + off, phi + off,
                            eta_x + off, eta_y + off, eta_z + off,
                            zeta_x + off, zeta_y + off, zeta_z + off,
                            vel_x + off, vel_y + off, vel_z + off,
                            global_face_start + depth - 1, j,
                            depth, height, width, timestep,
                            LAST_FACE, INNER_ROW, tmp);
        }
        uint64_t off = height * width * (depth - 1) + width * (height - VLEN);
        solve_vtile_row(k + off, w + off, p + off, phi + off,
                        eta_x + off, eta_y + off, eta_z + off,
                        zeta_x + off, zeta_y + off, zeta_z + off,
                        vel_x + off, vel_y + off, vel_z + off,
                        global_face_start + depth - 1, height - VLEN,
                        depth, height, width, timestep,
                        LAST_FACE, LAST_ROW, tmp);

    } else if (t_id < num_remainder_faces + 1) {
        /* Solve remainder faces except the last. */
        for (uint32_t j = 0; j < height - VLEN; j += VLEN) {
            /* Solving each row of tiles, except the last. */
            uint64_t off = height * width * (depth - (t_id + 1)) + width * j;
            solve_vtile_row(k + off, w + off, p + off, phi + off,
                            eta_x + off, eta_y + off, eta_z + off,
                            zeta_x + off, zeta_y + off, zeta_z + off,
                            vel_x + off, vel_y + off, vel_z + off,
                            global_face_start + depth - (t_id + 1), j,
                            depth, height, width, timestep,
                            INNER_FACE, INNER_ROW, tmp);
        }
        /* Solving each tile of the last row. */
        uint64_t off = height * width * (depth - (t_id + 1)) +
                       width * (height - VLEN);
        solve_vtile_row(k + off, w + off, p + off, phi + off, 
                        eta_x + off, eta_y + off, eta_z + off,
                        zeta_x + off, zeta_y + off, zeta_z + off,
                        vel_x + off, vel_y + off, vel_z + off,
                        global_face_start + depth - (t_id + 1), height - VLEN,
                        depth, height, width, timestep,
                        INNER_FACE, LAST_ROW, tmp);
    }
}

static inline __attribute__((always_inline))
void gauss_reduce_row(const ftype *restrict w,
                      const ftype *restrict f_x,
                      const ftype *restrict f_y,
                      const ftype *restrict f_z,
                      const ftype *restrict u_x,
                      const ftype *restrict u_y,
                      const ftype *restrict u_z,
                      uint32_t width,
                      uint32_t row_stride,
                      ftype *restrict upper,
                      ftype *restrict dst_f_x,
                      ftype *restrict dst_f_y,
                      ftype *restrict dst_f_z)
{
    for (uint32_t i = 0; i < width; i += VLEN) {
        vftype ws = vload(w + i);
        vftype upper_prevs = vload(upper - row_stride + i);
        vftype f_x_prevs = vload(dst_f_x - row_stride + i);
        vftype f_y_prevs = vload(dst_f_y - row_stride + i);
        vftype f_z_prevs = vload(dst_f_z - row_stride + i);

        vftype fs_x = vsub(vload(f_x + i), vload(u_x + i));
        vftype fs_y = vsub(vload(f_y + i), vload(u_y + i));
        vftype fs_z = vsub(vload(f_z + i), vload(u_z + i));

        vftype norm_coefs = vfmadd(ws, upper_prevs, vadd(ONES, vadd(ws, ws)));
        vstore(upper + i, vdiv(vneg(ws), norm_coefs));
        vstore(dst_f_x + i, vdiv(vfmadd(ws, f_x_prevs, fs_x), norm_coefs));
        vstore(dst_f_y + i, vdiv(vfmadd(ws, f_y_prevs, fs_y), norm_coefs));
        vstore(dst_f_z + i, vdiv(vfmadd(ws, f_z_prevs, fs_z), norm_coefs));
    }
}

static inline __attribute__((always_inline))
void backward_sub_row(const ftype *restrict f_x,
                      const ftype *restrict f_y,
                      const ftype *restrict f_z,
                      const ftype *restrict upper,
                      uint32_t width,
                      ftype *restrict tmp_u_x,
                      ftype *restrict tmp_u_y,
                      ftype *restrict tmp_u_z,
                      ftype *restrict u_x,
                      ftype *restrict u_y,
                      ftype *restrict u_z)
{
    for (int k = 0; k < width; k += VLEN) {
        vftype fs_x = vload(f_x + k);
        vftype fs_y = vload(f_y + k);
        vftype fs_z = vload(f_z + k);
        vftype uppers = vload(upper + k);

        vftype u_x_prevs = vload(tmp_u_x + k);
        vftype u_y_prevs = vload(tmp_u_y + k);
        vftype u_z_prevs = vload(tmp_u_z + k);

        vftype us_x = vfmadd(vneg(uppers), u_x_prevs, fs_x);
        vftype us_y = vfmadd(vneg(uppers), u_y_prevs, fs_y);
        vftype us_z = vfmadd(vneg(uppers), u_z_prevs, fs_z);

        vstore(tmp_u_x + k, us_x);
        vstore(tmp_u_y + k, us_y);
        vstore(tmp_u_z + k, us_z);

        vstore(u_x + k, vadd(vload(u_x + k), us_x));
        vstore(u_y + k, vadd(vload(u_y + k), us_y));
        vstore(u_z + k, vadd(vload(u_z + k), us_z));
    }
}

static inline __attribute__((always_inline))
void apply_top_bc(uint32_t x,
                  uint32_t y,
                  uint32_t z,
                  uint32_t t,
                  ftype *restrict upper,
                  ftype *restrict f_x,
                  ftype *restrict f_y,
                  ftype *restrict f_z)
{
    vftype u0_x, u0_y, u0_z;
    _get_top_bc_u_delta(x, y, z, t, &u0_x, &u0_y, &u0_z);

    apply_start_bc(u0_y, u0_x, u0_z, upper, f_y, f_x, f_z);
}

static inline __attribute__((always_inline))
void apply_bottom_bc(const ftype *restrict w,
                     const ftype *restrict upper_prev,
                     const ftype *restrict f_x_prev,
                     const ftype *restrict f_z_prev,
                     const ftype *restrict f_x,
                     const ftype *restrict f_z,
                     uint32_t x,
                     uint32_t y,
                     uint32_t z,
                     uint32_t t,
                     ftype *restrict tmp_u_x,
                     ftype *restrict tmp_u_y,
                     ftype *restrict tmp_u_z,
                     ftype *restrict u_x,
                     ftype *restrict u_y,
                     ftype *restrict u_z)
{
    vftype un_x, un_y, un_z;
    _get_bottom_bc_u_delta(x, y, z, t, &un_x, &un_y, &un_z);

    vftype _un_x, _un_y, _un_z;
    apply_end_bc(w, upper_prev, f_x_prev, f_z_prev, f_x, f_z,
                 u_x, u_z, un_y, un_x, un_z, &_un_y, &_un_x, &_un_z);

    vstore(tmp_u_x, _un_x);
    vstore(tmp_u_y, _un_y);
    vstore(tmp_u_z, _un_z);

    vstore(u_x, vadd(vload(u_x), _un_x));
    vstore(u_y, vadd(vload(u_y), _un_y));
    vstore(u_z, vadd(vload(u_z), _un_z));
}

/* Solves the block diagonal system (I - wDyy)(u_n+1 - u_n) = f - u_n. */
static void solve_Dyy_blocks(const ftype *restrict w,
                             uint32_t depth,
                             uint32_t height,
                             uint32_t width,
                             uint32_t timestep,
                             ftype *restrict tmp,
                             const ftype *restrict f_x,
                             const ftype *restrict f_y,
                             const ftype *restrict f_z,
                             ftype *restrict u_x,
                             ftype *restrict u_y,
                             ftype *restrict u_z,
                             uint32_t t_id,
                             uint32_t num_threads,
                             int proc_rank)
{
    ZEROS = vbroadcast(0.0);
    ONES = vbroadcast(1.0);
    SIGN_MASK = vbroadcast(-0.0f);

    /* TODO: Can I prefetch something? */

    /* TODO: Try to perform cache blocking, reducing only a column
     * of tiles at a time, it should be better in terms of cache
     * locality of the intermediate coefficients but worse in terms
     * of TLB misses of the f and u. Additionally, it would reduce
     * the tmp storage required. Same holds for Dzz solver. */

    ftype *restrict tmp_f_x = tmp + width * height;
    ftype *restrict tmp_f_y = tmp + width * height * 2;
    ftype *restrict tmp_f_z = tmp + width * height * 3;

    ftype *restrict tmp_u_x = tmp + width * height * 4;
    ftype *restrict tmp_u_y = tmp + width * height * 4 + width;
    ftype *restrict tmp_u_z = tmp + width * height * 4 + width * 2;

    /* WARNING: Assumes depth is multiple of num_threads. */
    uint32_t block_depth = depth / num_threads;

    uint32_t global_face_start = proc_rank * depth;

    /* We solve for each face of the domain, one at a time. */
    for (int i = block_depth * t_id;
             i < block_depth * (t_id + 1); ++i) {
        /* Gauss reduce the first row. */
        uint64_t face_offset = height * width * i;

        /* Apply BCs on the first row of the domain. */
        for (uint32_t k = 0; k < width; k += VLEN) {
            apply_top_bc(k, 0, global_face_start + i, timestep, tmp + k,
                         tmp_f_x + k, tmp_f_y + k, tmp_f_z + k);
        }
        /* Gauss reduce the remaining face, one row at a time,
         * except the last one. */
        for (uint32_t j = 1; j < height - 1; ++j) {
            gauss_reduce_row(w + face_offset + width * j,
                             f_x + face_offset + width * j,
                             f_y + face_offset + width * j,
                             f_z + face_offset + width * j,
                             u_x + face_offset + width * j,
                             u_y + face_offset + width * j,
                             u_z + face_offset + width * j,
                             width,
                             width,
                             tmp + width * j,
                             tmp_f_x + width * j,
                             tmp_f_y + width * j,
                             tmp_f_z + width * j);
        }
        /* Apply BCs on the last row. */
        for (uint32_t k = 0; k < width; k += VLEN) {
            apply_bottom_bc(w + face_offset + width * (height - 1) + k,
                            tmp + width * (height - 2) + k,
                            tmp_f_x + width * (height - 2) + k,
                            tmp_f_z + width * (height - 2) + k,
                            f_x + face_offset + width * (height - 1) + k,
                            f_z + face_offset + width * (height - 1) + k,
                            k, height - 1, global_face_start + i, timestep,
                            tmp_u_x + k,
                            tmp_u_y + k,
                            tmp_u_z + k,
                            u_x + face_offset + width * (height - 1) + k,
                            u_y + face_offset + width * (height - 1) + k,
                            u_z + face_offset + width * (height - 1) + k);
        }

        /* Backward substitute the remaining face, one row at a time. */
        for (int j = 1; j < height; ++j) {
            uint64_t row_offset = face_offset + width * (height - j - 1);
            backward_sub_row(tmp_f_x + width * (height - j - 1),
                             tmp_f_y + width * (height - j - 1),
                             tmp_f_z + width * (height - j - 1),
                             tmp + row_offset - face_offset,
                             width,
                             tmp_u_x,
                             tmp_u_y,
                             tmp_u_z,
                             u_x + row_offset,
                             u_y + row_offset,
                             u_z + row_offset);
        }
    }
}

static inline __attribute__((always_inline))
void apply_front_bc(uint32_t x,
                    uint32_t y,
                    uint32_t z,
                    uint32_t t,
                    ftype *restrict upper,
                    ftype *restrict f_x,
                    ftype *restrict f_y,
                    ftype *restrict f_z)
{
    vftype u0_x, u0_y, u0_z;
    _get_front_bc_u_delta(x, y, z, t, &u0_x, &u0_y, &u0_z);

    apply_start_bc(u0_z, u0_x, u0_y, upper, f_z, f_x, f_y);
}

static inline __attribute__((always_inline))
void apply_back_bc(const ftype *restrict w,
                   const ftype *restrict upper_prev,
                   const ftype *restrict f_x_prev,
                   const ftype *restrict f_y_prev,
                   const ftype *restrict f_x,
                   const ftype *restrict f_y,
                   uint32_t x,
                   uint32_t y,
                   uint32_t z,
                   uint32_t t,
                   ftype *restrict tmp_u_x,
                   ftype *restrict tmp_u_y,
                   ftype *restrict tmp_u_z,
                   ftype *restrict u_x,
                   ftype *restrict u_y,
                   ftype *restrict u_z)
{
    vftype un_x, un_y, un_z;
    _get_back_bc_u_delta(x, y, z, t, &un_x, &un_y, &un_z);

    vftype _un_x, _un_y, _un_z;
    apply_end_bc(w, upper_prev, f_x_prev, f_y_prev, f_x, f_y,
                 u_x, u_y, un_z, un_x, un_y, &_un_z, &_un_x, &_un_y);

    vstore(tmp_u_x, _un_x);
    vstore(tmp_u_y, _un_y);
    vstore(tmp_u_z, _un_z);

    vstore(u_x, vadd(vload(u_x), _un_x));
    vstore(u_y, vadd(vload(u_y), _un_y));
    vstore(u_z, vadd(vload(u_z), _un_z));

}

#define NUM_COMM_CHUNKS 4

/* Solves the block diagonal system (I - wDzz)(u_n+1 - u_n) = f - u_n. */
static void solve_Dzz_blocks(const ftype *restrict w,
                             uint32_t depth,
                             uint32_t height,
                             uint32_t width,
                             uint32_t timestep,
                             ftype *restrict f_x,
                             ftype *restrict f_y,
                             ftype *restrict f_z,
                             ftype *restrict u_x,
                             ftype *restrict u_y,
                             ftype *restrict u_z,
                             uint32_t t_id,
                             uint32_t num_threads,
                             int proc_rank,
                             int num_procs,
                             MPI_Comm comm,
                             ArenaAllocator *restrict arena)
{
    /* TODO: Handle num_threads > 1. */
    assert(num_threads == 1);

    arena_enter(arena);

    /* WARNING: Assumes height is multiple of num_threads. */
    const uint32_t block_height = height / num_threads;

    const uint64_t block_size = depth * block_height * width;
    const uint64_t slice_size = block_height * width;

    const uint64_t send_size = 2 * slice_size;
    const uint64_t inner_size = block_size - send_size;

    /* Buffers for the inner local systems coefficients. */
    ftype *restrict low = arena_push_count(arena, ftype, inner_size);
    ftype *restrict upp = arena_push_count(arena, ftype, inner_size);
    ftype *restrict rhs_x = arena_push_count(arena, ftype, inner_size);
    ftype *restrict rhs_y = arena_push_count(arena, ftype, inner_size);
    ftype *restrict rhs_z = arena_push_count(arena, ftype, inner_size);

    /* Buffers for the reduced local systems coefficients. */
    ftype *restrict send_low = arena_push_count(arena, ftype, send_size);
    ftype *restrict send_upp = arena_push_count(arena, ftype, send_size);
    ftype *restrict send_rhs_x = arena_push_count(arena, ftype, send_size);
    ftype *restrict send_rhs_y = arena_push_count(arena, ftype, send_size);
    ftype *restrict send_rhs_z = arena_push_count(arena, ftype, send_size);

    const uint32_t block_row_start = block_height * t_id;
    const uint32_t block_row_end = block_height * (t_id + 1);

    const uint32_t global_face_start = proc_rank * depth;

    const int has_no_front_bc = proc_rank > 0;
    const int has_back_bc = proc_rank == num_procs - 1;

    /* TODO: Move sweeps code to individual functions. */

    /* Start the forward sweep on the first slice of the block,
     * applying front boundary conditions if needed. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            /* Index wrt locally owned subdomain. */
            uint64_t loc_idx = width * j + k;
            /* Index wrt send system block. */
            uint64_t send_idx = width * (j - block_row_start) + k;

            vftype w_loc = vload(w + loc_idx) * has_no_front_bc;
            vftype d_loc = 1 + 2 * w_loc;

            vstore(send_low + send_idx, -w_loc / d_loc);
            vstore(send_upp + send_idx, -w_loc / d_loc);

            vftype bc_u_x, bc_u_y, bc_u_z;
            _get_front_bc_u_delta(k, j, global_face_start, timestep,
                                  &bc_u_x, &bc_u_y, &bc_u_z);

            vftype rhs_x_loc = (vload(f_x + loc_idx) -
                                vload(u_x + loc_idx)) * has_no_front_bc +
                                              bc_u_x * !has_no_front_bc;
            vftype rhs_y_loc = (vload(f_y + loc_idx) -
                                vload(u_y + loc_idx)) * has_no_front_bc +
                                              bc_u_y * !has_no_front_bc;
            vftype rhs_z_loc = (vload(f_z + loc_idx) -
                                vload(u_z + loc_idx)) * has_no_front_bc +
                                              bc_u_z * !has_no_front_bc;

            /* I think that using an if on proc_rank outside the two
             * loops would basically perform the same. */
            vstore(send_rhs_x + send_idx, rhs_x_loc / d_loc);
            vstore(send_rhs_y + send_idx, rhs_y_loc / d_loc);
            vstore(send_rhs_z + send_idx, rhs_z_loc / d_loc);
        }
    }

    /* Continue forward sweep on the second slice of the block. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            /* Index wrt locally owned subdomain. */
            uint64_t loc_idx = height * width + width * j + k;
            /* Index wrt inner block, ignoring first/last slices. */
            uint64_t inner_idx = width * (j - block_row_start) + k;

            vftype w_loc = vload(w + loc_idx);
            vftype d_loc = 1 + 2 * w_loc;

            vstore(low + inner_idx, -w_loc / d_loc);
            vstore(upp + inner_idx, -w_loc / d_loc);

            vftype rhs_x_loc = vload(f_x + loc_idx) - vload(u_x + loc_idx);
            vftype rhs_y_loc = vload(f_y + loc_idx) - vload(u_y + loc_idx);
            vftype rhs_z_loc = vload(f_z + loc_idx) - vload(u_z + loc_idx);

            vstore(rhs_x + inner_idx, rhs_x_loc / d_loc);
            vstore(rhs_y + inner_idx, rhs_y_loc / d_loc);
            vstore(rhs_z + inner_idx, rhs_z_loc / d_loc);
        }
    }

    /* Complete forward sweep on the rest of the block, except
     * the last slice. */
    for (uint32_t i = 2; i < depth - 1; ++i) {
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t loc_idx = height * width * i + width * j + k;
                uint64_t inner_idx = block_height * width * (i - 1) +
                                     width * (j - block_row_start) + k;

                vftype w_loc = vload(w + loc_idx);
                vftype d_loc = 1 + 2 * w_loc;
                vftype low_loc = -w_loc;
                vftype upp_loc = -w_loc;

                /* r = 1 / (b[i] - a[i] * c[i - 1]) */
                vftype r_loc = (ftype) 1.0 /
                    (d_loc - low_loc * vload(upp + inner_idx - slice_size));

                /* a[i] *= r * -a[i - 1] */
                vstore(low + inner_idx,
                       r_loc * low_loc *
                       -vload(low + inner_idx - slice_size));
                /* c[i] *= r */
                vstore(upp + inner_idx, r_loc * upp_loc);

                vftype rhs_x_loc = vload(f_x + loc_idx) -
                                   vload(u_x + loc_idx);
                vftype rhs_y_loc = vload(f_y + loc_idx) -
                                   vload(u_y + loc_idx);
                vftype rhs_z_loc = vload(f_z + loc_idx) -
                                   vload(u_z + loc_idx);

                /* f[i] = r * (f[i] - f[i - 1] * a[i]) */
                vstore(rhs_x + inner_idx,
                       r_loc * (rhs_x_loc - low_loc *
                                vload(rhs_x + inner_idx - slice_size)));
                vstore(rhs_y + inner_idx,
                       r_loc * (rhs_y_loc - low_loc *
                                vload(rhs_y + inner_idx - slice_size)));
                vstore(rhs_z + inner_idx,
                       r_loc * (rhs_z_loc - low_loc *
                                vload(rhs_z + inner_idx - slice_size)));
            }
        }
    }

    /* Complete forward sweep by processing last slice, applying
     * back boundary conditions if needed. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            uint64_t loc_idx = height * width * (depth - 1) + width * j + k;
            uint64_t send_idx = block_height * width +
                                width * (j - block_row_start) + k;
            uint64_t inner_idx = block_height * width * (depth - 2) +
                                 width * (j - block_row_start) + k;

            vftype w_loc = vload(w + loc_idx);

            /* Only tangent component bc coefficients will be communicated.
             * Nevertheless the last row of the reduced system has no impact
             * on the values of the previous coefficients, since the backward
             * sweep starts at the third-last row. */

            /* Note that for tangent components, the equation with bc is:
             * -4/3 w_n u_n-1 + (1 + 4 w_n) u_n = 8/3 w_n u_ex + f_n */

            vftype d_loc = 1 + 2 * w_loc + 2 * w_loc * has_back_bc;
            vftype low_loc = -w_loc - (ftype) 1.0 / 3 * w_loc * has_back_bc;
            vftype upp_loc = -w_loc * !has_back_bc;

            vftype bc_u_x, bc_u_y, bc_u_z;
            _get_back_bc_u_delta(k, j, global_face_start + depth - 1,
                                 timestep, &bc_u_x, &bc_u_y, &bc_u_z);

            vftype rhs_x_loc = vload(f_x + loc_idx) -
                               vload(u_x + loc_idx) +
                               has_back_bc *
                               (ftype) 8.0 / 3 * w_loc * bc_u_x;

            vftype rhs_y_loc = vload(f_y + loc_idx) -
                               vload(u_y + loc_idx) +
                               has_back_bc *
                               (ftype) 8.0 / 3 * w_loc * bc_u_y;

            vftype rhs_z_loc = vload(f_z + loc_idx) -
                               vload(u_z + loc_idx);

            /* r = 1 / (b[i] - a[i] * c[i - 1]) */
            vftype r_loc = (ftype) 1.0 /
                (d_loc - low_loc * vload(upp + inner_idx - slice_size));

            /* a[i] *= r * -a[i - 1] */
            vstore(send_low + send_idx,
                   r_loc * low_loc * -vload(low + inner_idx - slice_size));
            /* c[i] *= r */
            vstore(send_upp + send_idx, r_loc * upp_loc);

            /* f[i] = r * (f[i] - f[i - 1] * a[i]) */
            vstore(send_rhs_x + send_idx,
                   r_loc * (rhs_x_loc - low_loc *
                            vload(rhs_x + inner_idx - slice_size)));

            vstore(send_rhs_y + send_idx,
                   r_loc * (rhs_y_loc - low_loc *
                            vload(rhs_y + inner_idx - slice_size)));

            vstore(send_rhs_z + send_idx,
                   !has_back_bc * r_loc * (rhs_z_loc -
                       low_loc * vload(rhs_z + inner_idx - slice_size)) +
                   has_back_bc * bc_u_z);
        }
    }

    /* TODO: I could already start communicating the local back
     * coefficients, since the backward sweep does not affect them.
     * I guess this would go against having a contiguous buffer
     * for the reduced system coefficients. */

    /* Perform backward sweep on the block, except for the first slice. */
    for (uint32_t i = 2; i < depth - 1; ++i) {
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t inner_idx = block_height * width * (depth - 2 - i) +
                                     width * (j - block_row_start) + k;

                vftype upp_loc = vload(upp + inner_idx);

                /* f[n - 1 - i] -= c[n - 1 - i] * f[n - i] */
                vstore(rhs_x + inner_idx,
                       vload(rhs_x + inner_idx) -
                       vload(rhs_x + inner_idx + slice_size) * upp_loc);
                vstore(rhs_y + inner_idx,
                       vload(rhs_y + inner_idx) -
                       vload(rhs_y + inner_idx + slice_size) * upp_loc);
                vstore(rhs_z + inner_idx,
                       vload(rhs_z + inner_idx) -
                       vload(rhs_z + inner_idx + slice_size) * upp_loc);

                /* a[n - 1 - i] -= c[n - 1 - i] * a[n - i] */
                vstore(low + inner_idx,
                       vload(low + inner_idx) -
                       vload(low + inner_idx + slice_size) * upp_loc);

                /* c[n - 1 - i] *= -c[n - i] */
                vstore(upp + inner_idx,
                       upp_loc * -vload(upp + inner_idx + slice_size));
            }
        }
    }

    /* Complete backward sweep on the first slice of the block. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            /* NOTE: The same idx can be used to access the coefficients
             * of the second slice (first slice of the inner block). */
            uint64_t send_idx = width * (j - block_row_start) + k;

            vftype upp_loc = vload(send_upp + send_idx);
            /* r = 1.0 / (1.0 - c[0] * a[1]) */
            vftype r_loc = (ftype) 1.0 / ((ftype) 1.0 -
                            upp_loc * vload(low + send_idx));

            /* f[0] = r * (f[0] - c[0] * f[1]) */
            vstore(send_rhs_x + send_idx,
                   r_loc * (vload(send_rhs_x + send_idx) -
                            vload(rhs_x + send_idx) * upp_loc));
            vstore(send_rhs_y + send_idx,
                   r_loc * (vload(send_rhs_y + send_idx) -
                            vload(rhs_y + send_idx) * upp_loc));
            vstore(send_rhs_z + send_idx,
                   r_loc * (vload(send_rhs_z + send_idx) -
                            vload(rhs_z + send_idx) * upp_loc));

            /* a[0] *= r */
            vstore(send_low + send_idx, r_loc * vload(send_low + send_idx));

            /* c[0] *= r * -c[1] */
            vstore(send_upp + send_idx,
                   r_loc * upp_loc * -vload(upp + send_idx));
        }
    }

    const uint64_t redu_size = num_procs * send_size;
    /* Buffers for the reduced systems coefficients. */
    ftype *restrict redu_low = arena_push_count(arena, ftype, redu_size);
    ftype *restrict redu_upp = arena_push_count(arena, ftype, redu_size);
    ftype *restrict redu_rhs_x = arena_push_count(arena, ftype, redu_size);
    ftype *restrict redu_rhs_y = arena_push_count(arena, ftype, redu_size);
    ftype *restrict redu_rhs_z = arena_push_count(arena, ftype, redu_size);

    /* Assemble the reduced system by communicating with other procs. */

    /* TODO: Partition block into sub-blocks, overlap sub-block sweeps
     * with communication using MPI_Iallgather. */

    /* WARNING: This works as long as num_threads = 1 */
    MPI_Allgather(send_low, send_size, MPI_FTYPE,
                  redu_low, send_size, MPI_FTYPE, comm);
    MPI_Allgather(send_upp, send_size, MPI_FTYPE,
                  redu_upp, send_size, MPI_FTYPE, comm);
    MPI_Allgather(send_rhs_x, send_size, MPI_FTYPE,
                  redu_rhs_x, send_size, MPI_FTYPE, comm);
    MPI_Allgather(send_rhs_y, send_size, MPI_FTYPE,
                  redu_rhs_y, send_size, MPI_FTYPE, comm);
    MPI_Allgather(send_rhs_z, send_size, MPI_FTYPE,
                  redu_rhs_z, send_size, MPI_FTYPE, comm);

    /* Solve the reduced system, note that reduced system
     * has 1s on the diagonal by construction. */

    /* Forward sweep on the reduced system, except last slice. */
    for (uint32_t i = 1; i < num_procs * 2 - 1; ++i) {
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t redu_idx = block_height * width * i +
                                    width * (j - block_row_start) + k;

                vftype low_loc = vload(redu_low + redu_idx);
                /* r = 1.0 / (b[i] - c[i - 1] * a[i]) */
                vftype r_loc = (ftype) 1.0 / ((ftype) 1.0 - low_loc *
                               vload(redu_upp + redu_idx - slice_size));

                /* f[i] = r * (f[i] - f[i - 1] * a[i]) */
                vstore(redu_rhs_x + redu_idx,
                       r_loc * (vload(redu_rhs_x + redu_idx) - low_loc *
                                vload(redu_rhs_x + redu_idx - slice_size)));
                vstore(redu_rhs_y + redu_idx,
                       r_loc * (vload(redu_rhs_y + redu_idx) - low_loc *
                                vload(redu_rhs_y + redu_idx - slice_size)));
                vstore(redu_rhs_z + redu_idx,
                       r_loc * (vload(redu_rhs_z + redu_idx) - low_loc *
                                vload(redu_rhs_z + redu_idx - slice_size)));

                /* c[i] *= r */
                vstore(redu_upp + redu_idx,
                       r_loc * vload(redu_upp + redu_idx));
            }
        }
    }

    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            uint64_t redu_idx = block_height * width *
                                (num_procs * 2 - 1) +
                                width * (j - block_row_start) + k;

            vftype low_loc = vload(redu_low + redu_idx);
            /* r = 1.0 / (b[i] - c[i - 1] * a[i]) */
            vftype r_loc = (ftype) 1.0 / ((ftype) 1.0 - low_loc *
                           vload(redu_upp + redu_idx - slice_size));

            /* f[i] = r * (f[i] - f[i - 1] * a[i]) */
            vstore(redu_rhs_x + redu_idx,
                   r_loc * (vload(redu_rhs_x + redu_idx) - low_loc *
                            vload(redu_rhs_x + redu_idx - slice_size)));
            vstore(redu_rhs_y + redu_idx,
                   r_loc * (vload(redu_rhs_y + redu_idx) - low_loc *
                            vload(redu_rhs_y + redu_idx - slice_size)));
            /* Last row of the reduced system HAS back bc! */
            //vstore(redu_rhs_z + redu_idx, vload(redu_rhs_z + redu_idx));
        }
    }

    /* Backward sweep on the reduced system. */
    for (uint32_t i = 1; i < num_procs * 2; ++i) {
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t redu_idx = block_height * width *
                                    (num_procs * 2 - 1 - i) +
                                    width * (j - block_row_start) + k;

                vftype upp_loc = vload(redu_upp + redu_idx);
                /* f[n - 1 - i] -= c[n - 1 - i] * f[n - i] */
                vstore(redu_rhs_x + redu_idx,
                       vload(redu_rhs_x + redu_idx) - upp_loc *
                       vload(redu_rhs_x + redu_idx + slice_size));
                vstore(redu_rhs_y + redu_idx,
                       vload(redu_rhs_y + redu_idx) - upp_loc *
                       vload(redu_rhs_y + redu_idx + slice_size));
                vstore(redu_rhs_z + redu_idx,
                       vload(redu_rhs_z + redu_idx) - upp_loc *
                       vload(redu_rhs_z + redu_idx + slice_size));
            }
        }
    }

    /* Substitute reduced system solution in the local system. */

    /* NOTE: u is padded, so this is safe */
    ftype *restrict u_x_front_halo = u_x - height * width;
    ftype *restrict u_y_front_halo = u_y - height * width;
    ftype *restrict u_z_front_halo = u_z - height * width;

    if (proc_rank > 0) {
        /* Update front halo. */
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t halo_idx = width * j + k;
                uint64_t redu_idx = block_height * width *
                                    (2 * proc_rank - 1) +
                                    width * (j - block_row_start) + k;

                vstore(u_x_front_halo + halo_idx,
                       vload(u_x_front_halo + halo_idx) +
                       vload(redu_rhs_x + redu_idx));
                vstore(u_y_front_halo + halo_idx,
                       vload(u_y_front_halo + halo_idx) +
                       vload(redu_rhs_y + redu_idx));
                vstore(u_z_front_halo + halo_idx,
                       vload(u_z_front_halo + halo_idx) +
                       vload(redu_rhs_z + redu_idx));
            }
        }
    }

    /* Update first slice. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            uint64_t loc_idx = width * j + k;
            uint64_t redu_idx = block_height * width * 2 * proc_rank +
                                width * (j - block_row_start) + k;

            vstore(u_x + loc_idx,
                   vload(u_x + loc_idx) + vload(redu_rhs_x + redu_idx));
            vstore(u_y + loc_idx,
                   vload(u_y + loc_idx) + vload(redu_rhs_y + redu_idx));
            vstore(u_z + loc_idx,
                   vload(u_z + loc_idx) + vload(redu_rhs_z + redu_idx));
        }
    }

    /* Update inner block. */
    /* TODO: Consider sweeping backwards, last slices could
     * still be in cache. */

    for (uint32_t i = 1; i < depth - 1; ++i) {
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t loc_idx = height * width * i + width * j + k;
                uint64_t inner_idx = block_height * width * (i - 1) + 
                                     width * (j - block_row_start) + k;

                uint64_t redu_front_idx =
                    block_height * width * 2 * proc_rank +
                    width * (j - block_row_start) + k;
                uint64_t redu_back_idx =
                    block_height * width * (2 * proc_rank + 1) +
                    width * (j - block_row_start) + k;

                vftype low_loc = vload(low + inner_idx);
                vftype upp_loc = vload(upp + inner_idx);

                /* f[i] -= (u0 * a[i] + um * c[i]) */
                vstore(u_x + loc_idx,
                       vload(u_x + loc_idx) + vload(rhs_x + inner_idx) -
                       vload(redu_rhs_x + redu_front_idx) * low_loc -
                       vload(redu_rhs_x + redu_back_idx) * upp_loc);
                vstore(u_y + loc_idx,
                       vload(u_y + loc_idx) + vload(rhs_y + inner_idx) -
                       vload(redu_rhs_y + redu_front_idx) * low_loc -
                       vload(redu_rhs_y + redu_back_idx) * upp_loc);
                vstore(u_z + loc_idx,
                       vload(u_z + loc_idx) + vload(rhs_z + inner_idx) -
                       vload(redu_rhs_z + redu_front_idx) * low_loc -
                       vload(redu_rhs_z + redu_back_idx) * upp_loc);
            }
        }
    }

    /* Update last slice. */
    for (uint32_t j = block_row_start; j < block_row_end; ++j) {
        for (uint32_t k = 0; k < width; k += VLEN) {

            uint64_t loc_idx = height * width * (depth - 1) + width * j + k;
            uint64_t redu_idx = block_height * width * (2 * proc_rank + 1) +
                                width * (j - block_row_start) + k;

            vstore(u_x + loc_idx,
                   vload(u_x + loc_idx) + vload(redu_rhs_x + redu_idx));
            vstore(u_y + loc_idx,
                   vload(u_y + loc_idx) + vload(redu_rhs_y + redu_idx));
            vstore(u_z + loc_idx,
                   vload(u_z + loc_idx) + vload(redu_rhs_z + redu_idx));
        }
    }

    /* NOTE: u is padded, so this is safe */
    ftype *restrict u_x_back_halo = u_x + depth * height * width;
    ftype *restrict u_y_back_halo = u_y + depth * height * width;
    ftype *restrict u_z_back_halo = u_z + depth * height * width;

    if (proc_rank < num_procs - 1) {
        /* Update back halo, note that u is padded, so this is safe. */
        for (uint32_t j = block_row_start; j < block_row_end; ++j) {
            for (uint32_t k = 0; k < width; k += VLEN) {

                uint64_t halo_idx = width * j + k;
                uint64_t redu_idx = block_height * width *
                                    (2 * proc_rank + 2) +
                                    width * (j - block_row_start) + k;

                vstore(u_x_back_halo + halo_idx,
                       vload(u_x_back_halo + halo_idx) +
                       vload(redu_rhs_x + redu_idx));
                vstore(u_y_back_halo + halo_idx,
                       vload(u_y_back_halo + halo_idx) +
                       vload(redu_rhs_y + redu_idx));
                vstore(u_z_back_halo + halo_idx,
                       vload(u_z_back_halo + halo_idx) +
                       vload(redu_rhs_z + redu_idx));
            }
        }
    }

    arena_exit(arena);
}

void momentum_init(field_size size, field3 field)
{
    field3_fill(size, 0, field);

    /* WARNING: Dummy initialization for lid-driven cavity. */
    for (uint32_t i = 0; i < size.depth; ++i) {
        for (uint32_t k = 0; k < size.width; ++k) {
            uint64_t idx = field_idx(size, k, 0, i);
            field.z[idx] = 1.0;
        }
    }
}

void momentum_solve(const_field porosity,
                    const_field gamma,
                    const_field pressure,
                    const_field pressure_delta,
                    field_size size,
                    field3 velocity_Dxx,
                    field3 velocity_Dyy,
                    field3 velocity_Dzz,
                    uint32_t timestep,
                    Thread *thread,
                    DDecomp *ddecomp)
{
    ArenaAllocator *arena = thread_get_arena(thread);
    arena_enter(arena);

    uint32_t t_id = thread->t_id;
    uint32_t num_threads = thread_get_array_size(thread);

    /* TODO: Check actual size needed when num_threads > 1. */
    field_size tmp_size = { size.width,
                            size.height / num_threads,
                            size.depth * 4 + 3 };

    field tmp = field_alloc(tmp_size, arena);

    int proc_rank = get_proc_rank(ddecomp->comm_z);
    int num_procs = get_num_procs(ddecomp->comm_z);

    TIMER_CREATE(solve_momentum_Dxx_blocks_fused_rhs);
    TIMER_CREATE(solve_momentum_Dyy_blocks);
    TIMER_CREATE(solve_momentum_Dzz_blocks);

    TIMER_RESTART(solve_momentum_Dxx_blocks_fused_rhs);

    solve_Dxx_blocks_fused_rhs(
        porosity, gamma, pressure, pressure_delta,
        velocity_Dxx.x, velocity_Dxx.y, velocity_Dxx.z,
        velocity_Dyy.x, velocity_Dyy.y, velocity_Dyy.z,
        velocity_Dzz.x, velocity_Dzz.y, velocity_Dzz.z,
        size.depth, size.height, size.width, timestep, tmp,
        t_id, num_threads, proc_rank, num_procs);

    thread_wait_on_barrier(thread);
    TIMER_ELAPSED(solve_momentum_Dxx_blocks_fused_rhs,
                  t_id == 0 && proc_rank == 0);

    TIMER_RESTART(solve_momentum_Dyy_blocks);

    solve_Dyy_blocks(
        gamma, size.depth, size.height, size.width, timestep,
        tmp, velocity_Dxx.x, velocity_Dxx.y, velocity_Dxx.z,
        velocity_Dyy.x, velocity_Dyy.y, velocity_Dyy.z,
        t_id, num_threads, proc_rank);

    thread_wait_on_barrier(thread);
    TIMER_ELAPSED(solve_momentum_Dyy_blocks, t_id == 0 && proc_rank == 0);

    arena_exit(arena);

    TIMER_RESTART(solve_momentum_Dzz_blocks);

    solve_Dzz_blocks(
        gamma, size.depth, size.height, size.width, timestep,
        velocity_Dyy.x, velocity_Dyy.y, velocity_Dyy.z,
        velocity_Dzz.x, velocity_Dzz.y, velocity_Dzz.z,
        t_id, num_threads, proc_rank, num_procs, ddecomp->comm_z, arena);

    thread_wait_on_barrier(thread);
    TIMER_ELAPSED(solve_momentum_Dzz_blocks, t_id == 0 && proc_rank == 0);

}
