#include <math.h>

#include "test.h"
#include "ftype.h"
#include "alloc.h"
#include "field.h"
#include "consts.h"
#include "output.h"
#include "boundary.h"
#include "momentum.h"
#include "pressure.h"
#include "thread-array.h"
#include "convergence-test.h"

DEFINE_DX(1.0)
DEFINE_DT(1.0)
DEFINE_NU(1.0)

static ftype get_man_u_x(ftype x, ftype y, ftype z, ftype t)
{
    return sin(x) * cos(y + t) * sin(z);
}

static ftype get_man_u_y(ftype x, ftype y, ftype z, ftype t)
{
    return cos(x) * sin(y + t) * sin(z);
}

static ftype get_man_u_z(ftype x, ftype y, ftype z, ftype t)
{
    return 2 * cos(x) * cos(y + t) * cos(z);
}

static ftype get_porosity(ftype x, ftype y, ftype z, ftype t)
{
    // Values of k > 1e-4 are safe
    return 2 + sin(M_PI * x) * cos(M_PI * y + t) * sin(M_PI * z);
}

DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_LEFT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_RIGHT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_TOP)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_BOTTOM)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_FRONT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_BACK)

static ftype get_forcing_x(ftype x, ftype y, ftype z, ftype t)
{
    // NOTE: Currently k is evaluated at the center of the cell
    ftype k = get_porosity(x - _DX / 2, y, z, t);

    return sin(x) * (-sin(y + t) * sin(z) +
                      cos(y + t) * sin(z) * _NU * (3 + 1.0 / k) +
                      cos(y + t) * cos(z) * -3 * _NU);
}

static ftype get_forcing_y(ftype x, ftype y, ftype z, ftype t)
{
    ftype k = get_porosity(x, y - _DX / 2, z, t);

    return cos(x) * (cos(y + t) * sin(z) +
                     sin(y + t) * sin(z) * _NU * (3 + 1.0 / k) +
                     sin(y + t) * cos(z) * -3 * _NU);
}

static ftype get_forcing_z(ftype x, ftype y, ftype z, ftype t)
{
    ftype k = get_porosity(x, y, z - _DX / 2, t);

    return cos(x) * (sin(y + t) * cos(z) * -2 +
                     cos(y + t) * cos(z) * 2 * _NU * (3 + 1.0 / k) +
                     cos(y + t) * sin(z) * -3 * _NU);
}

DEFINE_FORCING(get_forcing_x, get_forcing_y, get_forcing_z)

static void compute_manufactured_solution(field_size size,
                                          uint32_t timestep,
                                          field3 dst)
{
    ftype time = timestep * _DT;
    for (uint32_t i = 0; i < size.depth; ++i) {
        for (uint32_t j = 0; j < size.height; ++j) {
            for (uint32_t k = 0; k < size.width; ++k) {
                uint64_t idx = size.height * size.width * i +
                               size.width * j + k;

                dst.x[idx] =
                    get_man_u_x(_DX * k + _DX / 2, _DX * j, _DX * i, time);
                dst.y[idx] =
                    get_man_u_y(_DX * k, _DX * j + _DX / 2, _DX * i, time);
                dst.z[idx] =
                    get_man_u_z(_DX * k, _DX * j, _DX * i + _DX / 2, time);
            }
        }
    }
}

static void compute_manufactured_pressure(field_size size,
                                          ftype timestep,
                                          field dst)
{
    ftype t = timestep * _DT - _DT / 2;

    for (uint32_t i = 0; i < size.depth; ++i) {
        for (uint32_t j = 0; j < size.height; ++j) {
            for (uint32_t k = 0; k < size.width; ++k) {
                uint64_t idx = size.height * size.width * i +
                               size.width * j + k;

                dst[idx] = 3 * _NU * cos(k * _DX)
                                   * cos(j * _DX + t)
                                   * cos(i * _DX);
            }
        }
    }
}

struct SolverData {
    int num_timesteps;
    field_size size;
    field porosity;
    field gamma;
    field pressure;
    field pressure_pred;
    field3 eta;
    field3 zeta;
    field3 vel;
    OutputVTK *output;
};

static void *solve_brinkman(void *thread)
{
    /*ArenaAllocator *arena = thread_get_arena(thread);*/
    struct SolverData *solver_data = thread_get_shared_data(thread);

    int num_timesteps = solver_data->num_timesteps;
    field_size size = solver_data->size;
    field porosity = solver_data->porosity;
    field gamma = solver_data->gamma;
    field pressure = solver_data->pressure;
    field pressure_pred = solver_data->pressure_pred;
    field3 eta = solver_data->eta;
    field3 zeta = solver_data->zeta;
    field3 vel = solver_data->vel;
    OutputVTK *output = solver_data->output;

    /*
    char output_file_name[32];
    sprintf(output_file_name, "output/solution-%.4f-%d.vtk", _DT, 0);
    output_vtk_write(output, output_file_name, thread);
    */

    for (uint32_t t = 1; t < num_timesteps + 1; ++t) {

        for (uint32_t z = 0; z < size.depth; ++z) {
            for (uint32_t y = 0; y < size.height; ++y) {
                for (uint32_t x = 0; x < size.width; ++x) {
                    // Updating time dependent porosity
                    uint64_t idx = field_idx(size, x, y, z);
                    ftype k = get_porosity(x * _DX, y * _DX, z * _DX, t * _DT);

                    porosity[idx] = k;
                    gamma[idx] = (_DT * _NU) / (2 + _DT * _NU / k) / (_DX * _DX);
                }
            }
        }

        momentum_solve(porosity, gamma, pressure_pred,
                       size, eta, zeta, vel, t, thread);

        /*
        arena_enter(arena);

        field3 vel_old = field3_alloc(size, arena);

        for (uint32_t i = 0; i < size.depth; ++i) {
            for (uint32_t j = 0; j < size.height; ++j) {
                for (uint32_t k = 0; k < size.width; ++k) {
                    uint64_t idx = size.height * size.width * i +
                                   size.width * j + k;

                    vel_old.x[idx] = vel.x[idx];
                    vel_old.y[idx] = vel.y[idx];
                    vel_old.z[idx] = vel.z[idx];

                }
            }
        }

        arena_exit(arena);
        */

        pressure_solve(to_const_field3(vel), size, pressure,
                       pressure_pred, t, thread);
 
        thread_wait_on_barrier(thread);

        /*
        pressure_correct_rot(to_const_field3(vel),
                             to_const_field3(vel_old),
                             size, pressure, 0.0, t);
        */

        /*
        char output_file_name[32];
        sprintf(output_file_name, "output/solution-%.4f-%d.vtk", _DT, t);
        output_vtk_write(output, output_file_name, thread);
        */
    }

    return 0;
}

DEF_TEST(test_convergence_space,
         ArenaAllocator *arena,
         int num_samples,
         int num_threads)
{
    arena_enter(arena);

    ThreadArray *t_array = thread_array_create(num_threads, arena);

    double *v_errors = arena_push_count(arena, double, num_samples);
    double *p_errors = arena_push_count(arena, double, num_samples);
    double *dxs = arena_push_count(arena, double, num_samples);

    double T = 5e-6;
    double dt = 1e-6;

    SET_DT(dt);

    for (int i = 0; i < num_samples; ++i) {
        arena_enter(arena);

        field_size size = { 16 << i, 16 << i, 16 << i };
        SET_DX(1.0 / (size.width - 0.5));

        dxs[i] = _DX;

        field porosity = field_alloc(size, arena);
        field gamma = field_alloc(size, arena);

        for (uint32_t z = 0; z < size.depth; ++z) {
            for (uint32_t y = 0; y < size.height; ++y) {
                for (uint32_t x = 0; x < size.width; ++x) {

                    uint64_t idx = field_idx(size, x, y, z);
                    ftype k = get_porosity(x * _DX, y * _DX, z * _DX, 0);

                    porosity[idx] = k;
                    gamma[idx] = (_DT * _NU) / (2 + _DT * _NU / k) / (_DX * _DX);
                }
            }
        }

        field3 eta = field3_alloc_pad(size, arena);
        field3 zeta = field3_alloc_pad(size, arena);
        field3 vel = field3_alloc_pad(size, arena);

        field3_fill(size, 0, eta);
        field3_fill(size, 0, zeta);
        field3_fill(size, 0, vel);

        field pressure = field_alloc(size, arena);
        field pressure_pred = field_alloc(size, arena);

        compute_manufactured_solution(size, 0, eta);
        compute_manufactured_solution(size, 0, zeta);
        compute_manufactured_solution(size, 0, vel);
        compute_manufactured_pressure(size, 0.5, pressure);
        compute_manufactured_pressure(size, 0.5, pressure_pred);

        field manufactured_pressure = field_alloc(size, arena);
        compute_manufactured_pressure(size, 0.5, manufactured_pressure);

        OutputVTK *output = output_vtk_create(size, _DX, arena);

        output_vtk_attach_field(output, pressure, "pressure", arena);
        //output_vtk_attach_field(output, manufactured_pressure,
        //                        "man_pressure", arena);
        output_vtk_attach_field3(output, to_const_field3(vel),
                                 "velocity", arena);
        output_vtk_attach_field(output, porosity, "porosity", arena);

        uint32_t num_timesteps = round(T / _DT);

        struct SolverData data = {
            num_timesteps,
            size,
            porosity,
            gamma,
            pressure,
            pressure_pred,
            eta,
            zeta,
            vel,
            output
        };

        thread_array_set_shared_data(t_array, &data);
        thread_array_run(t_array, solve_brinkman, arena);

        field3 manufactured_v = field3_alloc(size, arena);
        compute_manufactured_solution(size, num_timesteps, manufactured_v);
        v_errors[i] = field3_l2_norm_diff(size,
                                          _DX,
                                          to_const_field3(vel),
                                          to_const_field3(manufactured_v));

        ftype mean_pressure = 0;
        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            mean_pressure += pressure[i];
        }
        mean_pressure /= field_num_points(size);

        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            pressure[i] -= mean_pressure;
        }

        compute_manufactured_pressure(size, num_timesteps,
                                      manufactured_pressure);
        mean_pressure = 0;
        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            mean_pressure += manufactured_pressure[i];
        }
        mean_pressure /= field_num_points(size);

        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            manufactured_pressure[i] -= mean_pressure;
        }

        p_errors[i] = field_l2_norm_diff(size,
                                         _DX,
                                         pressure,
                                         manufactured_pressure);

        arena_exit(arena);
    }

    double *v_orders = arena_push_count(arena, double, num_samples);
    double *p_orders = arena_push_count(arena, double, num_samples);
    estimate_convergence_order(v_errors, dxs, num_samples, v_orders);
    estimate_convergence_order(p_errors, dxs, num_samples, p_orders);

    printf("%.8f %e   --  %e   --\n", dxs[0], v_errors[0], p_errors[0]);
    for (int i = 1; i < num_samples; ++i) {
        printf("%.8f %e %5.2f %e %5.2f\n",
               dxs[i], v_errors[i], v_orders[i], p_errors[i], p_orders[i]);
    }

    thread_array_destroy(t_array);

    arena_exit(arena);
}

DEF_TEST(test_convergence_time,
         ArenaAllocator *arena,
         int num_samples,
         int num_threads)
{
    arena_enter(arena);

    ThreadArray *t_array = thread_array_create(num_threads, arena);

    double *v_errors = arena_push_count(arena, double, num_samples);
    double *p_errors = arena_push_count(arena, double, num_samples);
    double *dts = arena_push_count(arena, double, num_samples);

    double T = 1.0;
    double dt = 0.1;

    field_size size = { 64, 64, 64 };
    SET_DX(M_PI / (size.width - 0.5));

    for (int i = 0; i < num_samples; ++i) {
        arena_enter(arena);

        SET_DT(dt);

        field porosity = field_alloc(size, arena);
        field gamma = field_alloc(size, arena);

        for (uint32_t z = 0; z < size.depth; ++z) {
            for (uint32_t y = 0; y < size.height; ++y) {
                for (uint32_t x = 0; x < size.width; ++x) {

                    uint64_t idx = field_idx(size, x, y, z);
                    ftype k = get_porosity(x * _DX, y * _DX, z * _DX, 0);

                    porosity[idx] = k;
                    gamma[idx] = (_DT * _NU) / (2 + _DT * _NU / k) / (_DX * _DX);
                }
            }
        }

        field3 eta = field3_alloc_pad(size, arena);
        field3 zeta = field3_alloc_pad(size, arena);
        field3 vel = field3_alloc_pad(size, arena);

        field3_fill(size, 0, eta);
        field3_fill(size, 0, zeta);
        field3_fill(size, 0, vel);

        field pressure = field_alloc(size, arena);
        field pressure_pred = field_alloc(size, arena);

        compute_manufactured_solution(size, 0, eta);
        compute_manufactured_solution(size, 0, zeta);
        compute_manufactured_solution(size, 0, vel);
        compute_manufactured_pressure(size, 0.5, pressure);
        compute_manufactured_pressure(size, 0.5, pressure_pred);

        field manufactured_pressure = field_alloc(size, arena);
        compute_manufactured_pressure(size, 0.5, manufactured_pressure);

        uint32_t num_timesteps = round(T / _DT);

        struct SolverData data = {
            num_timesteps,
            size,
            porosity,
            gamma,
            pressure,
            pressure_pred,
            eta,
            zeta,
            vel
        };

        thread_array_set_shared_data(t_array, &data);
        thread_array_run(t_array, solve_brinkman, arena);

        field3 manufactured_v = field3_alloc(size, arena);
        compute_manufactured_solution(size, num_timesteps, manufactured_v);
        v_errors[i] = field3_l2_norm_diff(size,
                                          _DX,
                                          to_const_field3(vel),
                                          to_const_field3(manufactured_v));

        ftype mean_pressure = 0;
        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            mean_pressure += pressure[i];
        }
        mean_pressure /= field_num_points(size);

        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            pressure[i] -= mean_pressure;
        }

        compute_manufactured_pressure(size, num_timesteps,
                                      manufactured_pressure);
        mean_pressure = 0;
        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            mean_pressure += manufactured_pressure[i];
        }
        mean_pressure /= field_num_points(size);

        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            manufactured_pressure[i] -= mean_pressure;
        }

        p_errors[i] = field_l2_norm_diff(size,
                                         _DX,
                                         pressure,
                                         manufactured_pressure);

        dts[i] = _DT;
        dt /= 2;

        arena_exit(arena);
    }

    double *v_orders = arena_push_count(arena, double, num_samples);
    double *p_orders = arena_push_count(arena, double, num_samples);
    estimate_convergence_order(v_errors, dts, num_samples, v_orders);
    estimate_convergence_order(p_errors, dts, num_samples, p_orders);

    printf("%.8f %e   --  %e   --\n", dts[0], v_errors[0], p_errors[0]);
    for (int i = 1; i < num_samples; ++i) {
        printf("%.8f %e %5.2f %e %5.2f\n",
               dts[i], v_errors[i], v_orders[i], p_errors[i], p_orders[i]);
    }

    thread_array_destroy(t_array);

    arena_exit(arena);
}

int main(void)
{
    ArenaAllocator arena;
    arena_init(&arena, 1ul << 33);

    RUN_TEST(test_convergence_space, &arena, 4, 4);
    RUN_TEST(test_convergence_time, &arena, 4, 4);

    arena_destroy(&arena);
}
