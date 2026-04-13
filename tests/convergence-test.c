#include <math.h>
#include <mpi.h>

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
#include "ddecomp.h"
#include "convergence-test.h"

DEFINE_DX(1.0)
DEFINE_DT(1.0)
DEFINE_NU(1.0)

static const ftype _K = 1.0;

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

DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_LEFT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_RIGHT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_TOP)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_BOTTOM)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_FRONT)
DEFINE_FUNCTION_BC_U(get_man_u_x, get_man_u_y, get_man_u_z, BC_BACK)

static ftype get_forcing_x(ftype x, ftype y, ftype z, ftype t)
{
    /* WARNING: _K here is constant! */

    return sin(x) * (-sin(y + t) * sin(z) +
                      cos(y + t) * sin(z) * _NU * (3 + 1.0 / _K) +
                     cos(y + t) * cos(z) * -3 * _NU);
}

static ftype get_forcing_y(ftype x, ftype y, ftype z, ftype t)
{
    return cos(x) * (cos(y + t) * sin(z) +
                     sin(y + t) * sin(z) * _NU * (3 + 1.0 / _K) +
                     sin(y + t) * cos(z) * -3 * _NU);
}

static ftype get_forcing_z(ftype x, ftype y, ftype z, ftype t)
{
    return cos(x) * (sin(y + t) * cos(z) * -2 +
                     cos(y + t) * cos(z) * 2 * _NU * (3 + 1.0 / _K) +
                     cos(y + t) * sin(z) * -3 * _NU);
}

DEFINE_FORCING(get_forcing_x, get_forcing_y, get_forcing_z)

static void compute_manufactured_solution(field_size size,
                                          uint32_t timestep,
                                          field3 dst)
{
    int proc_rank = get_proc_rank(MPI_COMM_WORLD);
    uint32_t global_i_start = size.depth * proc_rank;

    ftype time = timestep * _DT;
    for (uint32_t i = 0; i < size.depth; ++i) {
        for (uint32_t j = 0; j < size.height; ++j) {
            for (uint32_t k = 0; k < size.width; ++k) {
                uint64_t idx = size.height * size.width * i +
                               size.width * j + k;

                dst.x[idx] =
                    get_man_u_x(_DX * k + _DX / 2, _DX * j, _DX * (global_i_start + i), time);
                dst.y[idx] =
                    get_man_u_y(_DX * k, _DX * j + _DX / 2, _DX * (global_i_start + i), time);
                dst.z[idx] =
                    get_man_u_z(_DX * k, _DX * j, _DX * (global_i_start + i) + _DX / 2, time);
            }
        }
    }

    ftype *front_halo_x = dst.x - size.height * size.width;
    ftype *front_halo_y = dst.y - size.height * size.width;
    ftype *front_halo_z = dst.z - size.height * size.width;

    ftype *back_halo_x = dst.x + size.depth * size.height * size.width;
    ftype *back_halo_y = dst.y + size.depth * size.height * size.width;
    ftype *back_halo_z = dst.z + size.depth * size.height * size.width;

    for (uint32_t j = 0; j < size.height; ++j) {
        for (uint32_t k = 0; k < size.width; ++k) {
            uint64_t halo_idx = size.width * j + k;

            front_halo_x[halo_idx] =
                get_man_u_x(_DX * k + _DX / 2, _DX * j, _DX * (global_i_start - 1), time);
            front_halo_y[halo_idx] =
                get_man_u_y(_DX * k, _DX * j + _DX / 2, _DX * (global_i_start - 1), time);
            front_halo_z[halo_idx] =
                get_man_u_z(_DX * k, _DX * j, _DX * (global_i_start - 1) + _DX / 2, time);

            back_halo_x[halo_idx] =
                get_man_u_x(_DX * k + _DX / 2, _DX * j, _DX * (global_i_start + size.depth), time);
            back_halo_y[halo_idx] =
                get_man_u_y(_DX * k, _DX * j + _DX / 2, _DX * (global_i_start + size.depth), time);
            back_halo_z[halo_idx] =
                get_man_u_z(_DX * k, _DX * j, _DX * (global_i_start + size.depth) + _DX / 2, time);
        }
    }


}

static void compute_manufactured_pressure(field_size size,
                                          ftype timestep,
                                          field dst)
{
    ftype t = timestep * _DT - _DT / 2;

    int proc_rank = get_proc_rank(MPI_COMM_WORLD);
    uint32_t global_i_start = size.depth * proc_rank;

    for (uint32_t i = 0; i < size.depth; ++i) {
        for (uint32_t j = 0; j < size.height; ++j) {
            for (uint32_t k = 0; k < size.width; ++k) {
                uint64_t idx = size.height * size.width * i +
                               size.width * j + k;

                dst[idx] = 3 * _NU * cos(k * _DX)
                                   * cos(j * _DX + t)
                                   * cos((global_i_start + i) * _DX);
            }
        }
    } 

    ftype *front_halo = dst - size.height * size.width;
    ftype *back_halo = dst + size.depth * size.height * size.width;

    for (uint32_t j = 0; j < size.height; ++j) {
        for (uint32_t k = 0; k < size.width; ++k) {
            uint64_t idx = size.width * j + k;

            front_halo[idx] = 3 * _NU * cos(k * _DX)
                                * cos(j * _DX + t)
                                * cos((global_i_start - 1) * _DX);

            back_halo[idx] = 3 * _NU * cos(k * _DX)
                               * cos(j * _DX + t)
                               * cos((global_i_start + size.depth) * _DX);
        }
    }
}

struct SolverData {
    int num_timesteps;
    field_size size;
    field porosity;
    field gamma;
    field pressure;
    field phi;
    field3 eta;
    field3 zeta;
    field3 vel;
    OutputVTK *output;
    DDecomp *ddecomp;
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
    field phi = solver_data->phi;
    field3 eta = solver_data->eta;
    field3 zeta = solver_data->zeta;
    field3 vel = solver_data->vel;
    OutputVTK *output = solver_data->output;
    DDecomp *ddecomp = solver_data->ddecomp;

    char output_file_name[32];
    sprintf(output_file_name, "solution-%.4f-%d-%d.vtk",
            _DT, get_proc_rank(ddecomp->comm_z), 0);
    output_vtk_write(output, output_file_name, thread);


    for (uint32_t t = 1; t < num_timesteps + 1; ++t) {

        momentum_solve(porosity, gamma, pressure, phi,
                       size, eta, zeta, vel, t, thread, ddecomp);

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

        pressure_solve(to_const_field3(vel), size,
                       pressure, phi, t, thread, ddecomp);

       /*
       ftype mean_pressure = 0;
       ftype mean_phi = 0;
        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            mean_pressure += pressure[i];
            mean_phi += phi[i];
        }
        mean_pressure /= field_num_points(size);
        mean_phi /= field_num_points(size);

        for (uint64_t i = 0; i < field_num_points(size); ++i) {
            pressure[i] -= mean_pressure;
            phi[i] -= mean_phi;
        }
        */

        thread_wait_on_barrier(thread);

        /*
        pressure_correct_rot(to_const_field3(vel),
                             to_const_field3(vel_old),
                             size, pressure, 0.0, t);
        */

        char output_file_name[32];
        sprintf(output_file_name, "solution-%.4f-%d-%d.vtk",
                _DT, get_proc_rank(ddecomp->comm_z), t);
        output_vtk_write(output, output_file_name, thread);
    }

    return 0;
}

DEF_TEST(test_convergence_time_splitting_brinkman,
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

    DDecomp *ddecomp = ddecomp_create(64, 64, 64, arena);

    field_size size = ddecomp->local_size;
    SET_DX(M_PI / (size.width - 0.5));

    for (int i = 0; i < num_samples; ++i) {
        arena_enter(arena);

        SET_DT(dt);

        field porosity = field_alloc(size, arena);
        field_fill(size, _K, porosity);

        field gamma = field_alloc(size, arena);
        field_fill(size, (_K * _NU * _DT) /
                         ((2.0 * _K + _NU * _DT) * _DX * _DX), gamma);

        field3 eta = field3_alloc_pad(size, arena);
        field3 zeta = field3_alloc_pad(size, arena);
        field3 vel = field3_alloc_pad(size, arena);

        field3_fill(size, 0, eta);
        field3_fill(size, 0, zeta);
        field3_fill(size, 0, vel);

        field pressure = field_alloc_pad(size, arena);
        field phi = field_alloc_pad(size, arena);

        /* WARNING: Fill halos with manufactured solution! */

        compute_manufactured_solution(size, 0, eta);
        compute_manufactured_solution(size, 0, zeta);
        compute_manufactured_solution(size, 0, vel);
        compute_manufactured_pressure(size, 0.5, pressure);
        field_fill(size, 0, phi);

        field manufactured_pressure = field_alloc(size, arena);
        compute_manufactured_pressure(size, 0.5, manufactured_pressure);

        OutputVTK *output = output_vtk_create(size, _DX, arena);

        output_vtk_attach_field(output, pressure, "pressure", arena);
        output_vtk_attach_field(output, manufactured_pressure,
                                "man_pressure", arena);
        output_vtk_attach_field3(output, to_const_field3(vel),
                                 "velocity", arena);


        if (get_proc_rank(MPI_COMM_WORLD) == 0) {
            printf("%d/%d\n", i + 1, num_samples);
        }
        uint32_t num_timesteps = round(T / _DT);

        struct SolverData data = {
            num_timesteps,
            size,
            porosity,
            gamma,
            pressure,
            phi,
            eta,
            zeta,
            vel,
            output,
            ddecomp
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

    if (get_proc_rank(MPI_COMM_WORLD) == 0) {
        printf("%.8f %e   --  %e   --\n", dts[0], v_errors[0], p_errors[0]);
        for (int i = 1; i < num_samples; ++i) {
            printf("%.8f %e %5.2f %e %5.2f\n",
                   dts[i], v_errors[i], v_orders[i],
                   p_errors[i], p_orders[i]);
        }
    }

    thread_array_destroy(t_array);

    arena_exit(arena);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    ArenaAllocator arena;
    arena_init(&arena, 1ul << 33);

    RUN_TEST(test_convergence_time_splitting_brinkman, &arena, 3, 1);

    arena_destroy(&arena);

    MPI_Finalize();
}
