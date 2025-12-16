#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "interactions.h"
#include "simulation.h"

static inline uint32_t lcg_rand(uint32_t *state) {
    *state = (*state * 1664525u + 1013904223u);
    return *state;
}

static double random_uniform(uint32_t *state) {
    // Use the full 32-bit output range to sample [0, 1) without clipping
    return lcg_rand(state) / 4294967296.0;  // 2^32
}

static double random_normal(uint32_t *state) {
    static int has_spare = 0;
    static double spare;

    if (has_spare) {
        has_spare = 0;
        return spare;
    }

    double u, v, s;
    do {
        u = 2.0 * random_uniform(state) - 1.0;
        v = 2.0 * random_uniform(state) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    has_spare = 1;
    return u * s;
}

static void read_input(const char *path, SimulationParams *params) {
    FILE *file = fopen(path, "r");
    if (!file) {
        perror("Unable to open input file");
        exit(EXIT_FAILURE);
    }

    // defaults for optional parameters
    params->output_stride = 1000;
    params->stats_stride = 1000;

    char line[128];
    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "epsilon = %lf", &params->epsilon) == 1) continue;
        if (sscanf(line, "kappa = %lf", &params->kappa) == 1) continue;
        if (sscanf(line, "dt = %lf", &params->dt) == 1) continue;
        if (sscanf(line, "tmax = %d", &params->tmax) == 1) continue;
        if (sscanf(line, "Nbb = %d", &params->Nbb) == 1) continue;
        if (sscanf(line, "Nsc = %d", &params->Nsc) == 1) continue;
        if (sscanf(line, "D = %lf", &params->D) == 1) continue;
        if (sscanf(line, "output_stride = %d", &params->output_stride) == 1) continue;
        if (sscanf(line, "stats_stride = %d", &params->stats_stride) == 1) continue;
    }

    fclose(file);
    params->box_length = 2.0 * (double)params->Nbb;
}

void initialize_chains(Bead *chain_a, Bead *chain_b, const SimulationParams *params, int total_beads) {
    const int Nbb = params->Nbb;
    const int Nsc = params->Nsc;
    const double D = params->D;

    for (int i = 0; i < Nbb; ++i) {
        chain_a[i].x = 0.0;
        chain_a[i].y = 0.0;
        chain_a[i].z = 1.0 + 2.0 * (double)i;

        chain_b[i].x = D;
        chain_b[i].y = 0.0;
        chain_b[i].z = 1.0 + 2.0 * (double)i;
    }

    int sidechain_count = total_beads - Nbb;
    for (int i = 0; i < sidechain_count; ++i) {
        int offset = i + Nbb;
        int col = i % Nsc;
        int row = i / Nsc;

        chain_a[offset].x = -2.0 - 2.0 * (double)col;
        chain_a[offset].y = 0.0;
        chain_a[offset].z = 1.0 + 2.0 * (double)row;

        chain_b[offset].x = D + 2.0 + 2.0 * (double)col;
        chain_b[offset].y = 0.0;
        chain_b[offset].z = 1.0 + 2.0 * (double)row;
    }
}

void zero_forces(Bead *beads, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        beads[i].fx = 0.0;
        beads[i].fy = 0.0;
        beads[i].fz = 0.0;
    }
}

static inline double wrap_z(double value, double box_length) {
    if (value >= box_length) {
        value -= box_length;
    } else if (value < 0.0) {
        value += box_length;
    }
    return value;
}

void apply_brownian_kick(Bead *beads, size_t count, double dt, double sqrt_2dt, uint32_t *rng_state, double box_length) {
    for (size_t i = 0; i < count; ++i) {
        double noise_x = random_normal(rng_state);
        double noise_y = random_normal(rng_state);
        double noise_z = random_normal(rng_state);

        beads[i].x += dt * beads[i].fx + sqrt_2dt * noise_x;
        beads[i].y += dt * beads[i].fy + sqrt_2dt * noise_y;
        beads[i].z += dt * beads[i].fz + sqrt_2dt * noise_z;

        beads[i].z = wrap_z(beads[i].z, box_length);
    }
}

double chain_center_distance(const Bead *chain_a, const Bead *chain_b, size_t count, double box_length) {
    double ax = 0.0, ay = 0.0;
    double bx = 0.0, by = 0.0;

    // Use first bead as reference for periodic-aware z averaging
    double ref_z_a = chain_a[0].z;
    double ref_z_b = chain_b[0].z;
    double sum_dz_a = 0.0, sum_dz_b = 0.0;

    for (size_t i = 0; i < count; ++i) {
        ax += chain_a[i].x; ay += chain_a[i].y;
        bx += chain_b[i].x; by += chain_b[i].y;

        double dz_a = chain_a[i].z - ref_z_a;
        if (dz_a > 0.5 * box_length) dz_a -= box_length;
        if (dz_a < -0.5 * box_length) dz_a += box_length;
        sum_dz_a += dz_a;

        double dz_b = chain_b[i].z - ref_z_b;
        if (dz_b > 0.5 * box_length) dz_b -= box_length;
        if (dz_b < -0.5 * box_length) dz_b += box_length;
        sum_dz_b += dz_b;
    }

    ax /= (double)count; ay /= (double)count;
    bx /= (double)count; by /= (double)count;
    double az = ref_z_a + sum_dz_a / (double)count;
    double bz = ref_z_b + sum_dz_b / (double)count;

    double dz = az - bz;
    if (dz > 0.5 * box_length) dz -= box_length;
    if (dz < -0.5 * box_length) dz += box_length;

    double dx = ax - bx;
    double dy = ay - by;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double radius_of_gyration(const Bead *chain, size_t count, double box_length) {
    // Use first bead as reference for periodic-aware z averaging
    double cx = 0.0, cy = 0.0;
    double ref_z = chain[0].z;
    double sum_dz = 0.0;

    for (size_t i = 0; i < count; ++i) {
        cx += chain[i].x;
        cy += chain[i].y;

        double dz = chain[i].z - ref_z;
        if (dz > 0.5 * box_length) dz -= box_length;
        if (dz < -0.5 * box_length) dz += box_length;
        sum_dz += dz;
    }
    cx /= (double)count;
    cy /= (double)count;
    double cz = ref_z + sum_dz / (double)count;

    double rg2 = 0.0;
    for (size_t i = 0; i < count; ++i) {
        double dx = chain[i].x - cx;
        double dy = chain[i].y - cy;
        double dz = chain[i].z - cz;

        if (dz > 0.5 * box_length) dz -= box_length;
        if (dz < -0.5 * box_length) dz += box_length;

        rg2 += dx * dx + dy * dy + dz * dz;
    }

    return sqrt(rg2 / (double)count);
}

static void write_stats_header(const char *path) {
    FILE *file = fopen(path, "w");
    if (!file) {
        perror("Unable to open stats file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "step,time,pmf_force,pmf_energy,center_distance,lj_energy,spring_energy,rg_a,rg_b\n");
    fclose(file);
}

static void append_stats(const char *path, int step, double time, double pmf_force, double pmf_energy, double center_distance, double lj_energy, double spring_energy, double rg_a, double rg_b) {
    FILE *file = fopen(path, "a");
    if (!file) {
        perror("Unable to write to stats file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", step, time, pmf_force, pmf_energy, center_distance, lj_energy, spring_energy, rg_a, rg_b);
    fclose(file);
}

static void append_xyz(const char *path, int step, double time, const Bead *chain_a, const Bead *chain_b, int total_beads) {
    FILE *file = fopen(path, "a");
    if (!file) {
        perror("Unable to open trajectory file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%d\nstep=%d time=%.6f\n", total_beads * 2, step, time);
    for (int i = 0; i < total_beads; ++i) {
        fprintf(file, "A\t%.6f\t%.6f\t%.6f\n", chain_a[i].x, chain_a[i].y, chain_a[i].z);
    }
    for (int i = 0; i < total_beads; ++i) {
        fprintf(file, "B\t%.6f\t%.6f\t%.6f\n", chain_b[i].x, chain_b[i].y, chain_b[i].z);
    }
    fclose(file);
}

int main(void) {
    SimulationParams params = {0};
    read_input("input.txt", &params);

    const int total_beads = params.Nbb * (params.Nsc + 1);
    const double sqrt_2dt = sqrt(2.0 * params.dt);
    uint32_t rng_state = (uint32_t)time(NULL);

    Bead *chain_a = calloc((size_t)total_beads, sizeof(Bead));
    Bead *chain_b = calloc((size_t)total_beads, sizeof(Bead));
    if (!chain_a || !chain_b) {
        fprintf(stderr, "Failed to allocate bead arrays\n");
        return EXIT_FAILURE;
    }

    initialize_chains(chain_a, chain_b, &params, total_beads);

    char trajectory_path[64];
    char stats_path[64];
    snprintf(trajectory_path, sizeof(trajectory_path), "trajectory_D%d_Nbb%d_Nsc%d.xyz", (int)(params.D * 10), params.Nbb, params.Nsc);
    snprintf(stats_path, sizeof(stats_path), "stats_D%d_Nbb%d_Nsc%d.csv", (int)(params.D * 10), params.Nbb, params.Nsc);
    write_stats_header(stats_path);

    for (int step = 0; step < params.tmax; ++step) {
        zero_forces(chain_a, (size_t)total_beads);
        zero_forces(chain_b, (size_t)total_beads);

        ForceEnergyStats energy = compute_forces(chain_a, chain_b, &params, total_beads);
        apply_brownian_kick(chain_a, (size_t)total_beads, params.dt, sqrt_2dt, &rng_state, params.box_length);
        apply_brownian_kick(chain_b, (size_t)total_beads, params.dt, sqrt_2dt, &rng_state, params.box_length);

        int next_step = step + 1;
        if (next_step % params.stats_stride == 0) {
            double time = next_step * params.dt;
            double center_distance = chain_center_distance(chain_a, chain_b, (size_t)total_beads, params.box_length);
            double rg_a = radius_of_gyration(chain_a, (size_t)total_beads, params.box_length);
            double rg_b = radius_of_gyration(chain_b, (size_t)total_beads, params.box_length);
            append_stats(stats_path, next_step, time, energy.pmf_force, energy.pmf_energy, center_distance, energy.total_lj_energy, energy.spring_energy, rg_a, rg_b);
        }

        if (next_step % params.output_stride == 0) {
            double time = next_step * params.dt;
            append_xyz(trajectory_path, next_step, time, chain_a, chain_b, total_beads);
        }
    }

    free(chain_a);
    free(chain_b);
    return EXIT_SUCCESS;
}
