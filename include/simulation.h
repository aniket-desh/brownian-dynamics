#ifndef SIMULATION_H
#define SIMULATION_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    double x, y, z;
    double fx, fy, fz;
} Bead;

typedef struct {
    double epsilon;
    double kappa;
    double dt;
    int tmax;
    int Nbb;
    int Nsc;
    double D;
    int output_stride;
    int stats_stride;
    double box_length;
} SimulationParams;

typedef struct {
    double pmf_force;
    double pmf_energy;
    double total_lj_energy;
    double spring_energy;
} ForceEnergyStats;

void initialize_chains(Bead *chain_a, Bead *chain_b, const SimulationParams *params, int total_beads);
void zero_forces(Bead *beads, size_t count);
void apply_brownian_kick(Bead *beads, size_t count, double dt, double sqrt_2dt, uint32_t *rng_state, double box_length);
double chain_center_distance(const Bead *chain_a, const Bead *chain_b, size_t count, double box_length);
double radius_of_gyration(const Bead *chain, size_t count, double box_length);

#endif // SIMULATION_H
