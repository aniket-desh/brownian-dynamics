#include <math.h>
#include <stddef.h>

#include "interactions.h"

static inline double min_image(double dz, double box_length) {
    if (dz > 0.5 * box_length) {
        dz -= box_length;
    } else if (dz < -0.5 * box_length) {
        dz += box_length;
    }
    return dz;
}

static void apply_spring_forces(Bead *chain, const SimulationParams *params, int total_beads, ForceEnergyStats *stats) {
    const int Nbb = params->Nbb;
    const int Nsc = params->Nsc;
    const double kappa = params->kappa;
    const double L = params->box_length;

    for (int i = 0; i < total_beads - Nbb; ++i) {
        int current = Nbb + i;
        int partner = (i % Nsc == 0) ? (i / Nsc) : current - 1;

        double dx = chain[current].x - chain[partner].x;
        double dy = chain[current].y - chain[partner].y;
        double dz = min_image(chain[current].z - chain[partner].z, L);

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);
        if (r == 0.0) {
            continue;
        }

        double displacement = r - 2.0;
        double force_scalar = -kappa * displacement / r;
        double fx = force_scalar * dx;
        double fy = force_scalar * dy;
        double fz = force_scalar * dz;

        chain[current].fx += fx;
        chain[current].fy += fy;
        chain[current].fz += fz;

        chain[partner].fx -= fx;
        chain[partner].fy -= fy;
        chain[partner].fz -= fz;

        stats->spring_energy += 0.5 * kappa * displacement * displacement;
    }
}

static void apply_lj_forces(Bead *chain, const SimulationParams *params, int total_beads, ForceEnergyStats *stats) {
    const double epsilon = params->epsilon;
    const double cutoff_sq = 5.04; // squared cutoff
    const double L = params->box_length;

    for (int i = 0; i < total_beads; ++i) {
        for (int j = i + 1; j < total_beads; ++j) {
            double dx = chain[i].x - chain[j].x;
            double dy = chain[i].y - chain[j].y;
            double dz = min_image(chain[i].z - chain[j].z, L);

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= cutoff_sq || r2 == 0.0) {
                continue;
            }

            double inv_r2 = 1.0 / r2;
            double ratio = 4.0 * inv_r2;
            double r6 = ratio * ratio * ratio;
            double force_coeff = (epsilon * inv_r2) * (48.0 * r6 * r6 - 24.0 * r6);
            double fx = force_coeff * dx;
            double fy = force_coeff * dy;
            double fz = force_coeff * dz;

            chain[i].fx += fx;
            chain[i].fy += fy;
            chain[i].fz += fz;

            chain[j].fx -= fx;
            chain[j].fy -= fy;
            chain[j].fz -= fz;

            stats->total_lj_energy += 4.0 * epsilon * (r6 * r6 - r6) + epsilon;
        }
    }
}

static void apply_cross_forces(Bead *chain_a, Bead *chain_b, const SimulationParams *params, int total_beads, ForceEnergyStats *stats) {
    const double epsilon = params->epsilon;
    const double cutoff_sq = 5.04; // squared cutoff
    const double L = params->box_length;

    for (int i = 0; i < total_beads; ++i) {
        for (int j = 0; j < total_beads; ++j) {
            double dx = chain_a[i].x - chain_b[j].x;
            double dy = chain_a[i].y - chain_b[j].y;
            double dz = min_image(chain_a[i].z - chain_b[j].z, L);

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= cutoff_sq || r2 == 0.0) {
                continue;
            }

            double inv_r2 = 1.0 / r2;
            double ratio = 4.0 * inv_r2;
            double r6 = ratio * ratio * ratio;
            double force_coeff = (epsilon * inv_r2) * (48.0 * r6 * r6 - 24.0 * r6);
            double fx = force_coeff * dx;
            double fy = force_coeff * dy;
            double fz = force_coeff * dz;

            chain_a[i].fx += fx;
            chain_a[i].fy += fy;
            chain_a[i].fz += fz;

            chain_b[j].fx -= fx;
            chain_b[j].fy -= fy;
            chain_b[j].fz -= fz;

            stats->pmf_force -= fx;
            double potential = 4.0 * epsilon * (r6 * r6 - r6) + epsilon;
            stats->pmf_energy += potential;
            stats->total_lj_energy += potential;
        }
    }
}

ForceEnergyStats compute_forces(Bead *chain_a, Bead *chain_b, const SimulationParams *params, int total_beads) {
    ForceEnergyStats stats = {0};
    apply_spring_forces(chain_a, params, total_beads, &stats);
    apply_spring_forces(chain_b, params, total_beads, &stats);

    apply_lj_forces(chain_a, params, total_beads, &stats);
    apply_lj_forces(chain_b, params, total_beads, &stats);
    apply_cross_forces(chain_a, chain_b, params, total_beads, &stats);
    return stats;
}
