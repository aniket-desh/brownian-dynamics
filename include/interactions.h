#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "simulation.h"

ForceEnergyStats compute_forces(Bead *chain_a, Bead *chain_b, const SimulationParams *params, int total_beads);

#endif // INTERACTIONS_H
