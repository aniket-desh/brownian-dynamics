// external var declaration
extern int N, Nbb, Nsc; // chain length, number of backbone beads, number of sidechain beads
extern int tmax, k, t, i, j, n, nb, PMFcount; 
extern double epsilon, kappa, dt, D, L; // epsilon, kappa, timestep, diffusion constant, box length
extern double dx, dy, dz, Fs, rr, r, p, r6, ratio, coeff, c, c2, c3, PMF_force, PMF_energy, PMFave_force, PMFave_energy; // variables for force calculation

// spring force
for (i = 0; i < (N - Nbb); i++) {
    if (i % Nsc > 0) {

    }
}