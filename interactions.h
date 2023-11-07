// external var declaration
// extern int N, Nbb, Nsc; // chain length, number of backbone beads, number of sidechain beads
// extern int tmax, k, t, i, j, n, nb, PMFcount; 
// extern double epsilon, kappa, dt, D, L; // epsilon, kappa, timestep, diffusion constant, box length
// extern double dx, dy, dz, Fs, rr, r, p, r6, ratio, coeff, c, c2, c3, PMF_force, PMF_energy, PMFave_force, PMFave_energy; // variables for force calculation

// spring force
for (int i = 0; i < (N - Nbb); i++) {
    if (i % Nsc > 0) {
        dx = chain_a[Nbb + i].x - chain_a[Nbb + i - 1].x;
        dy = chain_a[Nbb + i].y - chain_a[Nbb + i - 1].y;
        dz = chain_a[Nbb + i].z - chain_a[Nbb + i - 1].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }

        rr = dx * dx + dy * dy + dz * dz;
        r = sqrt(rr);
        Fs = -kappa * (r - 2.0);

        chain_a[Nbb + i].fx += Fs * dx / r;
        chain_a[Nbb + i].fy += Fs * dy / r;
        chain_a[Nbb + i].fz += Fs * dz / r;
        chain_a[Nbb + i - 1].fx += -Fs * dx / r;
        chain_a[Nbb + i - 1].fy += -Fs * dy / r;
        chain_a[Nbb + i - 1].fz += -Fs * dz / r;

    } else if (i % Nsc == 0) {
        nb = i / Nsc;
        dx = chain_a[Nbb + i].x - chain_a[nb].x;
        dy = chain_a[Nbb + i].y - chain_a[nb].y;
        dz = chain_a[Nbb + i].z - chain_a[nb].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }

        rr = dx * dx + dy * dy + dz * dz;
        r = sqrt(rr);
        Fs = -kappa * (r - 2.0);

        chain_a[Nbb + i].fx += Fs * dx / r;
        chain_a[Nbb + i].fy += Fs * dy / r;
        chain_a[Nbb + i].fz += Fs * dz / r;
    }
}

for (i = 0; i < (N - Nbb); i++) {

    if (i % Nsc > 0) {
        dx = chain_b[Nbb + i].x - chain_b[Nbb + i - 1].x;
        dy = chain_b[Nbb + i].y - chain_b[Nbb + i - 1].y;
        dz = chain_b[Nbb + i].z - chain_b[Nbb + i - 1].z;

        if (dz > L / 2) {
            dz -= L;
        }
        else if (dz < -L / 2) {
            dz += L;
        }

        rr = dx * dx + dy * dy + dz * dz;
        r = sqrt(rr);
        Fs = -kappa * (r - 2.0);

        chain_b[Nbb + i].fx += Fs * dx / r;
        chain_b[Nbb + i].fy += Fs * dy / r;
        chain_b[Nbb + i].fz += Fs * dz / r;
        chain_b[Nbb + i - 1].fx += -Fs * dx / r;
        chain_b[Nbb + i - 1].fy += -Fs * dy / r;
        chain_b[Nbb + i - 1].fz += -Fs * dz / r;
    } else if (i % Nsc == 0) {
        nb = i / Nsc;
        dx = chain_b[Nbb + i].x - chain_b[nb].x;
        dy = chain_b[Nbb + i].y - chain_b[nb].y;
        dz = chain_b[Nbb + i].z - chain_b[nb].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }

        rr = dx * dx + dy * dy + dz * dz;
        r = sqrt(rr);
        Fs = -kappa * (r - 2.0);

        chain_b[Nbb + i].fx += Fs * dx / r;
        chain_b[Nbb + i].fy += Fs * dy / r;
        chain_b[Nbb + i].fz += Fs * dz / r;
    }
}

// LJ force
for (int i = 0; i < N; i++) {

    for (int j = i + 1; j < N; j++) {
        dx = chain_a[i].x - chain_a[j].x;
        dy = chain_a[i].y - chain_a[j].y;
        dz = chain_a[i].z - chain_a[j].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }
        rr = dx * dx + dy * dy + dz * dz;
        if (rr < 5.04) {
            ratio = 4.00 / rr;
            r6 = ratio * ratio * ratio;
            coeff = (epsilon / rr) * (48.0 * r6 * r6 - 24.0 * r6);
            chain_a[i].fx += coeff * dx;
            chain_a[i].fy += coeff * dy;
            chain_a[i].fz += coeff * dz;
            chain_a[j].fx -= coeff * dx;
            chain_a[j].fy -= coeff * dy;
            chain_a[j].fz -= coeff * dz;
        }
    }
}

for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
        dx = chain_b[i].x - chain_b[j].x;
        dy = chain_b[i].y - chain_b[j].y;
        dz = chain_b[i].z - chain_b[j].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }
        rr = dx * dx + dy * dy + dz * dz;

        if (rr < 5.04) {
            ratio = 4.00 / rr;
            r6 = ratio * ratio * ratio;
            coeff = (epsilon / rr) * (48.0 * r6 * r6 - 24.0 * r6);
            chain_b[i].fx += coeff * dx;
            chain_b[i].fy += coeff * dy;
            chain_b[i].fz += coeff * dz;
            chain_b[j].fx -= coeff * dx;
            chain_b[j].fy -= coeff * dy;
            chain_b[j].fz -= coeff * dz;
        }
    }
}

for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
        dx = chain_a[i].x - chain_b[j].x;
        dy = chain_a[i].y - chain_b[j].y;
        dz = chain_a[i].z - chain_b[j].z;

        if (dz > L / 2) {
            dz -= L;
        }
        if (dz < -L / 2) {
            dz += L;
        }

        rr = dx * dx + dy * dy + dz * dz;
        if (rr < 5.04) {
            ratio = 4.00 / rr;
            r6 = ratio * ratio * ratio;
            coeff = (epsilon / rr) * (48.0 * r6 * r6 - 24.0 * r6);
            chain_a[i].fx += coeff * dx;
            chain_a[i].fy += coeff * dy;
            chain_a[i].fz += coeff * dz;
            chain_b[j].fx -= coeff * dx;
            chain_b[j].fy -= coeff * dy;
            chain_b[j].fz -= coeff * dz;

            PMF_force -= coeff * dx;
            PMF_energy += (4 * epsilon) * (r6 * r6 - r6) + epsilon;
        }
    }
}