#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <interactions.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)

void cholesky(double *mm, int chain_length, double *cc, int t);
float ran1(long *idum);
float gasdev(long *idum);
long initRan();
float stochastic();

struct Bead {
    double x, y, z;
    double fx, fy, fz;
};

int main(int argc, const char *argv[]) {

    // variables
    int N, Nbb, Nsc; // chain length, number of backbone beads, number of sidechain beads
    int tmax, k, t, i, j, n, nb, PMFcount; 
    double epsilon, kappa, dt, D, L; // epsilon, kappa, timestep, diffusion constant, box length
    double dx, dy, dz, Fs, rr, r, p, r6, ratio, coeff, c, c2, c3, PMF_force, PMF_energy, PMFave_force, PMFave_energy; // variables for force calculation

    // file input
    FILE *inputfile;
    inputfile = fopen("input.txt", "r");
    fscanf(inputfile, "epsilon = %lf\n", &epsilon);
    fscanf(inputfile, "kappa = %lf\n", &kappa);
    fscanf(inputfile, "dt = %lf\n", &dt);
    fscanf(inputfile, "tmax = %d\n", &tmax);
    fscanf(inputfile, "Nbb = %d\n", &Nbb);
    fscanf(inputfile, "Nsc = %d\n", &Nsc);
    fscanf(inputfile, "D = %lf\n", &D);
    fclose(inputfile);

    N = Nbb * (Nsc + 1); // total number of beads (nbb * nsc + nbb)
    L = 2.0 * (double) Nbb; // box length
    PMFcount = 0; // counter for PMF calculation
    PMFave_force = 0.0; // average force for PMF calculation
    PMFave_energy = 0.0; // average energy for PMF calculation

    time_t the_time;
    time(&the_time);

    // open files
    char *str1 = malloc(sizeof(char) * 30);
    char *str2 = malloc(sizeof(char) * 30);
    char *str3 = malloc(sizeof(char) * 30);

    sprintf(str1, "RFD%d_%d_edot%d.xyz", (int) (D * 10), Nbb, Nsc);
    sprintf(str2, "RFD%d_%d_edot%d.dat", (int) (D * 10), Nbb, Nsc);
    sprintf(str3, "stats.txt");

    FILE *output_file, *output, *stats;
    output_file = fopen(str1, "w");

    // fprintf(output_file, ctime(&the_time));
    fprintf(output_file, "epsilon = %lf\n", epsilon);
    fprintf(output_file, "kappa = %lf\n", kappa);
    fprintf(output_file, "dt = %lf\n", dt);
    fprintf(output_file, "tmax = %d\n", tmax);
    fprintf(output_file, "Nbb = %d\n", Nbb);
    fprintf(output_file, "Nsc = %d\n", Nsc);
    fprintf(output_file, "D = %lf\n", D);
    fprintf(output_file, "\n");
    fprintf(output_file, "\n");
    fclose(output_file);

    struct Bead *chain_a = malloc(N * sizeof(struct Bead)); 
    struct Bead *chain_b = malloc(N * sizeof(struct Bead));

    double *random_a = calloc(N, sizeof(double)); // random numbers for chain a
    double *random_b = calloc(N, sizeof(double)); // random numbers for chain b

    p = sqrt(2.0 * dt); // sqrt(2 * dt) for random number generation
    long *idum = malloc(sizeof(long)); // random number seed

}

float random(long *idum) {
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; --j)
        {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    if (*idum < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

float gasdev(long *idum) {

    float random(long *idum);
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;

    if (*idum < 0)
        iset = 0;
    if (iset == 0) {
        do {
            v1 = 2.0 * random(idum) - 1.0;
            v2 = 2.0 * random(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }

}

void cholesky(double *mm, int chain_length, double *cc, int t) {

    int i, j, k;
    int n = chain_length * 3;
    double sum;
    int error = 0;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            cc[i + n * j] = 0;
    }
    for (i = 0; i < n; ++i) {
        for (j = i; j < n; ++j) {
            sum = mm[i + n * j];
            for (k = i - 1; k >= 0; --k)
                sum -= cc[i + n * k] * cc[j + n * k];
            if (i == j) {
                if (sum <= 0) {
                    error = 1;
                    break;
                }
                cc[i + n * i] = sqrt(sum);
            } else
                cc[j + n * i] = sum / cc[i + n * i];
        }
    }
    if (error == 1) {
        printf("Cholesky decomposition failed at t = %d\n", t);
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                cc[i + n * j] = 0;
                if (i == j)
                    cc[i + n * j] = 1;
            }
        }
    }

}

long init_random() {
    time_t seconds;
    time(&seconds);
    return (long) seconds;
}
