#include "../include/utils.h"
#include <omp.h>
#include <errno.h>

/* Compilar usando el Makefile */

int main(void)
{
    int N_THREADS = 0, N_PART = 0, BINS = 0, steps[50], retake = 0, dump = 0;
    unsigned int Ntandas = 0u;
    char inputFilename[255], saveFilename[255];
    double DT = 0.0, M = 0.0, sigmaL = 0.0;

    int X0 = 1;
    char filename[32];

    double d = 1.0e-72, alfa = 1.5E+42; // alfa = 1.4E+43
    int evolution = 0;
    //    double pmin075 = 7.20843424240426E-020, pmax075 = 1.2818610191887E-017;
    double pmin = 3.0E-026, pmax = 3.0E-023;

    char data_filename[] = "datos.in";

    load_parameters_from_file(data_filename, &N_PART, &BINS, &DT, &M, &N_THREADS, &Ntandas, steps, inputFilename,
                              saveFilename, &retake, &dump, &sigmaL);

    omp_set_dynamic(0);
    omp_set_num_threads(N_THREADS);

    uint32_t base_seed = (uint32_t)time(NULL);
    const char *seed_text = getenv("KINETICBOX_SEED");
    if (seed_text != NULL)
    {
        char *end;
        errno = 0;
        unsigned long parsed = strtoul(seed_text, &end, 10);
        if (errno != 0 || seed_text == end || *end != '\0' || (seed_text[0] < '0' || seed_text[0] > '9') || parsed > UINT32_MAX)
        {
            fprintf(stderr, "Error: KINETICBOX_SEED must be an unsigned 32-bit integer\n");
            return EXIT_FAILURE;
        }
        base_seed = (uint32_t)parsed;
    }
    printf("Random seed = %u; requested threads = %d\n", (unsigned int)base_seed, N_THREADS);

    double *x = malloc(sizeof(double) * N_PART);
    double *p = malloc(sizeof(double) * N_PART);
    double *DxE = malloc(sizeof(double) * (2 * BINS + 5));
    double *DpE = malloc(sizeof(double) * (2 * BINS + 1));
    int *h = malloc(sizeof(int) * (2 * BINS + 5));
    int *g = malloc(sizeof(int) * (2 * BINS + 1));
    int *hg = malloc(sizeof(int) * (size_t)(2 * BINS + 5) * (size_t)(2 * BINS + 1));

    bool memory_allocations = check_memory_allocations(x, p, DxE, DpE, h, g, hg);
    if (!memory_allocations)
    {
        return 1;
    }

    uint32_t *seeds = malloc(sizeof(*seeds) * (size_t)N_THREADS);
    if (seeds == NULL)
    {
        fprintf(stderr, "Error: cannot allocate random states\n");
        free(x);
        free(p);
        free(DxE);
        free(DpE);
        free(h);
        free(g);
        free(hg);
        return EXIT_FAILURE;
    }
    for (int i = 0; i < N_THREADS; i++)
    {
        uint32_t seed = base_seed + UINT32_C(0x9e3779b9) * ((uint32_t)i + 1u);
        seed = (seed ^ (seed >> 16)) * UINT32_C(0x85ebca6b);
        seed = (seed ^ (seed >> 13)) * UINT32_C(0xc2b2ae35);
        seeds[i] = (seed ^ (seed >> 16));
        if (seeds[i] == 0)
            seeds[i] = UINT32_C(2463534242);
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i <= BINS << 1; i++)
    {
        double numerator = (3.0e-23 / BINS) * N_PART;
        double denominator = 5.24684E-24 * sqrt(2.0 * PI);
        double exponent = -pow(3.0e-23 * (1.0 * i / BINS - 1) / 5.24684E-24, 2) / 2;
        DpE[i] = (numerator / denominator) * exp(exponent);
    }

#pragma omp parallel for simd schedule(static)
    for (int i = 0; i <= (BINS + 2) << 1; i++)
    {
        DxE[i] = N_PART / (2.0 * BINS);
    }

    DxE[0] = 0.0;
    DxE[1] = 0.0;
    DxE[2] = DxE[2] * 0.5;
    DxE[2 * BINS + 2] = DxE[2 * BINS + 2] * 0.5;
    DxE[2 * BINS + 3] = 0.0;
    DxE[2 * BINS + 4] = 0.0;

    memset(h, 0, (2 * BINS + 5) * sizeof(int));
    memset(g, 0, (2 * BINS + 1) * sizeof(int));
    memset(hg, 0, (size_t)(2 * BINS + 5) * (size_t)(2 * BINS + 1) * sizeof(int));

    if (retake != 0)
    {
        while (X0 == 1)
        {
// initialize particles
#pragma omp parallel
            {
                uint32_t seed = seeds[omp_get_thread_num()];
#pragma omp for schedule(static)
                for (int i = 0; i < N_PART; i++)
                {
                    double randomValue = d_xorshift(&seed);
                    x[i] = randomValue * 0.5;
                }
#pragma omp for schedule(static)
                for (int i = 0; i < N_PART / 2 + N_PART % 2; i++)
                {
                    double randomValue1 = d_xorshift(&seed);
                    double randomValue2 = d_xorshift(&seed);

                    double xi1 = sqrt(-2.0 * log(randomValue1));
                    double xi2 = 2.0 * PI * randomValue2;

                    p[2 * i] = xi1 * cos(xi2) * 5.24684E-24;
                    if (2 * i + 1 < N_PART)
                        p[2 * i + 1] = xi1 * sin(xi2) * 5.24684E-24;
                }
                seeds[omp_get_thread_num()] = seed;
            }

            fill_hist(h, g, hg, x, p, N_PART, BINS);

            X0 = make_hist(h, g, hg, DxE, DpE, "X0000000.dat", BINS);
            if (X0 == 1)
            {
                printf("Falló algún chi2: X0=%1d\n", X0);
            }
        }
    }
    else
    {
        read_data(inputFilename, x, p, &evolution, N_PART);
        fill_hist(h, g, hg, x, p, N_PART, BINS);
        memset(h, 0, (2 * BINS + 5) * sizeof(int));
        memset(g, 0, (2 * BINS + 1) * sizeof(int));
        memset(hg, 0, (size_t)(2 * BINS + 5) * (size_t)(2 * BINS + 1) * sizeof(int));
    }

    int final_evolution = evolution;
    for (unsigned int j = 0; j < Ntandas; j++)
    {
        if (steps[j] > INT_MAX - final_evolution)
        {
            fprintf(stderr, "Error: step count exceeds the legacy checkpoint limit INT_MAX\n");
            exit(EXIT_FAILURE);
        }
        final_evolution += steps[j];
    }

    energy_sum(p, N_PART, evolution, M);
    printf("d=%12.9E  alfa=%12.9E\n", d, alfa);

    for (unsigned int j = 0; j < Ntandas; j++)
    {
        int invalid = 0;
#pragma omp parallel shared(x, p) reduction(| : invalid)
        {
            uint32_t seed = seeds[omp_get_thread_num()];
#pragma omp for schedule(static)
            for (int i = 0; i < N_PART; ++i)
            {
                double x_tmp = x[i];
                double p_tmp = p[i];
                int particle_invalid = 0;
                for (int step = 0; step < steps[j]; step++)
                {
                    x_tmp += p_tmp * DT / M;    // ¡OJO que p_tmp tiene un SIGNO!
                    int signop = (int)copysign(1.0, p_tmp);
                    double crossings = trunc(x_tmp + 0.5 * signop);
                    if (!isfinite(crossings) || fabs(crossings) >= (double)LONG_MAX)
                    {
                        particle_invalid = 1;
                        break;
                    }
                    long int k = (long int)crossings;
                    if (k != 0)
                    {
                        double randomValue = d_xorshift(&seed);
                        double xi1 = sqrt(-2.0 * log(randomValue));
                        randomValue = d_xorshift(&seed);
                        double xi2 = 2.0 * PI * randomValue;
                        double deltaX = sqrt((double)labs(k)) * xi1 * cos(xi2) * sigmaL;
                        if (!isfinite(deltaX))
                        {
                            particle_invalid = 1;
                            break;
                        }
                        deltaX = (fabs(deltaX) > 1.0 ? 1.0 * copysign(1.0, deltaX) : deltaX);
                        x_tmp = (k % 2 ? -1.0 : 1.0) * (x_tmp - k) + deltaX;
                        if (fabs(x_tmp) > 0.502)
                        {
                            x_tmp = 1.004 * copysign(1.0, x_tmp) - x_tmp;
                        }
                        p_tmp = fabs(p_tmp);    // <-- le saco el signo a p_tmp
                        for (long int l = 0; l < labs(k); l++)
                        {
                            double DeltaE = alfa * pow((p_tmp - pmin) * (pmax - p_tmp), 2);
                            randomValue = d_xorshift(&seed);
                            double p2 = p_tmp * p_tmp + DeltaE * (randomValue - 0.5);
                            if (!isfinite(p2) || p2 < 0)
                            {
                                particle_invalid = 1;
                                break;
                            }
                            p_tmp = sqrt(p2);
                        }
                        if (particle_invalid)
                            break;
                        p_tmp *= (k % 2 ? -1.0 : 1.0) * signop;
                    }
                }
                invalid |= particle_invalid;
                x[i] = x_tmp;
                p[i] = p_tmp;
            }
            seeds[omp_get_thread_num()] = seed;
        }
        if (invalid)
        {
            fprintf(stderr, "Error: non-finite trajectory, excessive crossings or negative momentum squared\n");
            exit(EXIT_FAILURE);
        }

        fill_hist(h, g, hg, x, p, N_PART, BINS);
        evolution += steps[j];
        snprintf(filename, sizeof(filename), "X%07d.dat", evolution);
        if (dump == 0)
        {
            save_data(saveFilename, x, p, evolution, N_PART);
        }
        make_hist(h, g, hg, DxE, DpE, filename, BINS);
        energy_sum(p, N_PART, evolution, M);
    }
    // End of Work code.

    printf("Completo evolution = %d\n", evolution);

    free(seeds);
    free(x);
    free(p);
    free(DxE);
    free(DpE);
    free(h);
    free(g);
    free(hg);

    return 0;
}

/* para graficar en el gnuplot: (sacando un archivo "hists.eps")
set terminal postscript enhanced color eps 20
set output "hists.eps"
# histograma de x  (las dos líneas siguientes alcanzan para graficar las x
dentro del gnuplot) set style fill solid 1.0 # o medio transparente: set style
fill transparent solid 0.5 noborder set key left ; set xrange[-0.5:0.5] p
'X1000000.dat' u 1:2 w boxes lc rgb "#dddddd" t 'X1000000.dat' , 'X2000000.dat'
u 1:2 w boxes lc rgb "#77ff77" t 'X2000000.dat' , 'X2000001.dat' u 1:2 w boxes
lc "#ffaaaa" t 'X2000001.dat' , 'X2000002.dat' u 1:2 w boxes lc "#dddd55" t
'X2000002.dat' , 'X2000003.dat' u 1:2 w boxes lc rgb "#ffdddd" t 'X2000003.dat'
, 'X2000008.dat' u 1:2 w boxes lc rgb "#cc44ff" t 'X2000008.dat' ,
'X2000018.dat' u 1:2 w boxes lc rgb "#888888" t 'X2000018.dat' , 'X2000028.dat'
u 1:2 w boxes lc rgb "#bbddbb" t 'X2000028.dat' , 'X2000038.dat' u 1:2 w boxes
lc rgb "#ffee00" t 'X2000038.dat' , 'X2000048.dat' u 1:2 w boxes lc rgb
"#8844ff" t 'X2000048.dat' , 'X2000058.dat' u 1:2 w boxes lc rgb "#cceeff" t
'X2000058.dat' , 'X2000068.dat' u 1:2 w boxes lc rgb "#44bb44" t 'X2000068.dat'
, 'X2000078.dat' u 1:2 w boxes lc rgb "#99ee77" t 'X2000078.dat' ,
'X2000088.dat' u 1:2 w boxes lc rgb "#ffdd66" t 'X2000088.dat' , 'X2000098.dat'
u 1:2 w boxes lc rgb "#4444ff" t 'X2000098.dat' # histograma de p  (las dos
líneas siguientes alcanzan para graficar las p dentro del gnuplot) set key left
; set xrange[-3e-23:3e-23] p 'X0000500.dat' u 3:4 w boxes lc rgb "#dddddd" t
'X0000500.dat' , 'X0001000.dat' u 3:4 w boxes lc rgb "#77ff77" t 'X0001000.dat'
, 'X0002000.dat' u 3:4 w boxes lc "#ffaaaa" t 'X0002000.dat' , 'X0005000.dat' u
3:4 w boxes lc "#dddd55" t 'X0005000.dat' , 'X0010000.dat' u 3:4 w boxes lc rgb
"#ffdddd" t 'X0010000.dat' , 'X0020000.dat' u 3:4 w boxes lc rgb "#cc44ff" t
'X0020000.dat' , 'X0050000.dat' u 3:4 w boxes lc rgb "#888888" t 'X0050000.dat'
, 'X0100000.dat' u 3:4 w boxes lc rgb "#bbddbb" t 'X0100000.dat' ,
'X0200000.dat' u 3:4 w boxes lc rgb "#ffee00" t 'X0200000.dat' , 'X0500000.dat'
u 3:4 w boxes lc rgb "#8844ff" t 'X0500000.dat' , 'X0995000.dat' u 3:4 w boxes
lc rgb "#cceeff" t 'X0995000.dat' , 'X0999000.dat' u 3:4 w boxes lc rgb
"#44bb44" t 'X0999000.dat' , 'X0999500.dat' u 3:4 w boxes lc rgb "#99ee77" t
'X0999500.dat' , 'X1000000.dat' u 3:4 w boxes lc rgb "#ffdd66" t 'X1000000.dat'
, 'X2000000.dat' u 3:4 w boxes lc rgb "#4444ff" t 'X2000000.dat' set terminal qt


p 'X0000001.dat' u 1:2 w boxes lc rgb "#dddddd" t 'X0000001.dat' ,
'X0000100.dat' u 1:2 w boxes lc rgb "#77ff77" t 'X0000100.dat' , 'X0001000.dat'
u 1:2 w boxes lc "#ffaaaa" t 'X0001000.dat' , 'X0001200.dat' u 1:2 w boxes lc
"#dddd55" t 'X0001200.dat' , 'X0001400.dat' u 1:2 w boxes lc rgb "#ffdddd" t
'X0001400.dat' , 'X0001500.dat' u 1:2 w boxes lc rgb "#cc44ff" t 'X0001500.dat'
, 'X0001600.dat' u 1:2 w boxes lc rgb "#888888" t 'X0001600.dat' ,
'X0001700.dat' u 1:2 w boxes lc rgb "#bbddbb" t 'X0001700.dat' , 'X0001800.dat'
u 1:2 w boxes lc rgb "#ffee00" t 'X0001800.dat' , 'X0001900.dat' u 1:2 w boxes
lc rgb "#8844ff" t 'X0001900.dat' , 'X0002000.dat' u 1:2 w boxes lc rgb
"#cceeff" t 'X0002000.dat'

*/