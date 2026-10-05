#include "../include/utils.h"
#include <ctype.h>
#include <errno.h>
#include <unistd.h>

static double d_rand()
{
    srand(time(NULL));
    return (double)rand() / (double)RAND_MAX;
}

static void fail(const char *message)
{
    fprintf(stderr, "Error: %s\n", message);
    exit(EXIT_FAILURE);
}

static bool next_line(FILE *file, char line[4096])
{
    while (fgets(line, 4096, file) != NULL)
    {
        if (strchr(line, '\n') == NULL && !feof(file))
            fail("configuration line is too long");
        char *comment = strchr(line, '#');
        if (comment != NULL)
            *comment = '\0';
        char *text = line;
        while (isspace((unsigned char)*text))
            text++;
        if (*text != '\0')
        {
            memmove(line, text, strlen(text) + 1);
            return true;
        }
    }
    if (ferror(file))
        fail("cannot read configuration");
    return false;
}

static char *value_after_colon(char *line)
{
    char *value = strchr(line, ':');
    if (value == NULL)
        fail("missing ':' in configuration");
    return value + 1;
}

static char *read_value(FILE *file, char line[4096])
{
    if (!next_line(file, line))
        fail("incomplete configuration");
    return value_after_colon(line);
}

static bool at_end(const char *text)
{
    while (isspace((unsigned char)*text))
        text++;
    return *text == '\0';
}

static int parse_int(const char *text)
{
    char *end;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno != 0 || text == end || !at_end(end) || value < 1 || value > INT_MAX)
        fail("expected an integer between 1 and INT_MAX");
    return (int)value;
}

static double parse_double(const char *text)
{
    char *end;
    errno = 0;
    double value = strtod(text, &end);
    if (errno != 0 || text == end || !at_end(end) || !isfinite(value))
        fail("expected a finite number");
    return value;
}

static int parse_switch(const char *text, char filename[255])
{
    char option[8], extra;
    if (sscanf(text, " %7s %254s %c", option, filename, &extra) != 2)
        fail("expected yes/no and a filename of at most 254 bytes without spaces");
    if (strcmp(option, "sí") == 0 || strcmp(option, "si") == 0 || strcmp(option, "yes") == 0)
        return 0;
    if (strcmp(option, "no") == 0)
        return 1;
    fail("expected sí, si, yes or no");
    return 1;
}

void load_parameters_from_file(char filename[], int *N_PART, int *BINS, double *DT, double *M, int *N_THREADS,
                               unsigned int *Ntandas, int steps[], char inputFilename[], char saveFilename[],
                               int *retake, int *dump, double *sigmaL)
{
    char line[4096];
    FILE *inputFile = fopen(filename, "r");
    if (inputFile == NULL)
        fail("cannot open configuration file");
    *N_PART = parse_int(read_value(inputFile, line));
    *BINS = parse_int(read_value(inputFile, line));
    *DT = parse_double(read_value(inputFile, line));
    *M = parse_double(read_value(inputFile, line));
    *N_THREADS = parse_int(read_value(inputFile, line));
    if (!at_end(read_value(inputFile, line)))
        fail("batch steps must start on the line after their heading");
    *Ntandas = 0;
    while (next_line(inputFile, line) && strchr(line, ':') == NULL)
    {
        char *text = line;
        while (!at_end(text))
        {
            char *end;
            errno = 0;
            long value = strtol(text, &end, 10);
            if (errno != 0 || text == end || value < 1 || value > INT_MAX ||
                (*end != '\0' && !isspace((unsigned char)*end)))
                fail("invalid batch step count");
            if (*Ntandas >= 50)
                fail("at most 50 batches are supported");
            steps[(*Ntandas)++] = (int)value;
            text = end;
        }
    }
    if (*Ntandas == 0 || strchr(line, ':') == NULL)
        fail("missing batches or checkpoint settings");
    *retake = parse_switch(value_after_colon(line), inputFilename);
    *dump = parse_switch(read_value(inputFile, line), saveFilename);
    *sigmaL = parse_double(read_value(inputFile, line));
    if (next_line(inputFile, line))
        fail("unexpected trailing configuration data");
    if (fclose(inputFile) != 0)
        fail("cannot close configuration file");
    if (*BINS < 2 || *BINS > 500 || *DT <= 0 || *M <= 0 || *sigmaL < 0)
        fail("BINS must be 2..500; time step and mass positive; wall uncertainty nonnegative");
    size_t nx = (size_t)(2 * *BINS + 5), np = (size_t)(2 * *BINS + 1);
    if ((size_t)*N_PART > SIZE_MAX / sizeof(double) || nx > SIZE_MAX / sizeof(double) ||
        np > SIZE_MAX / sizeof(double) || nx > SIZE_MAX / sizeof(int) / np ||
        (size_t)*N_THREADS > SIZE_MAX / sizeof(uint32_t))
        fail("requested arrays are too large");
    printf("sigma(L) = %le\n", *sigmaL);
}

void read_data(char filename[], double *x, double *p, int *evolution, int N_PART)
{
    FILE *readFile = fopen(filename, "rb");
    if (readFile == NULL)
        fail("cannot open checkpoint");
    bool valid = fread(evolution, sizeof(*evolution), 1, readFile) == 1 &&
                 fread(x, sizeof(*x), (size_t)N_PART, readFile) == (size_t)N_PART &&
                 fread(p, sizeof(*p), (size_t)N_PART, readFile) == (size_t)N_PART;
    if (!valid || fgetc(readFile) != EOF || ferror(readFile) || *evolution < 0)
        fail("checkpoint is truncated, has the wrong particle count or contains an invalid step count");
    if (fclose(readFile) != 0)
        fail("cannot close checkpoint");
}

void energy_sum(double *p, int N_PART, int evolution, double M)
{
    double sumEnergy = 0;
#pragma omp parallel for reduction(+ : sumEnergy) schedule(static)
    for (int i = 0; i < N_PART; i++)
    {
        sumEnergy += p[i] * p[i];
    }
    printf("N° de pasos %6d\tEnergía total = %12.9E\n", evolution, sumEnergy / (2 * M));
}

void save_data(char filename[], double *x, double *p, int evolution, int N_PART)
{
    char temporary[300];
    if (snprintf(temporary, sizeof(temporary), "%s.tmp.XXXXXX", filename) >= (int)sizeof(temporary))
        fail("checkpoint path is too long");
    int fd = mkstemp(temporary);
    if (fd < 0)
        fail("cannot create temporary checkpoint");
    FILE *saveFile = fdopen(fd, "wb");
    if (saveFile == NULL)
    {
        close(fd);
        remove(temporary);
        fail("cannot open temporary checkpoint stream");
    }
    bool valid = fwrite(&evolution, sizeof(evolution), 1, saveFile) == 1 &&
                 fwrite(x, sizeof(*x), (size_t)N_PART, saveFile) == (size_t)N_PART;
    int Npmod = (0 * N_PART) / (1 << 21);
    if (evolution % 1000000 == 0 && Npmod > 0)
    {
        double f = 0.7071; // fraccion de p+ que queda en p+'
        double *sqrtp2 = malloc(sizeof(double) * Npmod);
        int np = 0;
        int i0 = d_rand() * N_PART;
        int i = i0;
        while ((np < Npmod) && (i < N_PART))
        {
            if (fabs(p[i]) > (2.43 + 0.3 * np / Npmod) * 5.24684E-24)
            {
                sqrtp2[np] = sqrt(1.0 - f * f) * p[i];
                np++;
                p[i] *= f;
            }
            i++;
        }
        i = 0;
        while ((np < Npmod) && (i < i0))
        {
            if (fabs(p[i]) > (2.43 + 0.3 * np / Npmod) * 5.24684E-24)
            {
                sqrtp2[np] = sqrt(1.0 - f * f) * p[i];
                np++;
                p[i] *= f;
            }
            i++;
        }
        printf("np=%d   (2.43-2.73)sigma\n", np);
        np = 0; // repartimos 0.5*E+ en 2 partes iguales
        while ((np < Npmod) && (i < N_PART))
        {
            int signopr = copysign(1.0, sqrtp2[np]);
            if ((signopr * p[i] > 0) && (fabs(p[i]) > 0.15 * 5.24684E-24) && (fabs(p[i]) < 0.9 * 5.24684E-24))
            {
                p[i] = sqrt(p[i] * p[i] + sqrtp2[np] * sqrtp2[np] / 2.0);
                np++;
            }
            i++;
        }
        i = 0;
        while (np < Npmod)
        {
            int signopr = copysign(1.0, sqrtp2[np]);
            if ((signopr * p[i] > 0) && (fabs(p[i]) > 0.15 * 5.24684E-24) && (fabs(p[i]) < 0.9 * 5.24684E-24))
            {
                p[i] = sqrt(p[i] * p[i] + sqrtp2[np] * sqrtp2[np] / 2.0);
                np++;
            }
            i++;
        }
        np = 0; // otra vez la busqueda de p chicos, porque repartimos la otra
                // mitad de E+
        while ((np < Npmod) && (i < N_PART))
        {
            int signopr = copysign(1.0, sqrtp2[np]);
            if ((signopr * p[i] > 0) && (fabs(p[i]) > 0.15 * 5.24684E-24) && (fabs(p[i]) < 0.9 * 5.24684E-24))
            {
                p[i] = sqrt(p[i] * p[i] + sqrtp2[np] * sqrtp2[np] / 2.0);
                np++;
            }
            i++;
        }
        i = 0;
        while (np < Npmod)
        {
            int signopr = copysign(1.0, sqrtp2[np]);
            if ((signopr * p[i] > 0) && (fabs(p[i]) > 0.15 * 5.24684E-24) && (fabs(p[i]) < 0.9 * 5.24684E-24))
            {
                p[i] = sqrt(p[i] * p[i] + sqrtp2[np] * sqrtp2[np] / 2.0);
                np++;
            }
            i++;
        }
        free(sqrtp2);
    }
    if (valid)
        valid = fwrite(p, sizeof(*p), (size_t)N_PART, saveFile) == (size_t)N_PART;
    if (fclose(saveFile) != 0)
        valid = false;
    if (!valid || rename(temporary, filename) != 0)
    {
        remove(temporary);
        fail("cannot write or replace checkpoint");
    }
}

void fill_hist(int *h, int *g, int *hg, const double *x, const double *p, int N_PART, int BINS)
{
    for (int i = 0; i < N_PART; i++)
    {
        double h_idx = floor((2.0 * x[i] + 1) * BINS + 2.5);
        double g_idx = floor((p[i] / 3.0e-23 + 1) * BINS + 0.5);
        if (!isfinite(x[i]) || !isfinite(p[i]) || !isfinite(h_idx) || !isfinite(g_idx) ||
            h_idx < 0 || h_idx >= 2 * BINS + 5 || g_idx < 0 || g_idx >= 2 * BINS + 1)
        {
            fprintf(stderr, "Error: particle %d is outside histogram support: x=%.17g, p=%.17g\n", i, x[i], p[i]);
            exit(EXIT_FAILURE);
        }
        h[(int)h_idx]++;
        g[(int)g_idx]++;
        hg[(size_t)(2 * BINS + 1) * (size_t)h_idx + (size_t)g_idx]++;
    }
}

int make_hist(int *h, int *g, int *hg, double *DxE, double *DpE, const char *filename, int BINS)
{
    double chi2x = 0.0, chi2xr = 0.0, chi2p = 0.0, chiIp = 0.0, chiPp = 0.0, chiIx = 0.0, chiPx = 0.0;

    int borders = (int)(BORDES * (BINS / 500.0));

    if (strcmp(filename, "X0000000.dat") == 0)
    {
#pragma omp parallel for reduction(+ : chi2x) schedule(static)
        for (int i = BINS + 2; i <= 2 * BINS + 2; i++)
        {
            double expected = (i == BINS + 2) ? DxE[i] : 2.0 * DxE[i];
            chi2x += pow(h[i] - expected, 2) / expected;
        }
        chi2x /= BINS + 1;
    }
    else
    {
#pragma omp parallel for reduction(+ : chi2x) schedule(static)
        for (int i = 2; i <= 2 * (BINS + 1); i++)
        {
            chi2x += pow(h[i] - DxE[i], 2) / DxE[i];
        }
        chi2x = chi2x / (2.0 * BINS + 1);
        chi2xr = chi2x; // chi2xr = chi2x reducido
    }
#pragma omp parallel for reduction(+ : chi2p) schedule(static)
    for (int i = 0; i <= 2 * (BINS - borders); i++)
    {
        chi2p += pow(g[i + borders] - DpE[i + borders], 2) / DpE[i + borders];
    }
#pragma omp parallel for reduction(+ : chiIp, chiPp) schedule(static)
    for (int i = 0; i < (BINS - borders); i++)
    {
        chiIp += pow(g[i + borders] - g[2 * BINS - borders - i], 2) / DpE[i + borders];
        chiPp += pow(g[i + borders] + g[2 * BINS - borders - i] - 2.0 * DpE[i + borders], 2) / DpE[i + borders];
    }
#pragma omp parallel for reduction(+ : chiIx, chiPx) schedule(static)
    for (int i = 2; i < BINS + 2; i++)
    {
        chiIx += pow(h[i] - h[2 * BINS + 4 - i], 2) / DxE[i];
        chiPx += pow(h[i] + h[2 * BINS + 4 - i] - 2.0 * DxE[i], 2) / DxE[i];
    }
    chiIx = chiIx / (2.0 * BINS);
    chiPx = chiPx / (2.0 * BINS);
    chi2p = chi2p / (2.0 * (BINS - borders) + 1);
    chiIp = chiIp / (2.0 * (BINS - borders));
    chiPp = chiPp / (2.0 * (BINS - borders));

    FILE *hist = fopen(filename, "w");
    if (hist == NULL)
        fail("cannot open histogram output");
    fprintf(hist,
            "#   x    poblacion       p      poblacion    chi2x =%9.6f  chi2xr "
            "=%9.6f  chiIx =%9.6f  chiPx =%9.6f  chi2p =%9.6f  chiIp =%9.6f  "
            "chiPp =%9.6f\n",
            chi2x, chi2xr, chiIx, chiPx, chi2p, chiIp, chiPp);
    fprintf(hist, "%8.5f %6d %24.12E %6d\n", -0.5 - 1.0 / BINS, h[0], -3.0e-23, 0);
    fprintf(hist, "%8.5f %6d %24.12E %6d\n", -0.5 - 0.5 / BINS, h[1], -3.0e-23, 0);
    for (int i = 0; i <= BINS << 1; i++)
    {
        fprintf(hist, "%8.5f %6d %24.12E %6d\n", (0.5 * i / BINS - 0.5), h[i + 2], (3.0e-23 * i / BINS - 3.0e-23),
                g[i]);
    }
    fprintf(hist, "%8.5f %6d %24.12E %6d\n", 0.5 + 0.5 / BINS, h[2 * BINS + 3], 3.0e-23, 0);
    fprintf(hist, "%8.5f %6d %24.12E %6d\n", 0.5 + 1.0 / BINS, h[2 * BINS + 4], 3.0e-23, 0);

    bool failed = ferror(hist) != 0;
    if (fclose(hist) != 0 || failed)
        fail("cannot write histogram output");

    memset(h, 0, (2 * BINS + 5) * sizeof(int));
    memset(g, 0, (2 * BINS + 1) * sizeof(int));
    memset(hg, 0, (size_t)(2 * BINS + 5) * (size_t)(2 * BINS + 1) * sizeof(int));

    return 0; // avisa que se cumplió la condición sobre los chi2
}

bool check_memory_allocations(double *x, double *p, double *DxE, double *DpE, int *h, int *g, int *hg)
{
    if (x && p && DxE && DpE && h && g && hg)
        return true;
    fprintf(stderr, "Error: cannot allocate simulation arrays\n");
    free(x);
    free(p);
    free(DxE);
    free(DpE);
    free(h);
    free(g);
    free(hg);
    return false;
}