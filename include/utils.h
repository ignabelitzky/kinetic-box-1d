#ifndef UTILS_H
#define UTILS_H

#include "constants.h"
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static inline double d_xorshift(uint32_t *state)
{
    uint32_t x = *state ? *state : UINT32_C(2463534242);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return ((double)x + 0.5) / 4294967296.0;
}

void load_parameters_from_file(char filename[], int *N_PART, int *BINS, double *DT, double *M, int *N_THREADS,
                               unsigned int *Ntandas, int steps[], char inputFilename[], char saveFilename[],
                               int *retake, int *dump, double *sigmaL);

void read_data(char filename[], double *x, double *p, int *evolution, int N_PART);

void save_data(char filename[], double *x, double *p, int evolution, int N_PART);

void energy_sum(double *p, int N_PART, int evolution, double M);

void fill_hist(int *h, int *g, int *hg, const double *x, const double *p, int N_PART, int BINS);

int make_hist(int *h, int *g, int *hg, double *DxE, double *DpE, const char *filename, int BINS);

bool check_memory_allocations(double *x, double *p, double *DxE, double *DpE, int *h, int *g, int *hg);

#endif