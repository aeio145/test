#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#include <stdarg.h>
#include <ctype.h>
#include <string.h>


// BC UPdate
//left end (rank 0) T[0] = T_hot. Right end (rank size-1) compute ghost via convection.
void BC_Update(double *T, int local_n, int rank, int size, double T_hot, double T_inf, double h, double dx, double k) {
    if (rank == 0) {
        T[0] = T_hot;  // Fixed base temperature (ghost cell for left)
    }
    if (rank == size-1) {
        // Convective tip: set right ghost based on current tip value
        double T_tip = T[local_n];
        T[local_n+1] = T_tip - (h * dx / k) * (T_tip - T_inf);
    }
}


// temperature init condition: starting temperature is ambient in all points except the base (set by BC check previous function comment)
void Initial_C(double *T, int local_n, int rank, double T_hot, double T_inf) {
    // Set starting temp to ambient temperature
    for (int i = 0; i < local_n+2; i++) {
        T[i] = T_inf;
    }
    //base ghost cell init (if rank 0)
    if (rank == 0) {
        T[0] = T_hot;
    }
}

//calculate one Jacobi iteraton step on local block
void Solve_Step(double *T_old, double *T_new, int local_n, double a, double b, double T_inf) {
    for (int i = 1; i <= local_n; i++) {
        double left = T_old[i-1];
        double right = T_old[i+1];
        T_new[i] = (a*(left + right) + b*T_inf) / (2.0*a + b);
    }
}


// Define the material properties
void definematerial(const char *material, double *k, double *rho, double *cp, double *alpha) {
    if (strcmp(material, "Aluminum") == 0) {
        *k = 237;
        *rho = 2702;
        *cp = 903;
        *alpha = 97.1e-6;
    } else if (strcmp(material, "Copper") == 0) {
        *k = 401;
        *rho = 8933;
        *cp = 385;
        *alpha = 117e-6;
    } else if (strcmp(material, "Gold") == 0) {
        *k = 317;
        *rho = 19300;
        *cp = 129;
        *alpha = 127e-6;
    } else if (strcmp(material, "Iron") == 0) {
        *k = 80.2;
        *rho = 7870;
        *cp = 447;
        *alpha = 23.1e-6;
    } else if (strcmp(material, "Steel") == 0) {
        *k = 51.9;
        *rho = 7817;
        *cp = 446;
        *alpha = 14.9e-6;
    }
    else {
        // Unknown material: print error and terminate program
        fprintf(stderr, "Error: Unknown material '%s'\n", material);
        exit(EXIT_FAILURE);
    }
}

/*
    // Quick Reference table for generic values

    In this table material properties of different materials are specified in the following manner:
    // Material: 
    double k (thermal conductivity) = Watt/(meter * Kelvin), ro (density) = Kilogram/(cubic meter), cp (specific heat) = Joule/(Kilogram * Kelvin);
    // Material: 
    double k = W/(m · K), rho(density) = kg/(m³), cp = J/(kg * K);
  
    // Add more materials for reference if needed  
    //double alpha = k/(ro·cp);   // This could be calculated, defined with the other properties instead of calculated  
    
    */

void measure_time(double t0, MPI_Comm comm, double *avg_time) {
    double t1 = MPI_Wtime();
    double elapsed = t1 - t0;

    double min_time, max_time, sum_time;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Reduce(&elapsed, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&elapsed, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    if (rank == 0) {
        *avg_time = sum_time / size;
        printf("Min time: %.6f\n", min_time);
        printf("Max time: %.6f\n", max_time);
        printf("Avg time: %.6f\n", *avg_time);
    }
}

void write_values_csv(double T_hot, double T_cold, double k, double rho, double cp, double alpha, double h, double L, int N, double tol, int max_iter, int size, double avg_time) {
    FILE *fval = fopen("values.csv", "w");
    if (!fval) { perror("values.csv"); MPI_Abort(MPI_COMM_WORLD, 1); }
    fprintf(fval, "T_hot,T_cold,k,rho,cp,alpha,h,L,N,tol,max_iter,size,avg_time\n");
    fprintf(fval, "%.2f,%.2f,%.2f,%.2f,%.2f,%.5e,%.2f,%.2f,%d,%.5e,%d,%d,%.6f\n",
            T_hot, T_cold, k, rho, cp, alpha, h, L, N, tol, max_iter, size, avg_time);
    fclose(fval);
}
