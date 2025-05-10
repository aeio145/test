#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <mpi.h>

//equal spacing between them allows to see the comment on top of the function when hovering with the mouse in vscode, quick description useful!
//remember to be consistent as possible!

// Update ghost (halo) cells by exchanging w/ neighbors
void Halo_Update(double *T, int local_n, int rank, int size, MPI_Comm comm);

// Apply boundary conditions: fixed base and convective tip
void BC_Update(double *T, int local_n, int rank, int size, double T_hot, double T_inf, double h, double dx, double k);

// Create mesh or compute local cell count (not heavily used here)
void Mesh(int N, double L, int size, int rank, int *local_n, int *offset);

// Initialize temperature field
void Initial_C(double *T, int local_n, int rank, double T_hot, double T_inf);

// perform one Jacobi solve step
void Solve_Step(double *T_old, double *T_new, int local_n, double a, double b, double T_inf);

// Perform iterative loop until convergence
void Time_Loop(double *T_old, double *T_new, int local_n, int rank, int size, double dx, double h, double P, double T_inf, double k, int max_iter, double tol, MPI_Comm comm);

// Define the material properties
void definematerial(const char *material, double *k, double *rho, double *cp, double *alpha);

// THis is used at the end for the processes times
int compare_doubles(const void *a, const void *b);

// Gather and print results
void Print(double *T_local, int local_n, int rank, int size, double dx, int N, double k);

// Measure runtime
void measure_time(double t0, MPI_Comm comm, double *avg_time);

// Write values to csv file
void write_values_csv(double T_hot, double T_cold, double k, double rho, double cp, double alpha, double h, double L, int N, double tol, int max_iter, int size, double avg_time);



#endif // FUNCTIONS_H
