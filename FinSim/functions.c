#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
//remember2check if all libraries needed and then remove this cmment, las he arrastrado!


// Exchange halo cells between neighboring MPI ranks
void Halo_Update(double *T, int local_n, int rank, int size, MPI_Comm comm) {
    MPI_Request reqs[4];
    MPI_Status stats[4];
    int req_cnt = 0;

    // Send leftmost actual to left neighbor, receive from right neighbor
    if (rank > 0) {
        MPI_Irecv(&T[0], 1, MPI_DOUBLE, rank-1, 0, comm, &reqs[req_cnt++]);
        MPI_Isend(&T[1], 1, MPI_DOUBLE, rank-1, 1, comm, &reqs[req_cnt++]);
    }
    // Send rightmost actual to right neighbor, receive from left neighbor
    if (rank < size-1) {
        MPI_Irecv(&T[local_n+1], 1, MPI_DOUBLE, rank+1, 1, comm, &reqs[req_cnt++]);
        MPI_Isend(&T[local_n], 1, MPI_DOUBLE, rank+1, 0, comm, &reqs[req_cnt++]);
    }
    //cuidado que esto no tarde mucho al esperar?
    //medir tiempo?
    // Wait for all non-blocking ops to complete
    MPI_Waitall(req_cnt, reqs, stats);
}

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
    //i tried to make a linear interpolation instead, didn't work this is simpler and works
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

// Iterative time loop (Jacobi iterations) w/ convergence check
void Time_Loop(double *T_old, double *T_new, int local_n, int rank, int size, double dx, double h, double P, double T_inf, double k, int max_iter, double tol, MPI_Comm comm) {
    double a = k / (dx*dx);
    double b = h * P;
    double global_diff;
    int iter = 0;
//this is important, functions in a library can call other functions in the same libary!
    do {
        //exchange halo values with neighbors
        Halo_Update(T_old, local_n, rank, size, comm);
        //apply BC(updates ghost cells)
        BC_Update(T_old, local_n, rank, size, T_old[0], T_inf, h, dx, k);
        // Jacobi update
        Solve_Step(T_old, T_new, local_n, a, b, T_inf);

        //local max difference for convergence
        double local_diff = 0.0;
        for (int i = 1; i <= local_n; i++) {
            double diff = fabs(T_new[i] - T_old[i]);
            if (diff > local_diff) local_diff = diff;
        }
        // Global reduction to find worst difference
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, comm);

        //copy new values into old array for next iteration
        for (int i = 1; i <= local_n; i++) {
            T_old[i] = T_new[i];
        }

        iter++;


        //STATUS UPDATE FOR CONSOLE SANITY CHECK!!
        //helps to see if converging or reaching max_iter
        if (rank == 0 && (iter % 10 == 0 || global_diff < tol)) {
            if (global_diff < tol) {
            printf("Converged at iteration %d with residual %e\n", iter, global_diff);
            } else {
                printf("Iteration %d: Residual = %e\n", iter, global_diff);
            }
            fflush(stdout);
            }

        if (iter >= max_iter) break;
    } while (global_diff > tol);
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


// Gather local data and write output CSV (done by rank 0)
void Print(double *T_local, int local_n, int rank, int size, double dx, int N, double k) {
    int i;
    //Ccalculate heat flux thru conduction at each node: flux = -k * dT/dx
    //get left ghost value (base or from neighbor)
    double T_left = T_local[0];
    //assemble local arrays to send
    int send_count = (rank == 0) ? local_n+1 : local_n;
    double *T_send = malloc(send_count * sizeof(double));
    double *flux_send = malloc(send_count * sizeof(double));
    if (!T_send || !flux_send) {
        fprintf(stderr, "Rank %d: malloc failed in Print\n", rank);
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    //obtain position and flux for local nodes
    // Global start index (number of nodes before this rank's segment)
    int offset = 0;
    for (int r = 0; r < rank; r++) {
        offset += (N/size) + (r < N%size ? 1 : 0);
    }
    // If rank 0, include base at global index = 0
    if (rank == 0) {
        double x0 = 0.0;
        double T0 = T_local[0];  // base (ghost) = T_hot
        double T1 = T_local[1];
        double flux0 = -k * (T1 - T0) / dx;  // flux into fin at base
        T_send[0] = T0;
        flux_send[0] = flux0;
        // Actual nodes 1..local_n
        for (i = 1; i <= local_n; i++) {
            double Ti = T_local[i];
            double Ti_left = (i == 1) ? T0 : T_local[i-1];
            double flux_i = -k * (Ti - Ti_left) / dx;
            T_send[i] = Ti;
            flux_send[i] = flux_i;
        }
    } else {
        // Non-zero ranks: do not include overlapping node from previous rank
        for (i = 1; i <= local_n; i++) {
            double Ti = T_local[i];
            double Ti_left = (i == 1) ? T_local[0] : T_local[i-1];
            double flux_i = -k * (Ti - Ti_left) / dx;
            T_send[i-1] = Ti;
            flux_send[i-1] = flux_i;
        }
    }

    //prep gather
    int *recvcounts = NULL, *displs = NULL;
    int total_points = 0;
    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }
    int sendcnt = send_count;
    MPI_Gather(&sendcnt, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;
        total_points = recvcounts[0];
        for (int r = 1; r < size; r++) {
            displs[r] = displs[r-1] + recvcounts[r-1];
            total_points += recvcounts[r];
        }
    }

    // Gather temperature and flux to rank 0
    double *T_global = NULL, *flux_global = NULL;
    if (rank == 0) {
        T_global = malloc(total_points * sizeof(double));
        flux_global = malloc(total_points * sizeof(double));
    }
    MPI_Gatherv(T_send, sendcnt, MPI_DOUBLE, T_global, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(flux_send, sendcnt, MPI_DOUBLE, flux_global, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //rank 0 writes to data.csv
    //pay attention to this? may be a bottleneck?
    //is there a way to write data to same file in orderly fashion in parallel?
    //maybe write multiple files then merge them?
    //if too problematic runtime of this section could be measured

    if (rank == 0) {
        FILE *fout = fopen("data.csv", "w");
        if (!fout) { perror("data.csv"); MPI_Abort(MPI_COMM_WORLD,1); }
        fprintf(fout, "position,T,heat_flux\n");
        for (int idx = 0; idx < total_points; idx++) {
            double x = dx * idx;
            fprintf(fout, "%.6f,%.6f,%.6f\n", x, T_global[idx], flux_global[idx]);
        }
        fclose(fout);
    }

    //freememory
    free(T_send);
    free(flux_send);
    if (rank == 0) {
        free(T_global);
        free(flux_global);
        free(recvcounts);
        free(displs);
    }
}

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


