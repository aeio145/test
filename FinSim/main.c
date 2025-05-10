#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"


int main(int argc, char *argv[]) {
    // Material properties of the fin
    //definition in manual input comment
    double k, rho, cp, alpha;   
    // Select material and physical properties will be automatically retrieved "Aluminum", "Copper", "Gold", "Iron", "Steel".
    definematerial("Aluminum", &k, &rho, &cp, &alpha);

    //Uncomment for manual input of a prperty!
    //double k = 200.0;       // thermal conductivity [W/m*K]
    //double rho = 2700.0;    // density [kg/m^3]
    //double cp = 900.0;      // specific heat [J/kg*K]
    //double alpha = k/(rho*cp);  // thermal diffusivity [m^2/s]
    
    // Physical and numerical parameters
    double T_hot = 600.0;   // base temperature [ºC]
    double T_cold = 25.0;   // ambient temperature [ºC]
    double h = 100.0;       // convective heat transfer coefficient [W/m^2*K]
    double L = 2.5;           //fin length [m]
    int N = 100000;            // number of cells
    double tol = 1e-6;      // convergence tolerance
    int max_iter = 1e8;  // max Jacobi iterations

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);//obtain number of processes (size is num procs and wiil be stored in values.csv)
    
    //Measuring runtime
    //based upon the slides from Arnau, PDF7 page 40
    //start local time
    double t0 = MPI_Wtime();

    if (rank == 0) {
        // Show the retrieved properties (Sanity CHeck) and parameters
        printf("Material Properties:\n");
        printf("k  = %.2f W/m·K\n", k);
        printf("rho = %.2f kg/m³\n", rho);
        printf("cp = %.2f J/kg·K\n", cp);
        printf("alpha = %e m²/s\n", alpha);
        printf("T_hot = %.2f m²/s\n", T_hot);
        printf("T_cold = %.2f m²/s\n", T_cold);
        printf("h = %.2f m²/s\n", h);
        printf("L = %.2f m²/s\n", L);
        printf("N = %d m²/s\n", N);
        printf("tol = %.2f m²/s\n", tol);
        printf("max_iter = %d m²/s\n", max_iter);
    }
    
    //broadcast parameters to all ranks
    MPI_Bcast(&T_hot, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&T_cold, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&k, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&rho, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&cp, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&L, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&N, 1, MPI_INT, 0, comm);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&max_iter, 1, MPI_INT, 0, comm);

    // Compute grid spacing (length of fin div by number of cells)
    double dx = L / N;
    // Compute local number of cells for this rank (equal partition)
    int base = N / size;   // base number of cells per rank
    int rem = N % size;    // remainder for uneven division
    int local_n;
    if (rank < rem) {
        local_n = base + 1;
    } else {
        local_n = base;
    }

    // Compute global starting index of this rank's first cell
    int offset = 0;
    for (int r = 0; r < rank; r++) {
        if (r < rem) offset += base + 1;
        else offset += base;
    }

    // allocate arrays (w/2 ghost cells: West (left) and east (right))
    double *T_old = malloc((local_n + 2) * sizeof(double));
    double *T_new = malloc((local_n + 2) * sizeof(double));
    if (!T_old || !T_new) {
        fprintf(stderr, "Rank %d: malloc failed\n", rank);
        MPI_Abort(comm,1);
    }

    //initialize the temperature field (needs starting value!) (Initial guess is all cells at ambient temperature)
    //if the number of iterations is toolow it will be evident by not having a proper 
    //exponential since most cells will remain at a lower temp!(sharp curve instead of smooth) 
    Initial_C(T_old, local_n, rank, T_hot, T_cold);
    // We don't need to initialize T_new beyond the solver overwriting it.

    //exdcute iterative solver (Jacobi!)
    Time_Loop(T_old, T_new, local_n, rank, size, dx, h, 2.0, T_cold, k, max_iter, tol, comm);

    // Final boundary ghost update for tip convection using final T_old
    BC_Update(T_old, local_n, rank, size, T_hot, T_cold, h, dx, k);

    // Print results (gather and output by rank 0)
    Print(T_old, local_n, rank, size, dx, N, k);

    // free memory
    free(T_old);
    free(T_new);

    // Measure runtime
    double avg_time;
    measure_time(t0, comm, &avg_time);

    //store values to csv file
    write_values_csv(T_hot, T_cold, k, rho, cp, alpha, h, L, N, tol, max_iter, size, avg_time);
    
    MPI_Finalize();
    return 0;
}

/*

COmments to compile and run!
Add path of mpi to include path 

terminal
mpicc --showme:compile
-I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi

then in file c_cpp_properties.json
{
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
                "${workspaceFolder}/**",
                "/usr/lib/x86_64-linux-gnu/openmpi/include",
                "/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi"
            ],
            "defines": [],
            "compilerPath": "/usr/bin/mpicc",
            "cStandard": "c17",
            "cppStandard": "gnu++17",
            "intelliSenseMode": "linux-gcc-x64"
        }
    ],
    "version": 4
}

//IMPORTANT TO COMPILE ALONGSIDE functions.c ! OTherwise complains
mpicc main.c functions.c -o main
mpirun -np 4 ./main
*/