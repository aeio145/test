#include"functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>//I use it to store runtime in the file

#define N 100  // number of cells

int main(){

    // Measure time to run at beginning, to be used for amhdal parameters later
    //Is this correct or should it be done with time in the terminal when executing whole script?
    //ASK ARNAU

    clock_t start, end;     
    double cpu_time_used;
    start = clock();  // Start time

    // Create and prepare files for data storage and further plots in Matlab
    createfile("values.csv");
    writefileheader("values.csv", "T_hot", "T_cold", "k", "rho", "cp", "alpha", NULL);
    createfile("data.csv");
    writefileheader("data.csv", "Position", "T", "Heatflux", "val", NULL);
    //##################################
    // Let's define some values
    //##################################

    // Material properties of the fin
    double k, rho, cp, alpha;
    // Select material and physical properties will be automatically retrieved "Aluminum", "Copper", "Gold", "Iron", "Steel".
    definematerial("Aluminum", &k, &rho, &cp, &alpha);
    
    // Show the retrieved properties (Sanity CHeck)
    printf("Material Properties:\n");
    printf("k  = %.2f W/m·K\n", k);
    printf("rho = %.2f kg/m³\n", rho);
    printf("cp = %.2f J/kg·K\n", cp);
    printf("alpha = %e m²/s\n", alpha);
    
    //##################################
    
    // Fin Geometry
    double L = 1e-1;    // Fin length (m)
    
    //NOT IN USE RIGHT NOW
    //double e = 1e-3;    //Fin thickness (m)
    //double W = 1e-2;    //Fin Width (m)
    
    // Cell geometry
    
    double delta_x=L/N;   // Distance between nodes (m)
    printf("delta_x  = %e m\n", delta_x);   // Sanity check

    //##################################
    
    // Setup values
    // Temperature of the fin on the hot and cold side (ºC)
    double T_hot = 50.0;
    double T_cold = 20.0;
    double T_air =18.0; //Not in use right now



    // Until this point it works, below there's a simple attempt to compute conduction

    //##################################
    //TESTING PROTOTYPE
    //##################################

    double position[N];
    double heat_flux[N];
    double T[N];

    double cell_delta_T, q_1D_cell;//  cell Temperature gradien (K) and cell heat flux (J/m²)
    /*
    heat_cond(T_hot, T_cold, delta_x, k, &cell_delta_T, &q_1D_cell);
    printf("Temperature difference: %.2f K\n", cell_delta_T);
    printf("1D Heat flux: %.2f J/m²\n", q_1D_cell);
    */

   /*
    for (int i = 0; i < N; i++) {
        // Starting node is middle of 1st cell so add half delta_x at start
        heat_cond(T_hot, T_cold, 0.5*delta_x+i*delta_x, k, &cell_delta_T, &q_1D_cell);
        position[i] = i * delta_x;
        heat_flux[i] = q_1D_cell;  // define this function
        // Store T[i]

    }*/


    // Print results













    for (int i = 0; i < N; i++) {
        double x = i * delta_x;
        position[i] = x;
    
        // Temperature at each node (linear interpolation)
        T[i] = T_hot + (T_cold - T_hot) * (x / L);
    
        // Compute heat flux between this node and next (except for the last node)
        if (i < N - 1) {
            heat_cond(T[i], T[i+1], delta_x, k, &cell_delta_T, &q_1D_cell);
            heat_flux[i] = q_1D_cell;
        } else {
            heat_flux[i] = 0.0; // Or set to same as previous if desired
        }
    
        // (Optional) Write to data.csv if needed
        append_data("data.csv", position[i], T[i], heat_flux[i], 0.0);
    }



    //testmatmul
    int m, n, q;
    //memory alocation  
 double *A = (double *)malloc(m * n * sizeof(double));//these return the pointer (address) of each matrix as a double. The malloc allocates the required bytes for each matrix and returns their address (bytes in a double times total number of values within that matrix that are double) 
 double *B = (double *)malloc(n * q * sizeof(double)); 
 double *C = (double *)malloc(m * q * sizeof(double)); 
  
 if (A == NULL || B == NULL || C == NULL) { // If the memory of any of tthem is NULL, exit program 
     printf("Memory allocation failed!\n"); 
     return 1;//This exits the main, the value helps id the issue but could also be exit(1); return0 means success, return 1 or return -1 return n indicates failure and the meaning with n 
 }

    // Initialize A and B with example values, can't understand why can't be directly declared with double A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
 double init_A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
 double init_B[] = {9, 5, 25, 6, 2, 4, 0, 2, 10}; 






























    //##################################
    //TESTING PROTOTYPE
    //##################################

    // May be stored in values.csv?

    end = clock();    // End time

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("Elapsed time: %f seconds\n", cpu_time_used);
    
    
 
    return 0;
}

//gcc main.c functions.c -o main
