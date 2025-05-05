#include"functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define N 100  // number of cells

int main(){

    // Create and prepare files for data storage and further plots in Matlab
    createfile("values.csv");
    writefileheader("values.csv", "T_hot", "T_cold", "k", "ro", "cp", "alpha", NULL);
    createfile("data.csv");
    writefileheader("data.csv", "Position", "T", "Heatflux", "val", NULL);
    //##################################
    // Let's define some values
    //##################################

    // Material properties of the fin
    double k, ro, cp, alpha;
    // Select material and physical properties will be automatically retrieved "Aluminum", "Copper", "Gold", "Iron", "Steel".
    definematerial("Aluminum", &k, &ro, &cp, &alpha);
    
    // Show the retrieved properties (Sanity CHeck)
    printf("Material Properties:\n");
    printf("k  = %.2f W/m·K\n", k);
    printf("ro = %.2f kg/m³\n", ro);
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



























    //##################################
    //TESTING PROTOTYPE
    //##################################
    
    
 
    return 0;
}

//gcc main.c functions.c -o main
