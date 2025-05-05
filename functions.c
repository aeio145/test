#include "functions.h"
#include <stdio.h>
#include <stdlib.h>// Needed for exit()
#include <stdarg.h>
#include <ctype.h>
#include <string.h>

// Function implementations


//##################################
// Working Functions
//##################################

// File creation
void createfile(const char *filename) {
    // Check if the file exists. If it does, delete it.
    if (remove(filename) == 0) {
        printf("Existing file deleted: %s\n", filename);
    } else {
        // If the file doesn't exist, remove() will fail, and this message will not be shown.
        printf("No existing file found, proceeding to create a new one: %s\n", filename);
    }

    FILE *file = fopen(filename, "w");  // Open file for writing
    if (file == NULL) {
        printf("Failed to create the file: %s\n", filename);
        return;
    }
    printf("File created successfully: %s\n", filename);
    fclose(file);
}

// Function to create the header lines of the .csv files
void writefileheader(const char *filename, ...) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    va_list args;
    va_start(args, filename);

    const char *arg;
    int first = 1;

    while ((arg = va_arg(args, const char *)) != NULL) {
        if (!first) {
            fputc(',', fp);
        }
        fputs(arg, fp);
        first = 0;
    }

    fputc('\n', fp);  // Add newline at the end

    va_end(args);
    fclose(fp);
}

// Define the material properties
void definematerial(const char *material, double *k, double *ro, double *cp, double *alpha) {
    if (strcmp(material, "Aluminum") == 0) {
        *k = 237;
        *ro = 2702;
        *cp = 903;
        *alpha = 97.1e-6;
    } else if (strcmp(material, "Copper") == 0) {
        *k = 401;
        *ro = 8933;
        *cp = 385;
        *alpha = 117e-6;
    } else if (strcmp(material, "Gold") == 0) {
        *k = 317;
        *ro = 19300;
        *cp = 129;
        *alpha = 127e-6;
    } else if (strcmp(material, "Iron") == 0) {
        *k = 80.2;
        *ro = 7870;
        *cp = 447;
        *alpha = 23.1e-6;
    } else if (strcmp(material, "Steel") == 0) {
        *k = 51.9;
        *ro = 7817;
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
    double k = W/(m * K), ro(density) = kg/(m³), cp = J/(kg * K);
  
    // Add more materials for reference if needed    
    
    */

    //double alpha = k/(ro*cp);   // This could be a simple value, defined with the other properties instead of calculated







//##################################
// Test Staging Area
//##################################


// Matrix multiplication function
void matmul(double *C, double *A, double *B, int m, int n, int q) {

// Substitute AC_MAT with the stuff on the right
#define AC_MAT(A, cols, i, j) A[(i) * (cols) + (j)]


    for (int ii = 0; ii < m; ii++) {  // Iterate over rows of A
        for (int jj = 0; jj < q; jj++) {  // Iterate over columns of B
            AC_MAT(C, q, ii, jj) = 0;  // Initialize C(ii, jj)
            for (int kk = 0; kk < n; kk++) {  // product loop
                AC_MAT(C, q, ii, jj) += AC_MAT(A, n, ii, kk) * AC_MAT(B, q, kk, jj);
            }
        }
    }
}


// Function to compute temperature difference and heat flux
void heat_cond(double T_West, double T_east, double delta_x, double k, double *cell_delta_T, double *q_1D_cell) {
*cell_delta_T = T_east - T_West;
*q_1D_cell = -k * (*cell_delta_T) / delta_x;
}



void append_data(const char *filename, double pos, double T_val, double q, double val) {
    FILE *fp = fopen(filename, "a");
    if (!fp) {
        perror("Error opening file");
        return;
    }
    fprintf(fp, "%.5f,%.5f,%.5f,%.5f\n", pos, T_val, q, val);
    fclose(fp);
}



/*
This is the implementation file where you actually implement your functions. 
You include your header file to access its function prototypes and ensure that you don't have any type mismatches.
*/

