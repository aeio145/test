#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//##################################
//##################################
//##################################
// Function Prototypes
//##################################
//##################################
//##################################




//##################################
// Working Functions
//##################################


// Print a predefined message using an input string
void printGreeting(const char* name);

// Add two integer numbers
int add(int num1, int num2);

// File creation
void createfile(const char *filename);

// Function to create the header lines of the .csv files
void writefileheader(const char *filename, ...);

// Define material properties
void definematerial(const char *material, double *k, double *ro, double *cp, double *alpha);


//##################################
// Test Staging Area
//##################################


// Matrix multiplication function
void matmul(double *C, double *A, double *B, int m, int n, int q);



void heat_cond(double T_West, double T_east, double delta_x, double k, double *cell_delta_T, double *q_1D_cell);












/*
// Function to save values to a CSV file
void savevalues(const char *filename, ...);
*/




































#endif //FUNCTIONS_H

/*
This is just a header file. It includes function declarations and some basic guard macros to prevent multiple inclusions.
*/















void append_data(const char *filename, double pos, double T_val, double q, double val);