#include "functions.h"
#include<stdio.h>
// Function implementations
void printGreeting(const char* name) {
    printf("HAL9000: I'm sorry %s, I'm afraid I can't do that.\n", name);
}

int add(int num1, int num2) {
    return num1 + num2;
}










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









/*
This is the implementation file where you actually implement your functions. 
You include your header file to access its function prototypes and ensure that you don't have any type mismatches.
*/