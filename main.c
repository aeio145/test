#include"functions.h"
#include<stdio.h>
#include <stdlib.h>





// Substitute AC_MAT with the stuff on the right
#define AC_MAT(A, cols, i, j) A[(i) * (cols) + (j)]










int main(){
    printGreeting("Dave");

    int sum = add(5, 7);
    printf("The sum is: %d\n", sum);








    int m = 3, n = 3, q = 3;  // Matrix dims
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














    matmul(C, A, B, m, n, q);//Call matrix mult function

    // Print matrix C
    printf("Matrix C:\n");
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < q; jj++) {
            printf("%5.1f ", AC_MAT(C, q, ii, jj));
        }
        printf("\n");
    }

    // Free allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}

//gcc main.c functions.c -o main
