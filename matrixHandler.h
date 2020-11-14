//
// Created by alonso on 14-08-19.
//

#ifndef T1_MATRIXHANDLER_H
#define T1_MATRIXHANDLER_H

#include<iostream>
#include<math.h>

using namespace std;

// Function taken from https://www.tutorialspoint.com/cplusplus-program-to-compute-determinant-of-a-matrix
// Thanks to Chandu yadav for publishing this solution
// Modified by Alonso Utreras with Templates.
template<class T>
T determinant(int n, T *matrix) { //return the determinant of the given matrix
    //of dimension n up to size 10 (more than that would be too hard for the computer)
    T **new_matrix = (T **) malloc (sizeof(T*) * n);
    for (int i = 0; i < n; i++)
        new_matrix[i] = (T*) malloc (sizeof(T) * n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            new_matrix[i][j] = matrix[i + j*n];
        }
    }
    T det = 0;
    T *SubM = new T[n*n];
    if (n == 2) // Base case
        return ((new_matrix[0][0] * new_matrix[1][1]) - (new_matrix[1][0] * new_matrix[0][1]));
    else { //using cofactor's method.
        for (int x = 0; x < n; x++) { //going from row to row of the big matrix
            int subi = 0; //the step of the submatrix
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x) { // we do not get the submatrix of the row and column of the cofactor
                        continue; // so jump this row and column
                    }
                    //if the pivot is not in the col and row of the cofactor, add to submatrix
                    SubM[subi + subj*(n-1)] = new_matrix[i][j];
                    subj++; //count the following column for the submatrix
                }
                subi++; //count the following row for the submatrix
            } //the determinant is calculated recursively
            det += (pow(-1, x) * new_matrix[0][x] * determinant( n - 1, SubM));
        }
    }
    return det;
}


#endif //T1_MATRIXHANDLER_H
