//
// Created by alonso on 12-08-19.
// Changed in an new branch.
// Now we are going to use 
//

#ifndef T1_MATRIX_H
#define T1_MATRIX_H

// Importing libs
#include <iostream>
#include <string>
#include <stdexcept>
#include <cmath>
#include "matrixHandler.h"
#include <fstream>

#define fillMatrix(n, matrix) { for (int i = 0; i < n; i++) matrix[i] = 0; }

#define maxBetween(a, b) a >= b ? a : b;

#define epsilon 1e-16


template<class T>
class Matrix {
private:
    T *mat; // Store the matrix
    int n = 0; // Number of columns
    int m = 0; // Number of rows
public:
    Matrix(const Matrix &old_matrix); // Copy constructor
    Matrix(int n); // Constructor, vector like [1xn]
    Matrix(int n, int m);

    Matrix(T *mat);

    // Constructor [nxm]
    ~Matrix(); // Destructor

    int getN() const;

    int getM() const;

    // Setters and getters
    void set(int i, int j, T value);

    // Set value to (i,j) <row,column>
    T get(int i, int j) const; // Get value from (i,j) <row,column>
    void fill(T value); // Fill all the matrix with a value

    // Dimensions
    int *size() const; // Returns [n, m] as a int vector
    int length() const; // Return max dimension, usefull for vectors

    // Values
    T max() const; // Maximum value of the matrix
    T min() const; // Minimum value of the matrix
    float det() const; // Calculate determinant
    T norm() const; // Get norm of vector

    // Utilitary functions
    void disp() const; // Display matrix to console
    void save_to_file(std::string filename) const; // Save
    Matrix<T> *clone() const; // Clone the matrix

    // Booleans
    bool equals(Matrix<T> *mat) const; // Matrix has same values as another
    bool is_diag() const; // Matrix is diagonal
    bool is_identity() const; // Check if matrix is identity
    bool is_symmetric() const; // Check if matrix is symmetric
    bool is_square() const; // Check if matrix is square nxn
    bool is_vector() const; // Check if matrix is vector
    bool operator==(const Matrix<T> &matrix) const; // Equal operator
    bool operator!=(const Matrix<T> &matrix) const; // Not equal operator

    // Mathematical operation
    // Returns value of given location when asked in the form A(x,y)
    T &operator()(const unsigned &rowNo, const unsigned & colNo);
    Matrix<T> *transpose() const; // Transpose the matrix and return a new matrix
    Matrix<T> * getCofactorMatrix(int p, int q, int n) const;
    Matrix<T> * getAdjoint() const;
    Matrix<T> * inverse() const;
    Matrix<T> &operator*=(const Matrix<T> &matrix); // Multiplication with self
    Matrix<T> *operator*=(const Matrix<T> &matrix) const; // Return a new matrix
    Matrix<T> &operator*=(T a); // Multiply a constant with self
    Matrix<T> *operator*=(T a) const; // Multiply matrix by a constant and return new matrix
    Matrix<T> &operator+=(const Matrix<T> *matrix); // Add a matrix with self
    Matrix<T> *operator+(const Matrix<T> &matrix) const; // Add and return a new matrix
    Matrix<T> &operator-=(const Matrix<T> *matrix); // Substract a matrix with self
    Matrix<T> *operator-(const Matrix<T> &matrix) const; // Substract and return a new matrix
    void normalize(); // Normalize the matrix

    // Bonus (1pt), returns inverse of the Matrix

};

template<class T>
Matrix<T>::Matrix(const Matrix &old_matrix) //Copy constructor
{
    this->n = old_matrix.getN();
    this->m = old_matrix.getM();
    this->mat = new T[n*m];
    for (int i = 0; i <n; i++) {
        for (int j = 0; j < m; j++){
            this->mat[i + j * n] = old_matrix.mat[i + j * n];
        }
    }
}

template<class T>
Matrix<T>::Matrix(int n) //Vector like
{
    this->n = n;
    this->m = 1;
    this->mat = new T[n];
    fillMatrix(n, this->mat);
}

template<class T>
Matrix<T>::Matrix(int n, int m) //Matrix
{
    //assign parameters
    this->n = n;
    this->m = m;
    //initialize vector of 1 dimension
    this->mat = new T[n*m];
    fillMatrix(n*m, this->mat);
}


template<class T>
Matrix<T>::~Matrix() // Destructor
{
    std::cout << "Matrix erased" << std::endl;
}

// Setters and getters
template<class T>
int Matrix<T>::getN() const {
    return n;
}

template<class T>
int Matrix<T>::getM() const {
    return m;
}

template<class T>
void Matrix<T>::set(int i, int j, T value) // Set value to (i,j) <row,column>
{
    this->mat[j * this->n + i] = value;
}

template<class T>
T Matrix<T>::get(int i, int j) const { // Get value from (i,j) <row,column>
    return this->mat[j * this->n + i];
}

// Fill the whole matrix with a value
template<class T>
void Matrix<T>::fill(T value)
{
    for (int i = 0; i < this->n * this->m; i++){
        this->mat[i] = value;
    }
}

// Dimensions
template<class T>
int * Matrix<T>::size() const
{ // Returns [n, m] as an int vector
    int *respuesta = new int [2];
    respuesta[0] = this->n; respuesta[1] = this->m;
   // int respuesta[2] = new {this->n, this->m};
    return respuesta;
}

// Return max dimension, useful for vectors
template<class T>
int Matrix<T>::length() const
{
    return (this->n  >=  this-> m) ? n : m ;
}

// Values
template<class T>
T Matrix<T>::max() const { // Maximum value of the matrix
    T maxValue = this->get(0,0); //sets the initial value as the first pivot
    for (int i = 0; i < this->n * this->m; i++) {
        if (maxValue < this->mat[i]) {
            maxValue = this->mat[i];
        }
    }
    return maxValue;
}


template<class T>
T Matrix<T>::min() const { // Minimum value of the matrix
    T minValue = this->mat[0]; //sets the initial value as the first pivot
    for (int i = 0; i < this->n * this->m; i++){
        if (minValue > this->mat[i]) // if the minValue is bigger than the value found
            minValue = this->mat[i]; // then reassign this value.
    }
    return minValue;
}


template<class T>
float Matrix<T>::det() const { // Calculate determinant
    if (this->n != this->m) {
        throw std::logic_error("It must be a square matrix");
    }
    if (this->n == 2){ //base case
        return (float) this->mat[0]* (float) this->mat[3] - (float) this->mat[1]*(float)this->mat[2];
    }
    else if (this->n == 1){
        return (float) this->mat[0];
    }
    else { //we want to create a submatrix, and do this recursively
        return determinant(this->n, mat);
    }
}


// Function which returns cofactor of the matrix at [p][q] as answer. n is current dimension of A matrix. Returns the cofactor matrix
// Main idea extracted from GeeksForGeeks. Modified for this class, templates and for every size.
template<class T>
Matrix<T> * Matrix<T>::getCofactorMatrix(int p, int q, int n) const
{

    int i = 0, j = 0;
    Matrix *ans = new Matrix(this->n-1, this->m-1);
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                ans->set (i, j,  this->get(row, col) );
                j++;
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return ans;
}

// returns the adjoint of the given matrix. Taken from GeeksandGeeks and improved by me.
template<class T>
Matrix<T> * Matrix<T>::getAdjoint() const
{
    if ( this->getN() != this->getM() ) { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }

    int N = this->n;
    Matrix *ans = new Matrix(this->n, this->m);
    if (this->n == 1) //case size is 1
    {
        ans->set(0,0, 1);
        return ans;
    }


    int sign = 1;

    for (int i=0; i< N; i++)
    {
        for (int j=0; j< N; j++)
        { //cofactor is having wrong dimensions
            Matrix *Cofactor = this->getCofactorMatrix(i, j, N);
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            T value = (sign) * Cofactor->det(); //determinante del cofactor
            ans-> set(j, i, value ); //transposing value
            delete(Cofactor);
        }
    }
    return ans;
}


template<class T>
T Matrix<T>::norm() const{ // Get norm of vector
   // Using Frobenius's norm, for Hilbert's spaces. This works for both matrixes and vectors.
    //https://es.wikipedia.org/wiki/Norma_matricial
    T sum = 0;
    for (int i = 0; i < this->n * this->m; i++){
        sum += pow(this->mat[i],2);
    }
    return sqrt(sum);
}


// Utilitary functions

template<class T>
void Matrix<T>::disp() const { // Display matrix to console
    // printing rows
    std::cout << "/**" << std::endl;
    for (int i = 0; i < this->n; i++) {
        for (int j = 0; j < this->m; j++) {
            std::cout << this->mat[i + j * n] << "  ";
        }
        std::cout << endl;
    }
    std::cout << "**/" << std::endl;
}


template<class T>
void Matrix<T>::save_to_file(std::string filename) const {
    ofstream saved_matrix;
    saved_matrix.open(filename);
    saved_matrix << "/**";
    for (int i = 0; i < this->n*this->m; i++){
            if (i % this->n == 0) {//we have reached the beggining of a new rew
                saved_matrix << "\n";
            }
        saved_matrix << this->mat[i] << "  ";
    }
    saved_matrix << endl;
    saved_matrix << "/**" << endl;
    saved_matrix.close();
}// Save to file


template<class T>
Matrix<T> *Matrix<T>::clone() const
{
    int n = this->n;
    int m = this->m;
    Matrix *ans = new Matrix(n, m);
    //copy each value of the first matrix into the returning value

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            ans->set(i, j, this->get(i,j));
        }
    }
    return ans;

} // Clone the matrix


// Booleans
template<class T>
bool Matrix<T>::equals(Matrix<T> *mat) const{
    //first, let's see if both matrices have the same dimensions
    if ( mat->n != this->n ) //equal n
        return false;
    if ( mat->m != this->m ) //equal m
        return false;
    //now let's compare each value
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            T value = this->get(i,j) - mat->get(i,j);
            if (abs(value) > epsilon ) { //case using double or floats, precision is not reassured.
                return false;
            }
        }
    }
    return true; //if everything is equal.

} // Matrix has same values as another

template<class T>
bool Matrix<T>::is_diag() const{
    if (this->n != this->m){ //it needs to be squared.
        return false;
    }
    //check matrix
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if (i != j && this->get(i,j) != 0){ //case we are not in the diagonal part, then it must be zero
                return false;
            }
        }
    }
    return true;
} // Matrix is diagonal

template<class T>
bool Matrix<T>::is_identity() const { // Check if matrix is identity
    if (this->is_diag()==false) //a matrix needs to be diagonal to be the identity
        return false;
    for (int d = 0; d < this->n; d++){
        if (this->get(d,d) != 1 )
            return false;
    }
    return true;
}

template<class T>
bool Matrix<T>::is_symmetric() const { // Check if matrix is symmetric
    if (this->n != this->m){ //it needs to be squared.
        return false;
    }
    for (int i = 0; i < this->n; i ++){
        for (int j = 0; j < this-> m; j++){
            if (this->get(i,j) != this->get(j, i))
                return false;
        }
    }
    return true;
}

template<class T>
bool Matrix<T>::is_square() const {
    return (this->n == this->m);
}// Check if matrix is square nxn


template<class T>
bool Matrix<T>::is_vector() const{
    if (( this->n == 1 && this->m > 1) || (this->n > 1 && this-> m == 1))
        return true;
} // Check if matrix is vector



template<class T>
bool Matrix<T>::operator==(const Matrix<T> &matrix) const{
    //first, let's see if both matrices have the same dimensions
    if ( matrix.n != this->n) //equal n
        return false;
    if ( matrix.m != this->m) //equal m
        return false;
    //now let's compare each value
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (this->get(i, j) != matrix.get(i,j) ){
                return false;
            }
        }
    }
    return true; //if everything is equal.
} // Matrix has same values as another

template<class T>
bool Matrix<T>::operator!=(const Matrix<T> &matrix) const
{
    return !( this == matrix );
}// Not equal operator

template<class T>
T &Matrix<T>::operator()(const unsigned &rowNo, const unsigned & colNo){
    return this->mat[rowNo + this->n * colNo]; //(*this)(i,j)
}


// Mathematical operation
template<class T>
Matrix<T> * Matrix<T>::transpose() const {
    // Supuesto: La matriz nueva es la transpuesta. La que se entregó inicialmente queda intacta.
    Matrix *ans = new Matrix(this->n, this->m); //Create new matrix with the same dimensions
    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this->m ; j++){
            ans->set(j,i, this->get(i,j) );
        }
    }
    return ans;
} // Transpose the matrix and return a new matrix


template<class T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &matrix)
{
    if ( this->getM() != matrix.getN() ) { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }
    int limN = this->getN();
    int limM = this->getM();

    for (int i = 0; i < limN; i++) {
        for (int j = 0; j < limM; j++){
            T suma = 0;
            for (int k = 0; k < this->m; k++) {
                suma += this->get(i,k) * matrix.get(k,j);
            }
            this->set(i,j, suma);
        }
    }
    return *this;
}

// Multiplication with self
template<class T>
Matrix<T> *Matrix<T>::operator*=(const Matrix<T> &matrix) const //return a new matrix
    {
        if (this->getM() != matrix.getN()) { //case dimensions do not fit
            throw std::logic_error("Dimensiones no adecuadas");
        }
        Matrix<T> *ans = new Matrix(this->n, matrix.getM());
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < matrix.getM(); j++) {
                T suma = 0;
                for (int k = 0; k < this->m; k++) {
                    suma += this->get(i,k) * matrix.get(k, j);
                }
                ans->set(i, j, suma);
            }
        }
        return ans;
    } // Return a new matrix

template<class T>
Matrix<T> &Matrix<T>::operator*=(T a) {
    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this-> m; j++){
            this->set (i, j, a*this->get(i,j));
        }
    }
    return *this;
} // Multiply a constant with self

template<class T>
Matrix<T> *Matrix<T>::operator*=(T a) const {
    Matrix<T> *ans = new Matrix(this->n, this->m);
    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this-> m; j++) {
            ans->set(i, j, a * this->get(i, j));
        }
    }
    return ans;
} //


template<class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> *matrix){
    if ( (this->n != matrix->getN()) || (this->m != matrix->getM()) )  { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }

    for (int i = 0; i < this->n; i++) {
        for (int j = 0; j < this-> m; j++) {
            T value = this->get(i,j) + matrix->get(i,j);
            this->set(i,j, value);
        }
    }
    return *this;
} // Multiply a constant with self


template<class T>
Matrix<T> * Matrix<T>::operator+(const Matrix<T> &matrix) const {// Add and return a new matrix
    if ( (this->getN() != matrix.getN()) || (this->getM() != matrix.getM()) ) { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }

    Matrix<T>* ans = new Matrix<T>(this->n, this->m);

    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this-> m; j++){
            ans->set(i,j, this->get(i,j) + matrix.get(i,j) );
        }
    }
    return ans;
} // Add a matrix with self

template<class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> *matrix){ // Substract matrix with self
    if ( (this->n != matrix->getN()) || (this->m != matrix->getM()) )  { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }

    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this-> m; j++){
            T result = this->get(i,j) - matrix->get(i,j);
            this->set(i,j, result);
        }
    }
    return *this;
}

template<class T>
Matrix<T> *Matrix<T>::operator-(const Matrix<T> &matrix) const{
    if ( (this->n != matrix.getN()) || (this->m != matrix.getM()) )  { //case dimensions do not fit
        throw std::logic_error("Dimensiones no adecuadas");
    }

    Matrix<T> *ans = new Matrix(this->n, this->m);

    for (int i = 0; i < this->n; i++){
        for (int j = 0; j < this-> m; j++){
            ans->set(i,j, this->get(i,j) - matrix.get(i,j) );
        }
    }
    return ans;

} // Substract and return a new matrix



// Function to calculate and store inverse, returns false if
// matrix is singular
template<class T>
Matrix<T> *Matrix<T>::inverse() const
{
    // Find determinant of A[][]
    auto* inv = new Matrix(this->n, this->m);
    T det = this->det();
    if (det == 0 || this->is_square() == false)  {
        throw std::logic_error("there is not inverse matrix");
    }
    // Find adjoint
    auto* adj = this->getAdjoint();
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i< this->n; i++) {
        for (int j = 0; j < this->n; j++) {
            float adjent_value = (float) adj->get(i,j);
            float prevalue = (float) adj->get(i,j) / (float) det;
            T value = (float) adj->get(i,j) / (float) det;
            inv->set(i,j, value);
        }
    }
    delete(adj);
    return inv;
}

template<class T>
void Matrix<T>::normalize(){
    T max = this->max();
    for (int i = 0; i < this->n * this->m; i++){
        this->mat[i] = this->mat[i]/max;
    }
} // Normalize the matrix



#endif //T1_MATRIX_H
