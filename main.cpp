//
// Created by alonso on 13-08-19.
//


#include <cassert>
#include <iostream>
#include "matrix.h"

//Test module
int main() {
    //let's initialize some matrices
    auto *A = new Matrix<int>(3, 3);
    auto *B = new Matrix<int>(3, 3);
    auto *v = new Matrix<int>(3, 1); //vector
    auto *u = new Matrix<int>(1, 3); //another vector
    //Constructors are working fine

    A->fill(9);
    for (int i = 0; i < A->getN(); i++) {
        for (int j = 0; j < A->getM(); j++) {
            assert (A->get(i, j) == 9);
        }
    }
    B->fill(7);
    A->fill(9);
    for (int i = 0; i < B->getN(); i++) {
        for (int j = 0; j < B->getM(); j++) {
            assert (B->get(i, j) != 9);
            assert (B->get(i, j) == 7);
        }
    }
    //fill is working fine
    auto *C = *A + *B;
    auto *C2 = *A - *B;
    C->disp();

    A->disp();
    B->disp();
    auto *d = *A + *B;
    auto *d2 = *A - *B;
    d2->disp();
    C2->disp();
    assert(*C == *d); //sum is working fine
    assert(*C2 == *d2); //minus is working fine
    assert(C2->get(1,1) == 2);
    assert (C->equals(d));
    assert (!(*v == *u));
    delete v;
    delete u; //destructor is working fine

    //let us try the copy constructor

    auto *A_copied = new Matrix<int>(*A);
    assert (*A_copied == *A); //copy constructor working fine

    auto *a = new Matrix<int>(3, 1); //vector
    auto *b = new Matrix<int>(3); //vector, should be equal to v
    assert (*a == *b); //to create a vector of dim(v) = 3 is equal to a matrix (3,1)

    assert(a->getN() == 3);
    assert(a->getM() == 1);
    assert (a->get(0,0) == 0);
    assert (a->get(1,0) == 0);
    assert (a->get(2,0) == 0);

    a->disp(); //should display 0 0 0

    a->fill(5);
    a->disp(); //should display 5 5 5
    assert (a->get(0,0) == 5);
    assert (a->get(1,0) == 5);
    assert (a->get(2,0) == 5); //fill works fine

    //testing set and get
    a->set(1,0, 3);
    assert (a->get(1,0) == 3);
    //set and get are working fine

    assert(a->size()[0] == 3);
    assert(a->size()[1] == 1);
    assert(A->size()[0] == 3);
    assert(A->size()[1] == 3);
    //size is working fine

    //length is working fine
    assert(a->length() == 3);
    assert(A->length() == 3);

    assert(A->max() == 9);
    A->set(2,2, 14);
    assert(A->max() == 14);

    assert(A->min() == 9);
    A->set(1,2, 0);
    assert(A->min() == 0);


    auto *M = new Matrix<int>(2, 2);
    auto *N = new Matrix<int>(3, 3);
    M->fill(3);
    N->fill(2);
    assert (M->det() == 0);
    assert (N->det() == 0);
    auto *H = new Matrix<int>(3,3);
    int num = 1;
    for (int j = 0; j < 3; j++){
        for (int i = 0; i < 3; i++){
        H->set(i, j, num);
        num += 1;
        }
    }
    H->disp();
    assert (H->det() == 0);

    A->set(0,0, 7);
    A->set(0,1, 1);
    A->set(0,2, 3);

    A->set(1,0, 2);
    A->set(1,1, 4);
    A->set(1,2, 1);

    A->set(2,0, 1);
    A->set(2,1, 5);
    A->set(2,2, 1);

    A->disp();

    assert (A->det() == 10);
    //determinant is working fine
    A->fill(3);
    B->fill(3);
    auto *D = *A - *B;
    A->normalize();
    assert(A->get(0,0) == 1);
    assert(A->get(0,1) == 1);
    assert(A->get(0,2) == 1);
    assert(A->get(1,0) == 1);
    assert(A->get(2,0) == 1);
    D->set(1,1, 1);
    D->set(2,2, 2);
    D->set(0,0, 1);

    assert (D->is_diag());
    D->transpose();
    D->disp();
    assert (D->is_square());
    assert (A->is_square());
    assert (B->is_square());


    assert (!(a->is_square()));

    assert (a->is_vector());
    assert (!(a->is_symmetric()));

    assert((A->is_square()));

    auto *id = new Matrix<int>(3,3);
    id->set(0,0, 1);
    id->set(1,1,1);
    id->set(2,2,1);
    assert (id->is_identity());
    assert (id != A);

    auto *id_transp = id->transpose();
    auto *H_transp = H->transpose(); //transpose is working fine
    H->disp(); //transpose is working fine
    H_transp->disp();
    id_transp->disp();

    auto *D_Ht = *A - *H;
    A->disp();
    H->disp();
    D_Ht->disp();

    assert (D_Ht->get(0,2) == -6);
    assert (D_Ht->get(0,0) == 0);
    H->save_to_file("nombre.txt");
    A->save_to_file("intento.txt");
    //save to file is working well

    auto* H_copy = H->clone();
    assert (*H_copy == *H);

    auto* P = new Matrix<int> (2,2);
    P->set(0,0 ,8);
    P->set(1,0 ,19);
    P->set(0,1 ,-20);
    P->set(1,1 ,11);
    P->disp();

    auto* PF = new Matrix<float> (3, 3);
    PF->set(0,0 ,0);
    PF->set(1,0 ,0);
    PF->set(2,0 ,1);

    PF->set(0,1 ,2);
    PF->set(1,1 ,-1);
    PF->set(2,1 ,3);

    PF->set(0, 2, 1);
    PF->set(1, 2, 1);
    PF->set(2, 2, 4);

    cout << "PF\n";
    PF->disp();

    auto *F = PF->getAdjoint();
    F->disp(); //adjoint working fine
    F = PF->inverse(); //inverse is working!

    cout << "F\n";
    F->disp();
    //all boolean functions are working fine

    auto* PF2 = new Matrix<int> (3, 3);
    PF2->set(0,0 ,0);
    PF2->set(1,0 ,0);
    PF2->set(2,0 ,1);

    PF2->set(0,1 ,2);
    PF2->set(1,1 ,-1);
    PF2->set(2,1 ,3);

    PF2->set(0, 2, 1);
    PF2->set(1, 2, 1);
    PF2->set(2, 2, 4);

    auto* F2 = PF2->inverse();
    F2->disp(); //the inverse does not work so well when we do not use integers
    delete(F2);
    delete(PF2);
    delete(F);
    delete(PF);



    //let us test multiplication
    auto* K = new Matrix<float>(3,3);
    assert(K->get(0,0)==0);
    *K *= (float) 5.0; //mult by num
    assert(K->get(0,0)==0);
    //it shouldn't change

    K->set(0,0, 3);
    assert(K->get(0,0)==3);
    *K *= (float) 5.0; //mult by num
    assert(K->get(0,0)==15);
    assert(K->get(0,1)==0);

    K->set(0,1, 2);
    *K *= (float) 2.0; //mult by num
    assert(K->get(0,0)==30);
    assert(K->get(0,1)==4);
    //scalar mult is working fine.

    auto* K2 = new Matrix<float>(3,3);
    *K *=  *K2; // K should become 0
    assert (*K == *K2); //multiplication is working fine
    K->disp();
    assert(K->get(0,0)==0);
    assert(K->get(0,1)==0);
    //multiplication is working fine

    K2->set(0,1, 5.4);
    cout << K2->norm();
    assert(abs(K2->norm() - 5.4) < 0.000001); //norm is working fine


    return 0;

}