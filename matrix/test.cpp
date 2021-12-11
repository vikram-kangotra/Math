#include <iostream>
#include "matrix.h"

int main(){

    CTMatrix<double, 3, 3> matrix;

    matrix.set(10, 2, 10, 4, 3, 6, 5, 4, 18);
    
    std::cout << matrix.determinant();

    return 0;
}
