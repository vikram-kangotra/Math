#include "matrix.hpp"

int main() {
    Matrix<float, 4, 4> matrix;

    matrix.set({10, 2, 3, 4, 5, 60, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});

    std::cout << matrix.transpose().determinant();
}
