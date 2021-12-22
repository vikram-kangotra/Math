#pragma once 

#include <cstddef> 
#include <memory>
#include <iostream>
#include <vector>
#include <initializer_list>

namespace {
    template <typename T>
    class MatrixProxy;

    template <typename T>
    class BaseMatrix {
        public:
            BaseMatrix(const size_t Nrows, const size_t Ncols)
            : mRows{Nrows}
            , mCols{Ncols} {}

            virtual ~BaseMatrix() {}

            size_t size() const { return mRows * mCols; }

            virtual T& get(size_t index) = 0;
            virtual const T& get(size_t index) const = 0;

            const size_t rows() const { return mRows; }
            const size_t cols() const { return mCols; }

            T* begin() { return &get(0); }
            const T* begin() const { return &get(0); }
            T* end() { return begin() + size(); }
            const T* end() const { return begin() + size(); }

            auto operator[](size_t row) { return MatrixProxy<T>{this, row}; }
            const auto operator[](size_t row) const { return MatrixProxy<T>{this, row}; }

            void swapDimension() { std::swap(mRows, mCols); }

        private:
            size_t mRows;
            size_t mCols;
    };
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const BaseMatrix<T>& matrix) {
    for (auto i = 0; i < matrix.rows(); ++i) {
        for (auto j = 0; j < matrix.cols(); ++j) {
            out << matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

namespace {
    template <typename Matrix, typename T>
    class MatrixOperation;
}

template <typename T, size_t rows = 0, size_t cols = 0>
class Matrix : public BaseMatrix<T>, public MatrixOperation<Matrix<T, rows, cols>, T> {
    public:
        Matrix() 
        : BaseMatrix<T>{rows, cols}
        , MatrixOperation<Matrix<T, rows, cols>, T>{this} {}

        T& get(size_t index) override { return mData[index]; }
        const T& get(size_t index) const override { return mData[index]; }

    private:
        T mData[rows * cols];
};

template <typename T>
class Matrix<T, 0, 0> : public BaseMatrix<T>, public MatrixOperation<Matrix<T>, T> {
    public:
        Matrix(const size_t rows, const size_t cols)
        : BaseMatrix<T>{rows, cols} 
        , MatrixOperation<Matrix<T>, T>{this}
        , mData{std::make_unique<T[]>(rows * cols)} {}

        T& get(size_t index) { return mData[index]; }
        const T& get(size_t index) const { return mData[index]; }

    private:
        std::unique_ptr<T[]> mData;
};

namespace {

template <typename T>
class MatrixProxy {
    public:
        MatrixProxy(BaseMatrix<T>* matrix, size_t row)
        : mMatrix{matrix}
        , const_matrix{nullptr} 
        , mRow{row} {}

        MatrixProxy(const BaseMatrix<T>* const const_matrix, size_t row)
        : mMatrix{nullptr} 
        , const_matrix{const_matrix}
        , mRow{row} {}

        T& operator[](size_t cols) { 
            const auto index = mRow * mMatrix->cols() + cols;
            if (mRow >= mMatrix->rows() || cols >= mMatrix->cols()) {
                throw std::out_of_range("Failed attempt to access memory beyond assigned");
            }
            return mMatrix->get(index);
        }

        const T& operator[](size_t cols) const {
            const auto index = mRow * const_matrix->cols() + cols;
            if (mRow >= const_matrix->rows() || cols >= const_matrix->cols()) {
                throw std::out_of_range("Failed aatempt to access memory beyonf assigned");
            }
            return const_matrix->get(index);
        }

    private:
        BaseMatrix<T>* mMatrix;
        const BaseMatrix<T>* const const_matrix;
        size_t mRow;
};

namespace {
    template <typename GenMatrix, typename T>
    class MatrixOperation;

    template <typename Matrix, typename T>
        auto zip(Matrix& first, const Matrix& second) {
            std::vector<std::pair<T&, const T&>> container;
            if (MatrixOperation<Matrix, T>::haveEqualDimension(first, second)) {
                throw std::runtime_error("Cannot zip matrices because they have unequal size");
            }
            for (auto i = 0; i < first.rows(); ++i) {
                for (auto j = 0; j < second.cols(); ++j) {
                    container.push_back({first[i][j], second[i][j]});
                }
            }
            return container;
        }
    }

    template <typename GenMatrix, typename T>
    class MatrixOperation {
        public:
            MatrixOperation(GenMatrix* matrix)
            : mMatrix(*matrix) {}

            void identity() {
                if (!mMatrix.isSquareMatrix()) {
                    throw std::runtime_error("A Non-square matrix has no identity matrix representation");
                }

                for (auto i = 0; i < mMatrix.rows(); ++i) {
                    for (auto j = 0; j < mMatrix.cols(); ++j) {
                        mMatrix[i][j] = ( i == j );
                    }
                }
            }

            GenMatrix transpose() {
                GenMatrix matrix{mMatrix};
                matrix.swapDimension();
                for (auto i = 0; i < matrix.rows(); ++i) {
                    for (auto j = 0; j < matrix.cols(); ++j) {
                        if ( i == j )
                            continue;
                        matrix[i][j] = mMatrix[j][i];
                    }
                }
                return matrix;
            }

            GenMatrix& set(std::initializer_list<T> list) {
                const auto limit = mMatrix.size();

                for (size_t index = 0; auto l : list) {
                    if (index >= limit) {
                        throw std::out_of_range("Extra elements were provided as parameter");
                    }
                    mMatrix.get(index) = l;
                    ++index;
                }

                return mMatrix;
            }

            static bool haveEqualDimension(const GenMatrix& first, const GenMatrix& second) {
                return first.rows() == second.rows() && first.cols() == second.cols();
            }

            bool isSquareMatrix() const {
                return mMatrix.rows() == mMatrix.cols();
            }

            GenMatrix& operator-() {
                GenMatrix matrix{mMatrix};
                for (auto& elem : matrix) {
                    elem = -elem;
                }
                return matrix;
            }

            auto& operator*=(T scalar) {
                for (auto& elem : mMatrix) {
                   elem *= scalar; 
                }
                return mMatrix;
            }

            auto& operator+=(const GenMatrix& other) {
                if (!haveEqualDimension(mMatrix, other)) {
                    throw std::runtime_error("Matrices should have same dimension for the operation to proceed successfully");
                }
                for (auto& [a, b] : zip<GenMatrix, T>(mMatrix, other)) { a += b;}
                return mMatrix;
            }

            auto& operator-=(const GenMatrix& other) {
                if (!haveEqualDimension(mMatrix, other)) {
                    throw std::runtime_error("Matrices should have same dimension for the operation to proceed successfully");
                }
                for (auto& [a, b] : zip<GenMatrix, T>(mMatrix, other)) { a -= b; }
                return mMatrix;
            }

            auto& operator/=(T scalar) {
                if (scalar == 0) {
                    throw std::runtime_error("Cannot divide matrix by 0.");
                }
                for (auto& elem : mMatrix) { elem /= scalar; }
                return mMatrix;
            }

            auto& operator*(const GenMatrix& other) const {
                if (mMatrix.cols() != other.rows()) {
                    throw std::runtime_error("Can't multiplt matrices. GenMatrix's column does not the match other matrix's rows.");
                }
                GenMatrix matrix{mMatrix};

                for (auto i = 0; i < matrix.rows(); ++i) {
                    for (auto j = 0; j < matrix.cols(); ++j) {
                        for (auto k = 0; k < mMatrix.cols(); ++k) {
                            matrix[i][j] += mMatrix[i][k] + other[k][j];
                        }
                    }
                }
                return matrix;
            }

            auto operator*(T scalar) const { return GenMatrix{mMatrix} *= scalar; }
            auto operator+(const GenMatrix& other) const { return GenMatrix{mMatrix} += other; }
            auto operator-(const GenMatrix& other) const { return GenMatrix{mMatrix} -= other; }
            auto operator/(T scalar) { return GenMatrix{mMatrix} /= scalar; }

            auto operator*=(const GenMatrix& other) { return mMatrix *= other; }
            auto operator/(const GenMatrix& other) const { return mMatrix * other.inverse(); }
            auto operator/=(const GenMatrix& other) { return mMatrix *= other.inverse(); }

            T determinant() const {
                if (!mMatrix.isSquareMatrix()) {
                    throw std::runtime_error("Determinnant cannot be calculated for a Non-Sqaure matrix.");
                }
                if (mMatrix.rows() == 2 && mMatrix.cols() == 2) {
                    return mMatrix[0][0] * mMatrix[1][1] - mMatrix[0][1] * mMatrix[1][0];
                }
                T answer{};

                for (auto i = 0; i < mMatrix.rows(); ++i) {
                    answer += getSign(0,i) * mMatrix[0][i] * getMatrixPart(0,i).determinant();
                }
                return answer;
            }

            auto inverse() const { return adjoint() / determinant(); }

        private:
            auto cofactor() const {
                if (!mMatrix.isSquareMatrix()) {
                    throw std::runtime_error("GenMatrix needs to be a square matrix.");
                }
                GenMatrix matrix{mMatrix};

                for (auto i = 0; i < matrix.rows(); ++i) {
                    for (auto j = 0; j < matrix.cols(); ++j) {
                        matrix[i][j] = getSign(i,j) * getMatrixPart(i,j).determinant();
                    }
                }
                return matrix;
            }

            auto adjoint() const { return cofactor().transpose(); }

            auto getMatrixPart(size_t rows, size_t cols) const {
                Matrix<T> matrix{mMatrix.rows() - 1, mMatrix.cols() - 1};

                for (auto i = 0; i < mMatrix.rows(); ++i) {
                    for (auto j = 0; j < mMatrix.cols(); ++j) {
                        if (i == rows || j == cols)
                            continue;
                        matrix[i - (i>rows)][j - (j>cols)] = mMatrix[i][j];
                    }
                }
                return matrix;
            }

            auto getSign(size_t rows, size_t cols) const { 
                return ((rows+cols)%2) ? -1 : 1; 
            }

        private:
            GenMatrix& mMatrix;
    };
}
