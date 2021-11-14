#pragma once

#include <vector>

template <typename T>
class Matrix;

namespace {
template <typename T>
    class MatrixProxy {
        public:
            MatrixProxy(Matrix<T>* matrix, size_t row);
            MatrixProxy(const Matrix<T>* const const_matrix, size_t row);

            T& operator[](size_t cols);
            const T& operator[](size_t cols) const;
        protected:
        private:
            Matrix<T>* mMatrix;
            const Matrix<T>* const_matrix;
            size_t mRow;
    };
}

template <typename T>
class Matrix {
    public:
        Matrix(size_t rows, size_t cols)
        : mRows(rows), mCols(cols) {}

        virtual ~Matrix() {}

        virtual T& get(uint32_t index) = 0;
        virtual const T& get(uint32_t index) const = 0;

        const auto rows() const {
            return mRows;
        }

        const auto cols() const {
            return mCols;
        }

        template <typename... R>
        void set(R&&... elems) {
            size_t i = 0;
            uint32_t limit = mRows * mCols;

            static auto fillData = [&](T&& elem) {
                if (i >= limit) {
                    throw std::out_of_range("Extra elements were provied as parameter to set function");
                }
                this->get(i++) = elem;
            };

            (fillData(std::forward<T>(elems)), ...);
        }

        auto operator[](size_t row) {
            return MatrixProxy<T>(this, row);
        }

        const auto operator[](size_t row) const {
            return MatrixProxy<T>(this, row);
        }

        static auto haveEqualDimension(const Matrix<T>& first, const Matrix<T>& second) {
            return first.rows() == second.rows() && first.cols() == second.cols();
        }

        auto isSquareMatrix() const {
            return mRows == mCols;
        }

        auto begin() {
            return &get(0);
        }

        const auto begin() const {
            return &get(0);
        }

        auto end() {
            return &get(mRows * mCols - 1);
        }

        const auto end() const {
            return &get(mRows * mCols - 1);
        }

        void swapDimension() {
            std::swap(mRows, mCols);
        }

    protected:
    private:
        size_t mRows;
        size_t mCols;
};

template <typename T>
MatrixProxy<T>::MatrixProxy(Matrix<T>* matrix, size_t row)
: mMatrix(matrix), mRow(row) {}

template <typename T>
MatrixProxy<T>::MatrixProxy(const Matrix<T>* const const_matrix, size_t row)
: const_matrix(const_matrix), mRow(row) {}

template <typename T>
T& MatrixProxy<T>::operator[](size_t cols) {
    const auto index = mRow * mMatrix->cols() + cols;
    if (mRow >= mMatrix->rows() || cols >= mMatrix->cols()) {
        throw std::out_of_range("Failed attempt to access memory beyond assigned");
    }
    return mMatrix->get(index);
}

template <typename T>
const T& MatrixProxy<T>::operator[](size_t cols) const {
    const auto index = mRow * const_matrix->cols() + cols;
    if (mRow >= const_matrix->rows() || cols >= const_matrix->cols()) {
        throw std::out_of_range("Failed attempt to access memory beyond assigned");
    }
    return const_matrix->get(index);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix) {
    for(auto i=0; i<matrix.rows(); ++i) {
        for(auto j=0; j<matrix.cols(); ++j) {
            out << matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template <typename GenMatrix, typename T>
class MatrixOperation;

// CompileTimeMatrix (CTMatrix)

template <typename T, size_t rows, size_t cols>
class CTMatrix : public Matrix<T>, public MatrixOperation<CTMatrix<T, rows, cols>, T> {
    public:
        CTMatrix()
        : Matrix<T>(rows, cols),
          MatrixOperation<CTMatrix<T, rows, cols>, T>(this),
          mData{} {}

        virtual T& get(uint32_t index) {
            return mData[index];
        }

        virtual const T& get(uint32_t index) const {
            return mData[index];
        }
    protected:
    private:
        T mData[rows * cols];
};

// RunTimeMatrix (RTMatrix)

template <typename T>
class RTMatrix : public Matrix<T> , public MatrixOperation<RTMatrix<T>, T>{
    public:
        RTMatrix(size_t rows, size_t cols)
        : Matrix<T>(rows, cols),
          MatrixOperation<RTMatrix<T>, T>(this) {
            mData = new T[rows * cols]{};
        }

        virtual ~RTMatrix() {
            delete[] mData;
        }

        virtual T& get(uint32_t index) {
            return mData[index];
        }

        virtual const T& get(uint32_t index) const {
            return mData[index];
        }

    protected:
    private:
        T* mData;
};

namespace {

    template <typename GenMatrix, typename T>
    auto zip(GenMatrix& first, const GenMatrix& second) {
        std::vector<std::pair<T&, const T&>> container;
        if (!Matrix<T>::haveEqualDimension(first, second)) {
            throw std::runtime_error("Matrices should have same dimension for the operation to proceed successfully");
        }
        for (int i=0; i<first.rows(); ++i) {
            for(int j=0; j<first.cols(); ++j) {
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
            if ( !mMatrix.isSquareMatrix() ) {
                throw std::runtime_error("A Non-square matrix has no identity matrix equivalent. Can't convert matrix to identity matrix.");
            }
            
            for (int i=0 ; i < mMatrix.rows(); ++i) {
                for (int j=0 ; j < mMatrix.cols(); ++j) {
                    mMatrix[i][j] = ( i == j );
                }
            }
        }

        GenMatrix transpose() {
            GenMatrix matrix = mMatrix;
            matrix.swapDimension();
            for (int i=0 ; i < matrix.rows(); ++i){
                for (int j=0 ; j < matrix.cols(); ++j){
                    if ( i == j ) 
                        continue;
                    matrix[i][j] = mMatrix[j][i];
                }
            }
            return matrix;
        }

        T determinant() const {
            if ( !mMatrix.isSquareMatrix() ) {
                throw std::runtime_error("Determinant cannot be calculated for a Non-square matrix.");
            }

            if ( mMatrix.rows() == 2 && mMatrix.cols() == 2 ) {
                return mMatrix[0][0] * mMatrix[1][1] - mMatrix[0][1] * mMatrix[1][0];
            }

            T answer{};

            for (int i=0 ; i < mMatrix.rows(); ++i){
                answer += getSign(0,i) * (mMatrix[0][i] * getMatrixPart(0,i).determinant());
            }
            return answer;
        }

        auto inverse() const {
            return adjoint() / determinant();
        }

        auto& operator-() {
            for (auto& elem : mMatrix) {
                elem = -elem;
            }
            return mMatrix;
        }

        auto& operator*=(T scalar) {
            for (auto& elem : mMatrix) {
                elem *= scalar;
            }
            return mMatrix;
        }

        auto operator*(T scalar) const {
            return GenMatrix{mMatrix} *= scalar;
        }

        auto operator+(const GenMatrix& other) const {
            return GenMatrix{mMatrix} += other;
        }

        auto& operator+=(const GenMatrix& other) {
            if (!Matrix<T>::haveEqualDimension(mMatrix, other)) {
                throw std::runtime_error("Matrices should have same dimension for the operation to proceed successfully");
            }

            for (auto& [a, b] : zip<GenMatrix, T>(mMatrix, other)) { a += b; }
            return mMatrix;
        }

        auto operator-(const GenMatrix& other) const {
            return GenMatrix{mMatrix} -= other;
        }

        auto& operator-=(const GenMatrix& other) {
            if (!Matrix<T>::haveEqualDimension(mMatrix, other)) {
               throw std::runtime_error("Matrices should have same dimension for the operation to proceed successfully");
            }

            for (auto& [a, b] : zip<GenMatrix, T>(mMatrix, other)) { a -= b; }
            return mMatrix;
        }

        auto& operator/=(T scalar) {
            if ( scalar == 0 ) {
                throw std::runtime_error("Cannot divide matrix by 0.");
            }
            for (auto& elem : mMatrix) {
                elem /= scalar;
            }
            return mMatrix;
        }
        
        auto operator/(T scalar) {
            return GenMatrix{mMatrix} /= scalar;
        }

        auto& operator*(const GenMatrix& other) const {
            if (mMatrix.cols() != other.rows() ) {
                throw std::runtime_error("Can't multiply matrices. Matrix's column does not match other matrix's row");
            }

            GenMatrix matrix = mMatrix;

            for (int i=0 ; i < matrix.rows(); ++i){
                for (int j=0 ; j < matrix.cols(); ++j){
                    for (int k=0 ; k < mMatrix.cols(); ++k){
                        matrix[i][j] += mMatrix[i][k] * other[k][j];
                    }
                }
            }

            return matrix;
        }

        auto operator*=(const GenMatrix& other) {
            return mMatrix *= other;
        }

        auto operator/(const GenMatrix& other) const {
            return mMatrix * other.inverse();
        }

        auto& operator/=(const GenMatrix& other) {
            return mMatrix *= other.inverse();
        }

    private:
        auto cofactor() const {
            if ( !mMatrix.isSquareMatrix ) {
                throw std::runtime_error("Matrix needs to be a square matrix.");
            }

            GenMatrix matrix{mMatrix};

            for (auto i=0 ; i < matrix.rows(); ++i){
                for (auto j=0 ; j < matrix.cols() ; ++j){
                    matrix[i][j] = getSign(i,j) * getMatrixPart(i,j).determinant();
                }
            }
            return matrix;
        }

        auto adjoint() const {
            return cofactor().transpose();
        }

        auto getMatrixPart(size_t rows, size_t cols) const {
            RTMatrix<T> matrix{mMatrix.rows() - 1, mMatrix.cols() - 1};

            for (auto i=0 ; i < mMatrix.rows(); ++i){
                for (auto j=0 ; j < mMatrix.cols(); ++j){ if ( i == rows || j == cols) 
                        continue;
                    matrix[i - (i>rows)][j - (j>cols)] = mMatrix[i][j];
                }
            }

            return matrix;
        }

        auto getSign(size_t rows, size_t cols) const {
            return ((rows+cols)%2)?-1:1;
        }

    protected:
    private:
        GenMatrix& mMatrix;
};

