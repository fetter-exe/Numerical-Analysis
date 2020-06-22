#include "matrix.h"
#include <sstream>   // stringstream
#include <stdlib.h>  // srand(); frand();
#include <time.h>    // time;
#include <memory>    // unique_ptr
#include <cmath>     // pow()

    Matrix::Matrix(){}

    Matrix::Matrix(const Matrix& other){
        this->rows = other.rows;
        this->cols = other.cols;
        data.resize(rows);
        for(unsigned int i=0; i<rows ; ++i){ // Fills the Matrix
            copy(other.data[i].begin(), other.data[i].end(), back_inserter(this->data[i]));
        }
    }

    Matrix::Matrix(std::vector<std::vector<double>>& inputData){
        rows = inputData.size();
        cols = inputData[0].size();
        data.resize(rows);
        for(unsigned int i=0; i<rows ; ++i){ // Fills the Matrix
            copy(inputData[i].begin(), inputData[i].end(), back_inserter(data[i]));
        }
    }

    Matrix::Matrix(int numberOfRows,int numberOfCols){
        if(numberOfRows < 0 || numberOfCols < 0){
            throw std::invalid_argument("negative arguments are not allowed");
        }
        rows = numberOfRows;
        cols = numberOfCols;
        data.resize(rows);
        srand(time(NULL)); // Guarantees a random seed
        for(unsigned int i=0; i<rows; ++i){ // Fills the Matrix
            for(unsigned int j=0; j<cols; ++j){
                rand()%2 == 0  ?  data[i].push_back((double)rand())  :  data[i].push_back(-(double)rand());
            }
        }
    }

    Matrix::Matrix(int numberOfRows,int numberOfCols, double minLimit, double maxLimit){
        if(numberOfRows < 0 || numberOfCols < 0){
            throw std::invalid_argument("negative arguments are not allowed");
        }
        rows = numberOfRows;
        cols = numberOfCols;

        data.resize(rows);
        for(unsigned int i = 0; i < rows; ++i){
            data[i].resize(cols);
        }

        srand(time(NULL)); // Guarantees a random seed
        for(unsigned int i = 0; i < rows; ++i){ // Fills the Matrix
            for(unsigned int j = 0; j < cols; ++j){
                data[i][j] = minLimit + rand()%((int)(maxLimit - minLimit + 1));
            }
        }
    }

    Matrix::Matrix(int numberOfRows, int numberOfCols, std::string elements){
        std::stringstream ss(elements);
        rows = numberOfRows;
        cols = numberOfCols;

        data.resize(rows);
        for(unsigned int i = 0; i < rows; ++i){
            data[i].resize(cols);
        }

        for(unsigned int i = 0; i < rows; ++i){ // loop to fill in the data
            for(unsigned int j = 0; j < cols; ++j){
                if(!ss.eof()){
                    ss >> data[i][j];
                }else{ // if eof has been reached, fill the remaining data with 0s
                    data[i][j] = 0;
                }
            }
        }
    }

    unsigned int Matrix::Rows() const{
        return this->rows;
    }

    unsigned int Matrix::Cols() const{
        return this->cols;
    }

    double Matrix::operator()(const unsigned int& I, const unsigned int& J) const{
        if(I > this->rows || J > this->rows){
            throw std::invalid_argument("tried to access an element that does not exist");
        }
        return this->data[I][J];
    }

    Matrix Matrix::operator+(const Matrix& other) const{
        if(this->rows != other.rows || this->cols != other.cols){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        std::unique_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
        for(unsigned int i=0; i<returnedMatrix->rows; ++i){ // Fills the Matrix
           for(unsigned int j=0; j<returnedMatrix->cols; ++j){
               returnedMatrix->data[i][j] = this->data[i][j] + other.data[i][j];
           } 
        }
        return *returnedMatrix;
    }

    Matrix Matrix::operator-(const Matrix& other) const{
        if(this->rows != other.rows || this->cols != other.cols){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        std::unique_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
        for(unsigned int i=0; i<returnedMatrix->rows; ++i){ // Fills the Matrix
           for(unsigned int j=0; j<returnedMatrix->cols; ++j){
               returnedMatrix->data[i][j] = this->data[i][j] - other.data[i][j];
           } 
        }
        return *returnedMatrix;
    }

    Matrix Matrix::operator*(const Matrix& other) const{
        if(this->cols != other.rows){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        std::vector<std::vector<double>> returnedData;
        returnedData.resize(this->rows);
        for(unsigned int i=0; i<returnedData.size(); ++i){ // Fills the Matrix
            returnedData[i].resize(other.cols);
            for(unsigned int j=0; j<returnedData[0].size(); ++j){
                returnedData[i][j] = 0.0;
                for(unsigned int X=0; X<this->cols; ++X){
                    returnedData[i][j] += this->data[i][X]*other.data[X][j];
                }
            }
        }
        std::unique_ptr<Matrix> returnedMatrix (new Matrix(returnedData));
        return *returnedMatrix;
    }

    Matrix Matrix::operator*(const double& scalar) const{
        std::unique_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                returnedMatrix->data[i][j] = this->data[i][j] * scalar;
            }
        }
        return *returnedMatrix;
    }

    Matrix operator*(const double& scalar, const Matrix& matrix){
        std::unique_ptr<Matrix> returnedMatrix(new Matrix(matrix.rows, matrix.cols));
        for(unsigned int i=0; i<returnedMatrix->rows; ++i){
            for(unsigned int j=0; j<returnedMatrix->cols; ++j){
                returnedMatrix->data[i][j] = matrix.data[i][j] * scalar;
            }
        }
        return *returnedMatrix;
    }

    void Matrix::ShowContent(){
        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                std::cout<<"data["<<i<<"]["<<j<<"]: "<<data[i][j]<<std::endl;
            }
        }
    }

    void Matrix::ShowMatrix(){
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    bool Matrix::IsNull() const{
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                if(data[i][j] != 0){
                    return false;
                }
            }
        }
        return true;
    }

    bool Matrix::IsDiagonal() const{
        if(IsNull()){
            return false;
        }
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                if(i != j  &&  data[i][j] != 0){
                    return false;
                }
            }
        }
        return true;
    }

    bool Matrix::IsTriangular() const{
        bool isTriang, isUpper = true, isLower = true;
        if(rows != cols  ||  IsNull()){
            return false;
        }
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                if(i < j  &&  data[i][j] != 0){
                    isUpper = false;
                }
                if(i > j  &&  data[i][j] != 0){
                    isLower = false;
                }
            }
        }
        isTriang = (isUpper != isLower) ? true : false;
        return isTriang;
    }

    bool Matrix::IsSymmetric() const{
        if(rows != cols){
            return false;
        }
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                if(data[i][j] != data[j][i]){
                    return false;
                }
            }
        }
        return true;
    }

    double Matrix::Determinant() const{
        // Ensures only square matrixes can use this method
        if(rows != cols){
            throw std::length_error("non-square matrixes can't have a determinant");
        }

        double det = 1, total = 1;
        unsigned int idx, n = rows;
        std::vector<double> temp(rows);
        std::vector<std::vector<double>> mat = data;

        if(IsNull()){ // Optimizes for Null/Zero matrixes
            det = 0;
            return det;
        }
        if(IsDiagonal() || IsTriangular()){ // Optimizes for diagonal/triangular matrixes
            for(unsigned int k = 0; k < n; ++k){
                det *= mat[k][k];
            }
            return det;
        }

        if(n == 1){ // Optimizes for small matrixes
            det = data[0][0];
            return det;
        }else if(n == 2){
            det = data[0][0]*data[1][1] - data[0][1]*data[1][0];
            return det;
        }
        else if(n == 3){
            det = 0;
            det += data[0][0] * (data[1][1]*data[2][2] - data[1][2]*data[2][1]);
            det -= data[0][1] * (data[1][0]*data[2][2] - data[1][2]*data[2][0]);
            det += data[0][2] * (data[1][0]*data[2][1] - data[1][1]*data[2][0]);
            return det;
        }

        // Loop for traversing diagonal elements
        for(unsigned int k = 0; k < n; ++k){
            idx = k;

            // Finding the index that has a non-zero value
            while(mat[idx][k] == 0  &&  idx < n){
                ++idx;
            }

            // If there is non-zero element
            if(idx == n){
                continue; // The determinant of matrix as zero
            }

            if(idx != k){
                // Loop for swaping the index row and diagonal row
                for(unsigned int j = 0; j < n; ++j){
                    std::swap(mat[idx][j], mat[k][j]);
                }
                // Determinant changes sign when shifting rows
                det *= pow(-1, idx-k);
            }

            // Storing the values of the diagonal row elements
            for(unsigned int j = 0; j < n; ++j){
                temp[j] = mat[k][j];
            }

            // Traversing every row below the diagonal element
            for(unsigned int i = k+1; i < n; ++i){
                double diag_el = temp[k]; // Value of diagonal element
                double nxtr_el = mat[i][k]; // Value of next row element

                // Traversing every column of row and multiplying to every row
                for(unsigned int j = 0; j < n; ++j){
                    // Multiplying to make diagonal and next row elements equal
                    mat[i][j] = (diag_el * mat[i][j]) - (nxtr_el * temp[j]);
                }

                total *= diag_el; // det(kA) = k * det(A)
            }
        }

        // Multiplying the diagonal elements to get determinant
        for(unsigned int k = 0; k < n; ++k){
            det *= mat[k][k];
        }

        return (det/total); // det(A) = det(kA) / k
    }

    Matrix Matrix::Submatrix(int i, int j) const{
        if(rows <= 1  ||  cols <= 1){
            return Matrix(0,0);
        }

        std::vector<std::vector<double>> smData;
        unsigned int smi, smj;

        if(i < 0){ i = rows; }
        if(j < 0){ j = cols; }

        if((unsigned)i >= rows  &&  (unsigned)j >= cols){
            smData = std::vector<std::vector<double>>(rows,std::vector<double>(cols));
        }else if((unsigned)i >= rows  &&  (unsigned)j < cols){
            smData = std::vector<std::vector<double>>(rows,std::vector<double>(cols-1));
        }else if((unsigned)j >= cols  &&  (unsigned)i < rows){
            smData = std::vector<std::vector<double>>(rows-1,std::vector<double>(cols));
        }else{
            smData = std::vector<std::vector<double>>(rows-1,std::vector<double>(cols-1));
        }

        smi = 0;
        for(unsigned int i2 = 0; i2 < rows; ++i2){
            if(i2 != (unsigned)i){
                smj = 0;
                for(unsigned int j2 = 0; j2 < cols; ++j2){
                    if(j2 != (unsigned)j){
                        smData[smi][smj] = data[i2][j2];
                        ++smj;
                    }
                }
                ++smi;
            }
        }

        return Matrix(smData);
    }

    Matrix Matrix::Transpose() const{
        std::vector<std::vector<double>> tData;

        tData.resize(cols);
        for(unsigned int i = 0; i < cols; ++i){
            tData[i].resize(rows);
        }

        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                tData[j][i] = data[i][j];
            }
        }

        return Matrix(tData);
    }

    Matrix Matrix::Inverse() const{
        // Ensures only valid matrixes can use this method
        if(rows != cols){
            throw std::length_error("Non-square matrixes don't have an inverse");
        }
        if(Determinant() == 0){
            throw std::invalid_argument("Singular matrixes (det = 0) don't have an inverse");
        }

        // The Inverse of a Matrix can be calculated using the formula:
        // A^-1 = (1 / det(A)) * C^t;  C(nxn) being the cofactor matrix
        // of A(nxn), which is defined by the followig rules:
        // C_{i,j} = (-1)^(i+j) * det(M_{i,j});  M_{i,j} being a submatrix of A
        // of size (n-1)x(n-1) generated by deleting the ith row and jth column of A

        std::vector<std::vector<double>> cData, mData;

        cData.resize(rows); // Resizes the data for the C matrix
        for(unsigned int i = 0; i < rows; ++i){
            cData[i].resize(cols);
        }
        mData.resize(rows-1); // Resizes the data for the M matrixes
        for(unsigned int i = 0; i < (rows-1); ++i){
            mData[i].resize(cols-1);
        }

        if(rows == 1){
            cData[0][0] = 1/data[0][0];
            return Matrix(cData);
        }

        // Loop to construct the C_{i,j} element
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                unsigned int im, jm;
                // Loop to construct the M_{i,j} submatrix
                im = 0;
                for(unsigned int i2 = 0; i2 < rows; ++i2){
                    if(i2 != i){
                        jm = 0;
                        for(unsigned int j2 = 0; j2 < cols; ++j2){
                            if(j2 != j){
                                mData[im][jm] = data[i2][j2];
                                ++jm;
                            }
                        }
                        ++im;
                    }
                }
                Matrix m = Matrix(mData);
                cData[i][j] = pow(-1,i+j) * m.Determinant();
                if(cData[i][j] == -0){ cData[i][j] = 0; }
            }
        }

        Matrix inv = Matrix((1/Determinant()) * Matrix(cData).Transpose());

        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                if(inv.data[i][j] == -0){
                    inv.data[i][j] = 0;
                }
            }
        }

        return inv;
    }

    Matrix::~Matrix(){}


    // IdentityMatrix implementation
    IdentityMatrix::IdentityMatrix(int Order): Matrix(){
        if(Order < 0){
            throw std::invalid_argument("negative arguments are not allowed");
        }
        this->rows = Order;
        this->cols = Order;
        data.resize(rows);
        for(unsigned int i=0; i<rows; ++i){
            data[i].resize(cols);
        }
        for(unsigned int i=0; i<rows; ++i){
            data[i][i] =  1.0;
        }
    }

    IdentityMatrix::IdentityMatrix(const IdentityMatrix& other): Matrix(){
        this->rows = other.rows;
        this->cols = other.rows;
        data.resize(rows);
        for(unsigned int i=0; i<rows; ++i){
            data[i].resize(cols);
        }
        for(unsigned int i=0; i<rows; ++i){
            data[i][i] =  1.0;
        }
    }

    Matrix IdentityMatrix::operator*(const Matrix& other) const{
        if(this->cols != other.Rows()){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        return other;
    }

    double IdentityMatrix::Determinant() const{
        return 1.0;
    }

    IdentityMatrix::~IdentityMatrix(){}


    // VandermondeMatrix implementation
    VandermondeMatrix::VandermondeMatrix(std::vector<double>& inputData, int numberOfCols): Matrix(){
        rows = inputData.size();
        cols = numberOfCols;
        data.resize(rows);
        for(unsigned int i=0; i<rows; ++i){
            data[i].resize(cols);
            for(int j=0; j<numberOfCols; ++j){
                this->data[i][j] = pow(inputData[i], j);
            }
        }
    }

    VandermondeMatrix::VandermondeMatrix(std::vector<double>& inputData): Matrix(){
        rows = inputData.size();
        cols = inputData.size();
        data.resize(rows);
        for(unsigned int i=0; i<rows; ++i){
            data[i].resize(cols);
            for(unsigned int j=0; j<cols; ++j){
                this->data[i][j] = pow(inputData[i], j);
            }
        }
    }

    double VandermondeMatrix::Determinant() const{
        // Ensures only square matrixes can use this method
        if(rows != cols){
            throw std::length_error("non-square matrixes can't have a determinant");
        }

        double Product = 1.0;
        for(int i = rows-1; i > 0; --i){
            for(int j = i-1; j >= 0 ; --j){
                Product *= (data[i][1] - data[j][1]);
            }
        }
        return Product;
    }

    VandermondeMatrix::~VandermondeMatrix(){}
