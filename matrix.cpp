#include "matrix.h"
#include <sstream>   // stringstream
#include <stdlib.h>  // srand(); frand();
#include <time.h>    // time;
#include <memory>    // unique_ptr
#include <cmath>     // pow()
#include <iomanip>   // setprecision()

    // Matrix implementation
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

    Matrix::Matrix(int numberOfRows, int numberOfCols){
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

    Matrix::Matrix(int numberOfRows, int numberOfCols, double minLimit, double maxLimit){
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

    Matrix::Matrix(int numberOfRows, int numberOfCols, double value){
        this->rows = numberOfRows;
        this->cols = numberOfCols;
        this->data = std::vector<std::vector<double>>(rows,std::vector<double>(cols,value));
    }

    unsigned int Matrix::GetRows() const{
        return this->rows;
    }

    unsigned int Matrix::GetCols() const{
        return this->cols;
    }
    
    bool Matrix::IsEqual(const Matrix& other) const{
        if(this->rows != other.rows || this->cols != other.cols){
            return false;
        }
        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                if(data[i][j] != other.data[i][j]){
                    return false;
                }
            }
        }
        return true;
    }

    bool Matrix::operator==(const Matrix& other) const{
        return IsEqual(other);
    }

    bool Matrix::IsDifferent(const Matrix& other) const{
        if(this->rows != other.rows || this->cols != other.cols){
            return true;
        }
        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                if(data[i][j] != other.data[i][j]){
                    return true;
                }
            }
        }
        return false;
    }

    bool Matrix::operator!=(const Matrix& other) const{
        return IsDifferent(other);
    }

    double Matrix::operator()(const unsigned int& I, const unsigned int& J) const{
        if(I > this->rows || J > this->rows){
            throw std::invalid_argument("tried to access an element that does not exist");
        }
        return this->data[I][J];
    }

    Matrix Matrix::Add(const Matrix& other) const{
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

    Matrix Matrix::operator+(const Matrix& other) const{
        return Add(other);
    }

    Matrix Matrix::Subtract(const Matrix& other) const{
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

    Matrix Matrix::operator-(const Matrix& other) const{
        return Subtract(other);
    }

    Matrix Matrix::Multiply(const Matrix& other) const{
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

    Matrix Matrix::operator*(const Matrix& other) const{
        return Multiply(other);
    }

    Matrix Matrix::Multiply(const double& scalar) const{
        std::unique_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                returnedMatrix->data[i][j] = this->data[i][j] * scalar;
            }
        }
        return *returnedMatrix;
    }

    Matrix Matrix::operator*(const double& scalar) const{
        return Multiply(scalar);
    }

    Matrix operator*(const double& scalar, const Matrix& matrix){
        return matrix*scalar;
    }

    Matrix Matrix::Divide(const Matrix& other) const{
        return (this->Multiply(other.Inverse()));
    }

    Matrix Matrix::operator/(const Matrix& other) const{
        return Divide(other);
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
        if(GetRows() != GetCols()){ //Ensures only square matrixes can use this method
            throw std::length_error("non-square matrixes can't have a determinant");
        }

        double det = 1;
        unsigned int n = GetRows();
        std::vector<double> pivot(n,0);
        std::vector<std::vector<double>> aData = this->data;

        if(IsNull()){ //Optimizes for Null/Zero matrixes
            return 0;
        }
        if(IsDiagonal() || IsTriangular()){ //Optimizes for diagonal/triangular matrixes
            for(unsigned int i = 0; i < n; ++i){
                det *= aData[i][i];
            }
            return det;
        }

        if(n == 1){ //Optimizes for small matrixes
            return aData[0][0];
        }
        else if(n == 2){
            det = aData[0][0]*aData[1][1] - aData[0][1]*aData[1][0];
            return det;
        }
        else if(n == 3){
            det = aData[0][0] * (aData[1][1]*aData[2][2] - aData[1][2]*aData[2][1]);
            det -= aData[0][1] * (aData[1][0]*aData[2][2] - aData[1][2]*aData[2][0]);
            det += aData[0][2] * (aData[1][0]*aData[2][1] - aData[1][1]*aData[2][0]);
            return det;
        }

        // This algorithm was adapted from the book "Algoritmos Numéricos",
        // by professor Frederico Ferreira Campos filho, Ph.D.
        // All rights reserved. You can find more info about him and his work at
        // http://www2.dcc.ufmg.br/livros/algoritmosnumericos/

        for(unsigned int i = 0; i < n; ++i){
            pivot[i] = i;
        }
        for(unsigned int j = 0; j < (n-1); ++j){
            unsigned int p = j;
            double aMax = std::abs(aData[j][j]);
            for(unsigned int k = (j+1); k < n; ++k){
                if(std::abs(aData[k][j]) > aMax){
                    aMax = std::abs(aData[k][j]);
                    p = k;
                }
            }
            if(p != j){
                for(unsigned int k = 0; k < n; ++k){
                    std::swap(aData[j][k],aData[p][k]);
                }
                std::swap(pivot[j],pivot[p]);
                det *= -1;
            }
            det *= aData[j][j];
            if(std::abs(aData[j][j]) != 0){
                double r = 1/aData[j][j];
                for(unsigned int i = (j+1); i < n; ++i){
                    double mult = r * aData[i][j];
                    aData[i][j] = mult;
                    for(unsigned int k = (j+1); k < n; ++k){
                        aData[i][k] -= mult * aData[j][k];
                    }
                }
            }
        }
        det *= aData[n-1][n-1];

        return det;
    }

    Matrix Matrix::SubMatrix(int i, int j) const{
        if((rows==1 || cols==1) && (i==0 || j==0)){
            return Matrix(0,0);
        }

        std::vector<std::vector<double>> smData;
        unsigned int smi, smj;

        if(i < 0){ i = rows; }
        if(j < 0){ j = cols; }

        if((unsigned)i >= rows  &&  (unsigned)j >= cols){
            smData = std::vector<std::vector<double>>(rows,std::vector<double>(cols,0));
        }else if((unsigned)i >= rows  &&  (unsigned)j < cols){
            smData = std::vector<std::vector<double>>(rows,std::vector<double>(cols-1,0));
        }else if((unsigned)j >= cols  &&  (unsigned)i < rows){
            smData = std::vector<std::vector<double>>(rows-1,std::vector<double>(cols,0));
        }else{
            smData = std::vector<std::vector<double>>(rows-1,std::vector<double>(cols-1,0));
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

    void Matrix::Print() const{
        for(unsigned int i = 0; i < rows; ++i){
            for(unsigned int j = 0; j < cols; ++j){
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::ostream& operator<<(std::ostream& stream, const Matrix& matrix){
        matrix.Print();
        return stream;
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

    Matrix IdentityMatrix::Multiply(const Matrix& other) const{
        if(this->cols != other.GetRows()){
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
