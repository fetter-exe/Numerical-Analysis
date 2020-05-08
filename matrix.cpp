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
        if(rows == cols){
            DefineDeterminant();
        }
    }

    Matrix::Matrix(std::vector<std::vector<double>>& inputData){
        rows = inputData.size();
        cols = inputData[0].size();
        data.resize(rows);
        for(unsigned int i=0; i<rows ; ++i){ // Fills the Matrix
            copy(inputData[i].begin(), inputData[i].end(), back_inserter(data[i]));
        }
        if(rows == cols){
            DefineDeterminant();
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
                data[i].push_back((double)rand());
            }
        }
        if(rows == cols){
            DefineDeterminant();
        }
    }

    Matrix::Matrix(int numberOfRows,int numberOfCols, double minLimit, double maxLimit){
        if(numberOfRows < 0 || numberOfCols < 0){
            throw std::invalid_argument("negative arguments are not allowed");
        }
        rows = numberOfRows;
        cols = numberOfCols;
        data.resize(rows);
        srand(time(NULL)); // Guarantees a random seed
        for(unsigned int i=0; i<rows; ++i){ // Fills the Matrix
            for(unsigned int j=0; j<cols; ++j){
                data[i].push_back(minLimit + (double)rand()/maxLimit);
            }
        }
        if(rows == cols){
            DefineDeterminant();
        }
    }

    Matrix::Matrix(int numberOfRows, int numberOfCols, std::string numbers){
        std::stringstream ss(numbers);
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
                }else{ // if eof has been reached, fill the remaining data with 0
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

    double Matrix::Determinant() const{
        if(rows != cols){
            throw std::length_error("non-square matrixes can't have a determinant");
        }
        return this->determinant;
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

    void Matrix::DefineDeterminant(){
        double det = 0;

        if(rows != cols){  // Ensures only square matrixes can use this method
            throw std::length_error("non-square matrixes can't have a determinant");
        }

        if(rows == 1){  // Optimizes the function for small matrixes
            this->determinant = data[0][0];
            return;
        } else if(rows == 2){
            this->determinant = data[0][0]*data[1][1] - data[0][1]*data[1][0];
            return;
        } else if(rows == 3){
            det += data[0][0] * (data[1][1]*data[2][2] - data[2][1]*data[1][2]);
            det -= data[0][1] * (data[1][0]*data[2][2] - data[2][0]*data[1][2]);
            det += data[0][2] * (data[1][0]*data[2][1] - data[2][0]*data[1][1]);
            this->determinant = det;
            return;
        }

        // The determinant of a matrix is calculated by the following formula:
        // det(A) = a_{1,1}*det(M_{1,1}) - a_{1,2}*det(M_{1,2}) + ... + (-1^(n+1))*a_{1,n}*det(M_{1,n})
        // where a_{i,j} is the element at position [i,j] of the matrix A,
        // and M_{i,j} is a submatrix of size (m-1)x(n-1) generated by
        // removing the ith row and jth column of the matrix A(mxn)

        for(unsigned int j = 0; j < cols; ++j){
            unsigned int im, jm; // Indexes for the submatrix M
            std::vector<std::vector<double>> mData; // Data for the submatrix M

            mData.resize(rows-1); // Resizes data for the submatrix M
            for(unsigned int i = 0; i < rows; ++i){
                mData[i].resize(rows-1);
            }

            // Loops through matrix A to construct matrix M by ignoring the
            // data from the first row and jth column
            im = 0;
            for(unsigned int i2 = 0; i2 < rows; ++i2){
                if(i2 != 0){
                    jm = 0;
                    for(unsigned int j2 = 0; j2 < rows; ++j2){
                        if(j2 != j){
                            mData[im][jm] = data[i2][j2];
                            ++jm;
                        }
                    }
                    ++im;
                }
            }

            // Creates submatrix M and, implicitly through its constructor,
            // calculates its determinant (note: this is a recursive process)
            Matrix m = Matrix(mData);
            det += pow(-1,j) * data[0][j] * m.Determinant();
        }

        this->determinant = det;
    }

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
        DefineDeterminant();
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
        DefineDeterminant();
    }

    Matrix IdentityMatrix::operator*(const Matrix& other) const{
        if(this->cols != other.Rows()){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        return other;
    }

    void IdentityMatrix::DefineDeterminant(){
        determinant = 1.0;
    } 