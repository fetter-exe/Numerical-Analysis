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
                data[i].push_back((double)rand());
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
        srand(time(NULL)); // Guarantees a random seed
        for(unsigned int i=0; i<rows; ++i){ // Fills the Matrix
            for(unsigned int j=0; j<cols; ++j){
                data[i].push_back(minLimit + (double)rand()/maxLimit);
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

    double Matrix::Determinant() const{
        // Ensures only square matrixes can use this method
        if(rows != cols){
            throw std::length_error("non-square matrixes can't have a determinant");
        }

        double det = 1, total = 1;
        unsigned int index, n = rows;
        std::vector<double> temp(rows);
        std::vector<std::vector<double>> mat = data;

        // Loop for traversing diagonal elements
        for(unsigned int k = 0; k < n; ++k){
            index = k;

            // Finding the index that has a non-zero value
            while(mat[index][k] == 0  &&  index < n){
                ++index;
            }

            // If there is non-zero element
            if(index == n){
                continue; // The determinant of matrix as zero
            }

            if(index != k){
                // Loop for swaping the index row and diagonal row
                for(unsigned int j = 0; j < n; ++j){
                    std::swap(mat[index][j], mat[k][j]);
                }
                // Determinant changes sign when shifting rows
                det *= pow(-1, index-k);
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