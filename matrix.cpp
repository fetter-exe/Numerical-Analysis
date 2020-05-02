#include "matrix.h"
#include <stdlib.h>  // srand(); frand();
#include <time.h>    // time;
#include <memory>    // shared_ptr

    Matrix::Matrix(){}

    Matrix::Matrix(const Matrix& other){
        this->rows = other.rows;
        this->cols = other.cols;
        data.resize(rows);
        for(unsigned int i=0; i<rows ; ++i){ // Fills the Matrix
            copy(other.data[i].begin(), other.data[i].end(), back_inserter(this->data[i]));
        }
        DefineDeterminant();
    }

    Matrix::Matrix(std::vector<std::vector<double>>& inputData){
        rows = inputData.size();
        cols = inputData[0].size();
        data.resize(rows);
        for(unsigned int i=0; i<rows ; ++i){ // Fills the Matrix
            copy(inputData[i].begin(), inputData[i].end(), back_inserter(data[i]));
        }
        DefineDeterminant();
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
        DefineDeterminant();
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
        DefineDeterminant();
    }

    unsigned int Matrix::Rows() const{
        return this->rows;
    }

    unsigned int Matrix::Cols() const{
        return this->cols;
    }

    double Matrix::Determinant() const{
        return this->determinant;
    }

    Matrix Matrix::operator+(const Matrix& other) const{
        if(this->rows != other.rows || this->cols != other.cols){
            throw std::length_error("dimension mismatch"); // Checks for a required condition
        }
        std::shared_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
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
        std::shared_ptr<Matrix> returnedMatrix(new Matrix(this->rows, this->cols));
        for(unsigned int i=0; i<returnedMatrix->rows; ++i){ // Fills the Matrix
           for(unsigned int j=0; j<returnedMatrix->cols; ++j){
               returnedMatrix->data[i][j] = this->data[i][j] - other.data[i][j];
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
    // --------> Still needs implementation <--------  
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

    void IdentityMatrix::DefineDeterminant(){
        determinant = 1.0;
    } 