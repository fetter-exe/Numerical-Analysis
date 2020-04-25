#include "matrix.h"
#include <stdlib.h>  // srand(); frand();
#include <time.h>    // time

    Matrix::Matrix(){}

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

    unsigned int Matrix::Rows(){
        return this->rows;
    }

    unsigned int Matrix::Cols(){
        return this->cols;
    }

    double Matrix::Determinant(){
        return this->determinant;
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
        if(Order< 0){
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

    void IdentityMatrix::DefineDeterminant(){
        determinant = 1.0;
    } 