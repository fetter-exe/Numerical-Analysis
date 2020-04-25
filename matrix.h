#ifndef _MATRIX_H_
#define _MATRIX_H_

#include  <iostream>
#include  <vector>

class Matrix{
  protected:
    // Protected attributes
    std::vector<std::vector<double>> data;
    unsigned int rows;
    unsigned int cols;
    double determinant;

    // Protected constructors
    Matrix();
    // Required for the IdentityMatrix's constructors

    // Protected methods
    virtual void DefineDeterminant();
    // This method is called once in each constructor

  public:

    // Public constructors
    Matrix(std::vector<std::vector<double>>& inputData);
    // Instanciates a Matrix from a two dimension vector of double

    Matrix(int numberOfRows, int numberOfCols);
    // Instanciates a Matrix with random numbers

    Matrix(int numberOfRows, int numberOfCols, double minLimit, double maxLimit);
    // Instanciates a Matrix with random numbers between the
    // arguments minLimit and maxLimit

    // Public methods
    unsigned int Rows();
    // Return the number of rows

    unsigned int Cols();
    // Return the number of cols

    virtual double Determinant();
    // Return the Matrix's determinant


    void ShowContent();
};

class IdentityMatrix : public Matrix{
  public:

    // Protected constructors
    IdentityMatrix(unsigned int Order);
    // Instanciates a identity Matrix with order specified by the variable Order

    // Protected methods
    void DefineDeterminant() override;
    // Optimizes the definition of the determinant
};


#endif // _MATRIX_H_