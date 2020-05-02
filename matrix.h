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
    Matrix(const Matrix& other);
    // Instanciates a Matrix from a copy constructor

    Matrix(std::vector<std::vector<double>>& inputData);
    // Instanciates a Matrix from a two dimension vector of double

    Matrix(int numberOfRows,int numberOfCols);
    // Instanciates a Matrix with random positive numbers

    Matrix(int numberOfRows,int numberOfCols, double minLimit, double maxLimit);
    // Instanciates a Matrix with random positive numbers 
    // between the arguments minLimit and maxLimit

    // Public methods
    unsigned int Rows() const;
    // Return the number of rows

    unsigned int Cols() const;
    // Return the number of cols

    virtual double Determinant() const;
    // Return the Matrix's determinant

    virtual Matrix operator+(const Matrix& other) const;
    // Returns the result of a addition matrix operation

    virtual Matrix operator-(const Matrix& other) const;
    // Returns the result of a subtraction matrix operation

    virtual Matrix operator*(const Matrix& other) const;
    // Returns the result of a multiplication matrix operation

    void ShowContent();
};

class IdentityMatrix : public Matrix{
  protected:

    // Protected methods
    void DefineDeterminant() override;
    // Optimizes the determinant's definition

  public:

    // Public constructors
    IdentityMatrix(int Order);
    // Instanciates a IdentityMatrix with order specified by the variable Order

    IdentityMatrix(const IdentityMatrix& other);
    // Instanciates a IdentityMatrix  from a copy constructor

};


#endif // _MATRIX_H_