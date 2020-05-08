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

    Matrix(int numberOfRows, int numberOfCols, std::string numbers);
    // Instanciates a Matrix using the data from a string

    // Public methods
    unsigned int Rows() const;
    // Return the number of rows

    unsigned int Cols() const;
    // Return the number of cols

    virtual double Determinant() const;
    // Return the Matrix's determinant

    Matrix operator=(const Matrix& other) = delete;
    // Prevents Matrix data from being changed by the assignment operator

    double operator()(const unsigned int& I, const unsigned int& J) const;
    // Returns the Matrix's element indexed by row I and col J

    virtual Matrix operator+(const Matrix& other) const;
    // Returns the result of a addition matrix operation

    virtual Matrix operator-(const Matrix& other) const;
    // Returns the result of a subtraction matrix operation

    virtual Matrix operator*(const Matrix& other) const;
    // Returns the result of a multiplication matrix operation

    virtual Matrix operator*(const double& scalar) const;
    // Returns the result of a matrix scalar multiplication

    friend Matrix operator*(const double& scalar, const Matrix& matrix);
    // Returns the result of a matrix scalar multiplication

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

    virtual Matrix operator*(const Matrix& other) const override;
    // Optimizes the computation of a multiplication matrix operation
};


#endif // _MATRIX_H_