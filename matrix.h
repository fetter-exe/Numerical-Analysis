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

    // Protected constructors
    Matrix();
    // Required for the IdentityMatrix's constructors

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

    Matrix(int numberOfRows, int numberOfCols, std::string elements);
    // Instanciates a Matrix using the data from a string

    // Public methods
    unsigned int Rows() const;
    // Returns the number of rows

    unsigned int Cols() const;
    // Returns the number of cols

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

    void ShowMatrix();

    virtual double Determinant() const;
    // Calculates the matrix's Determinant

    virtual double DeterminantPALU() const;
    // Calculates the matrix's Determinant, uses LU decomposition

    Matrix Transpose() const;
    // Returns the Transpose of the matrix

    Matrix Inverse() const;
    // Returns the Inverse of the matrix

    Matrix InversePALU() const;
    // Returns the Inverse of the matrix, uses LU decomposition

    virtual ~Matrix();
};

class IdentityMatrix : public Matrix{
  public:

    // Public constructors
    IdentityMatrix(int Order);
    // Instanciates a IdentityMatrix with order specified by the variable Order

    IdentityMatrix(const IdentityMatrix& other);
    // Instanciates a IdentityMatrix  from a copy constructor

    virtual Matrix operator*(const Matrix& other) const override;
    // Optimizes the computation of a multiplication matrix operation

    double Determinant() const override;
    // Optimizes the determinant's calculus

    ~IdentityMatrix();
};

class VandermondeMatrix : public Matrix{
  public:

    // Public constructors
    VandermondeMatrix(std::vector<double>& inputData, int numberOfCols);
    // Instanciates a VandermondeMatrix with inputData's size of rows and cols specified by numberOfCols.

    VandermondeMatrix(std::vector<double>& inputData);
    // Instanciates a square VandermondeMatrix with order specified by inputData's size.

    double Determinant() const override;
    // Returns the matrix's determinant through optimized calculation methods.

    ~VandermondeMatrix();
};


#endif // _MATRIX_H_S