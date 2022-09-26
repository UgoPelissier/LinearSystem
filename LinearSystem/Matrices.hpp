#ifndef Matrices_hpp
#define Matrices_hpp

#include <stdio.h>
#include <iomanip> 
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <math.h>
#include <tuple>

using namespace std;

inline size_t sizeDouble(double const& d) {
    size_t j(0), n(0);
    string s = to_string(d);
    n = s.length();
    j = n-1;
    while ( s[j] == '0' or s[j] == '.' ) {
        if ( s[j] == '.' ) {
            return j;
        }
        j--;
    }
    j++;
    return j;
}

inline void printDouble(double const& d) {
    size_t j = sizeDouble(d);
    string s = to_string(d);
    for (size_t i(0); i<j; i++) {
        cout << s[i];
    }
}

class Matrix
{
public:
    Matrix();
    Matrix(size_t const& row, size_t const& col);
    Matrix(size_t const& row, size_t const& col, double const& d);
    Matrix(vector<vector<double>> matrix);
    
    virtual size_t getRow() const;
    virtual size_t getCol() const;
    
    double operator()(size_t const& i,size_t const& j) const;
    Matrix operator()(size_t const& i,size_t const& j, double const& v);
    Matrix operator()(vector<size_t> const& i,size_t const& j) const;
    Matrix operator()(size_t const& i,vector<size_t> const& j) const;
    Matrix operator()(vector<size_t> const& i,vector<size_t> const& j) const;
    
    Matrix operator+(double const& scalar) const;
    Matrix operator-(double const& scalar) const;
    Matrix operator*(double const& scalar) const;
    Matrix operator/(double const& scalar) const;
    
    Matrix operator+(Matrix const& B) const;
    Matrix operator-(Matrix const& B) const;
    Matrix operator*(Matrix const& B) const;
    Matrix transpose() const;
    
    virtual bool isEqual(Matrix const& B) const;
    
    virtual Matrix row(size_t const& i) const;
    virtual Matrix row(size_t const& i, size_t const& j_start, size_t const& j_end) const;
    virtual Matrix col(size_t const& j) const;
    virtual Matrix col(size_t const& i_start, size_t const& i_end, size_t const& j) const;
    virtual Matrix subMatrix(size_t const& i_start, size_t const& i_end, size_t const& j_start, size_t const& j_end) const;
    
    virtual void swapRow(size_t const& k, size_t const& l);
    virtual void swapCol(size_t const& k, size_t const& l);
    
    virtual pair<pair<size_t, size_t>, double> maxAbs() const;
    
    virtual vector<size_t> maxSizeChar() const;
    virtual void print() const;
    ~Matrix();
    
protected:
    size_t m_row;
    size_t m_col;
    vector<vector<double>> m_matrix;
};

class SquareMatrix : public Matrix
{
public:
    SquareMatrix();
    SquareMatrix(size_t const& n);
    SquareMatrix(size_t const& n, double const& d);
    SquareMatrix(vector<vector<double>> matrix);
    SquareMatrix(bool const& identity, size_t const& n);
    SquareMatrix(Matrix const& M);
    
    virtual size_t getSize() const;
    
    virtual tuple<SquareMatrix, SquareMatrix, SquareMatrix> partialLUFactorization() const;
    virtual tuple<SquareMatrix, SquareMatrix, SquareMatrix, SquareMatrix> totalLUFactorization() const;
    
    virtual void printPartialLUDecomposition() const;
    virtual void printTotalLUDecomposition() const;
    
    virtual Matrix backwardSubstitution(Matrix const& b) const;
    virtual Matrix forwardSubstitution(Matrix const& b) const;
    
    virtual Matrix solvePartialLU(Matrix const& b) const;
    virtual Matrix solveTotalLU(Matrix const& b) const;
    
    virtual SquareMatrix detSubMatrix(size_t const& i,size_t const& j) const;
    virtual double det() const;
    
    virtual vector<double> getPivot() const;
    virtual bool SPD() const;
    virtual SquareMatrix choleskyFactorization() const;
    virtual void printCholeskyDecomposition() const;
    
    ~SquareMatrix();
    
protected:
    size_t m_N;
};


#endif /* Matrices_hpp */
