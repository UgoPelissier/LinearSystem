#include "Matrices.hpp"

Matrix::Matrix()
    : m_row(0), m_col(0), m_matrix(0)
{}

Matrix::Matrix(size_t const& row, size_t const& col)
{
    m_row = row;
    m_col = col;
    m_matrix.resize(m_row);
    for (size_t i(0); i<m_row; i++) {
        m_matrix[i].resize(m_col);
    }
}

Matrix::Matrix(size_t const& row, size_t const& col, double const& d)
{
    m_row = row;
    m_col = col;
    m_matrix.resize(m_row);
    for (size_t i(0); i<m_row; i++) {
        m_matrix[i].resize(m_col, d);
    }
}

Matrix::Matrix(vector<vector<double>> matrix)
    : m_row(matrix.size()), m_col(matrix[0].size()), m_matrix(matrix)
{}

size_t Matrix::getRow() const
{
    return m_row;
}

size_t Matrix::getCol() const
{
    return m_col;
}

double Matrix::operator()(size_t const& i,size_t const& j) const
{
    return m_matrix[i][j];
}

Matrix Matrix::operator()(size_t const& i,size_t const& j, double const& v)
{
    m_matrix[i][j] = v;
    return *this;
}

Matrix Matrix::operator()(vector<size_t> const& i,size_t const& j) const
{
    Matrix X(i.size(), 1);
    for (size_t k(0); k<i.size(); k++) {
        X(k, 0, m_matrix[i[k]][j]);
    }
    return X;
}

Matrix Matrix::operator()(size_t const& i,vector<size_t> const& j) const
{
    Matrix X(1, j.size());
    for (size_t k(0); k<j.size(); k++) {
        X(0, k, m_matrix[i][j[k]]);
    }
    return X;
}

Matrix Matrix::operator()(vector<size_t> const& i,vector<size_t> const& j) const
{
    Matrix S(i.size(), j.size());
    for (size_t k(0); k<i.size(); k++) {
        for (size_t l(0); l<j.size(); l++) {
            S(k, l, m_matrix[i[k]][j[l]]);
        }
    }
    return S;
}

Matrix Matrix::operator+(double const& scalar) const
{
    Matrix T(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            T = T(i, j, m_matrix[j][i] + scalar);
        }
    }
    return T;
}

Matrix Matrix::operator-(double const& scalar) const
{
    Matrix T(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            T = T(i, j, m_matrix[j][i] - scalar);
        }
    }
    return T;
}

Matrix Matrix::operator*(double const& scalar) const
{
    Matrix T(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            T = T(i, j, m_matrix[j][i] * scalar);
        }
    }
    return T;
}

Matrix Matrix::operator/(double const& scalar) const
{
    Matrix T(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            T = T(i, j, m_matrix[j][i] / scalar);
        }
    }
    return T;
}

Matrix Matrix::operator+(Matrix const& B) const
{
    Matrix C(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            C = C(i, j, m_matrix[i][j]+B(i,j));
        }
    }
    return C;
}

Matrix Matrix::operator-(Matrix const& B) const
{
    Matrix C(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            C = C(i, j, m_matrix[i][j]-B(i,j));
        }
    }
    return C;
}

Matrix Matrix::operator*(Matrix const& B) const
{
    Matrix C(m_row, B.getCol(), 0);
    double sum(0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<B.getCol(); j++) {
            sum = 0;
            for (size_t k(0); k<m_col; k++) {
                sum += m_matrix[i][k]*B(k,j);
            }
            C = C(i, j, sum);
        }
    }
    return C;
}

Matrix Matrix::transpose() const
{
    Matrix T(m_row, m_col, 0);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            T = T(i, j, m_matrix[j][i]);
        }
    }
    return T;
}

bool Matrix::isEqual(Matrix const& B) const
{
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            if (m_matrix[i][j] != B(i,j)) {
                return false;
            }
        }
    }
    return true;
}

Matrix Matrix::row(size_t const& i) const
{
    vector<size_t> j;
    for (size_t k(0); k<m_col; k++) {
        j.push_back(k);
    }
    return (*this)(i,j);
}

Matrix Matrix::row(size_t const& i, size_t const& j_start, size_t const& j_end) const
{
    vector<size_t> l;
    for (size_t m(j_start); m<j_end+1; m++) {
        l.push_back(m);
    }
    return (*this)(i,l);
}

Matrix Matrix::col(size_t const& j) const
{
    vector<size_t> i;
    for (size_t k(0); k<m_row; k++) {
        i.push_back(k);
    }
    return (*this)(i,j);
}

Matrix Matrix::col(size_t const& i_start, size_t const& i_end, size_t const& j) const
{
    vector<size_t> l;
    for (size_t m(i_start); m<i_end+1; m++) {
        l.push_back(m);
    }
    return (*this)(l,j);
}

Matrix Matrix::subMatrix(size_t const& i_start, size_t const& i_end, size_t const& j_start, size_t const& j_end) const
{
    vector<size_t> m;
    for (size_t n(i_start); n<i_end+1; n++) {
        m.push_back(n);
    }
    vector<size_t> p;
    for (size_t q(j_start); q<j_end+1; q++) {
        p.push_back(q);
    }
    return (*this)(m,p);
}

void Matrix::swapRow(size_t const& k, size_t const& l)
{
    Matrix tempK = this->row(k);
    Matrix tempL = this->row(l);
    for (size_t j(0); j<m_col; j++) {
        m_matrix[k][j] = tempL(0,j);
        m_matrix[l][j] = tempK(0,j);
    }
}

void Matrix::swapCol(size_t const& k, size_t const& l)
{
    Matrix tempK = this->col(k);
    Matrix tempL = this->col(l);
    for (size_t i(0); i<m_row; i++) {
        m_matrix[i][k] = tempL(i,0);
        m_matrix[i][l] = tempK(i,0);
    }
}

pair<pair<size_t, size_t>, double> Matrix::maxAbs() const
{
    size_t k(m_row), l(m_col);
    double m(1e-10);
    for (size_t i(0); i<m_row; i++) {
        for (size_t j(0); j<m_col; j++) {
            if (abs(m_matrix[i][j])>m) {
                m = abs(m_matrix[i][j]);
                k = i;
                l = j;
            }
        }
    }
    return make_pair(make_pair(k, l), m);
}

vector<size_t> Matrix::maxSizeChar() const
{
    vector<size_t> v;
    size_t smax(0), s(0);
    for (size_t j(0); j<m_col; j++) {
        s = 0;
        smax = 0;
        for (size_t i(0); i<m_row; i++) {
            s = sizeDouble(m_matrix[i][j]);
            if (s>smax) {
                smax = s;
            }
        }
        v.push_back(s);
    }
    return v;
}

void Matrix::print() const
{
    vector<size_t> v = this->maxSizeChar();
    size_t s(0);
    for (size_t i(0); i<m_row; i++) {
        cout << "[";
        for (size_t j(0); j<m_col; j++) {
            s = sizeDouble(m_matrix[i][j]);
            if (j<m_col-1) {
                cout << string(v[j]-s,' ');
                printDouble(m_matrix[i][j]);
                cout << " ";
            } else {
                cout << string(v[j]-s,' ');
                printDouble(m_matrix[i][j]);
            }
        }
        cout << "]" << endl;
    }
}

Matrix::~Matrix()
{}

SquareMatrix::SquareMatrix()
    : Matrix(), m_N(0)
{}

SquareMatrix::SquareMatrix(size_t const& n)
    :Matrix(n, n), m_N(n)
{}

SquareMatrix::SquareMatrix(size_t const& n, double const& d)
    :Matrix(n, n, d), m_N(n)
{}

SquareMatrix::SquareMatrix(vector<vector<double>> matrix)
    : Matrix(matrix), m_N(matrix.size())
{}

SquareMatrix::SquareMatrix(bool const& identity, size_t const& n)
{
    m_row = n;
    m_col = n;
    m_N = n;
    m_matrix.resize(m_N);
    for (size_t i(0); i<m_N; i++) {
        m_matrix[i].resize(m_N);
        for (size_t j(0); j<m_N; j++) {
            if (i==j) {
                m_matrix[i][j] = 1;
            }
        }
    }
}

SquareMatrix::SquareMatrix(Matrix const& M)
    :Matrix(M), m_N(M.getCol())
{}

size_t SquareMatrix::getSize() const
{
    return m_N;
}

tuple<SquareMatrix, SquareMatrix, SquareMatrix> SquareMatrix::partialLUFactorization() const
{
    SquareMatrix P(true, m_N), L(m_N), U(*this), In(true, m_N);
    double m_ik(0);
    size_t p;
    for (size_t k(0); k<m_N-1; k++) {
        m_ik = 0;
        // double max(0);
        bool echelonForm = false;
        if (abs(U(k,k))<1e-10) {
            p = (this->col(k+1, m_N-1, k)).maxAbs().first.first;
            if (p==m_N) {
                echelonForm = true;
            }
            else {
                P.swapRow(k, p);
                L.swapRow(k, p);
                U.swapRow(k, p);
            }
        }
        if (not echelonForm) {
            for (size_t i(k+1); i<m_N; i++) {
                m_ik = U(i,k)/U(k,k);
                L(i, k, m_ik);
                for (size_t j(k); j<m_N; j++) {
                    U(i, j, U(i,j)-U(k,j)*m_ik);
                }
            }
        }
    }
    L = SquareMatrix(L+In);
    return make_tuple(P,L,U);
}

tuple<SquareMatrix, SquareMatrix, SquareMatrix, SquareMatrix> SquareMatrix::totalLUFactorization() const
{
    SquareMatrix P(true, m_N), Q(true, m_N), L(m_N), U(*this), In(true, m_N);
    double m_ik(0);
    size_t p, q;
    for (size_t k(0); k<m_N-1; k++) {
        m_ik = 0;
        bool echelonForm = false;
        if (abs(U(k,k))<1e-10) {
            tie(p,q) = (U.subMatrix(k, m_N-1, k, m_N-1)).maxAbs().first;
            p+=k;
            q+=k;
            if (p==m_N and q==m_N) {
                echelonForm = true;
            }
            else {
                P.swapRow(k, p);
                L.swapRow(k, p);
                U.swapRow(k, p);
                Q.swapCol(k, q);
                L.swapCol(k, q);
                U.swapCol(k, q);
            }
        }
        if (not echelonForm) {
            for (size_t i(k+1); i<m_N; i++) {
                m_ik = U(i,k)/U(k,k);
                L(i, k, m_ik);
                for (size_t j(k); j<m_N; j++) {
                    U(i, j, U(i,j)-U(k,j)*m_ik);
                }
            }
        }
    }
    L = SquareMatrix(L+In);
    return make_tuple(P,Q,L,U);
}

void SquareMatrix::printPartialLUDecomposition() const
{
    cout << "A = " << endl;
    this->print();
    cout << endl;
    
    SquareMatrix P, L, U;
    tie(P, L, U) = this->partialLUFactorization();
    
    cout << "P = " << endl;
    P.print();
    cout << endl;
    
    cout << "U = " << endl;
    U.print();
    cout << endl;
    
    cout << "L = " << endl;
    L.print();
}

void SquareMatrix::printTotalLUDecomposition() const
{
    cout << "A = " << endl;
    this->print();
    cout << endl;
    
    SquareMatrix P, Q, L, U;
    tie(P, Q, L, U) = this->totalLUFactorization();
    
    cout << "P = " << endl;
    P.print();
    cout << endl;
    
    cout << "Q = " << endl;
    Q.print();
    cout << endl;
    
    cout << "U = " << endl;
    U.print();
    cout << endl;
    
    cout << "L = " << endl;
    L.print();
}

Matrix SquareMatrix::backwardSubstitution(Matrix const& b) const
{
    Matrix X(m_N, 1);
    double sum(0);
    size_t i(0);
    X(m_N-1, 0, b(m_N-1, 0)/m_matrix[m_N-1][m_N-1]);
    for (size_t j(1); j<m_N+1; j++) {
        i = m_N - j;
        sum = 0;
        for (size_t k(i+1); k<m_N; k++) {
            sum += m_matrix[i][k]*X(k,0);
        }
        X(i, 0, (b(i,0)-sum)/m_matrix[i][i]);
    }
    return X;
}

Matrix SquareMatrix::forwardSubstitution(Matrix const& b) const
{
    Matrix X(m_N, 1);
    double sum(0);
    X(0, 0, b(0, 0)/m_matrix[0][0]);
    for (size_t i(1); i<m_N; i++) {
        sum = 0;
        for (size_t k(0); k<i; k++) {
            sum += m_matrix[i][k]*X(k,0);
        }
        X(i, 0, (b(i,0)-sum)/m_matrix[i][i]);
    }
    return X;
}

Matrix SquareMatrix::solvePartialLU(Matrix const& b) const
{
    SquareMatrix P, L, U;
    tie(P, L, U) = this->partialLUFactorization();
    Matrix c = P*b;
    
    Matrix Y = L.forwardSubstitution(c);
    Matrix X = U.backwardSubstitution(Y);
    
    return X;
}

Matrix SquareMatrix::solveTotalLU(Matrix const& b) const
{
    SquareMatrix P, Q, L, U;
    tie(P, Q, L, U) = this->totalLUFactorization();
    Matrix c = P*b;
    
    Matrix Y = L.forwardSubstitution(c);
    Matrix Z = U.backwardSubstitution(Y);
    Matrix X = Q*Z;
    
    return X;
}

SquareMatrix SquareMatrix::detSubMatrix(size_t const& i,size_t const& j) const
{
    SquareMatrix S(m_N-1);
    size_t m(0), n(0);
    for (size_t k(0); k<m_N; k++) {
        n = 0;
        for (size_t l(0); l<m_N; l++) {
            if (k!=i and l!=j) {
                S(m, n, m_matrix[k][l]);
                n++;
                if (n==m_N-1) {
                    m++;
                }
            }
        }
    }
    return S;
}

double SquareMatrix::det() const {
    double det(0);
    if (m_N==1) {
        return m_matrix[0][0];
    }
    else if (m_N==2) {
        return (m_matrix[0][0]*m_matrix[1][1] - m_matrix[0][1]*m_matrix[1][0]);
    }
    else {
        for (size_t i(0); i<m_N; i++) {
            if (abs(m_matrix[i][0])>1e-10) {
                SquareMatrix sub = this->detSubMatrix(i,0);
                det += pow(-1, i)*m_matrix[i][0]*(sub.det());
            }
        }
        return det;
    }
}

vector<double> SquareMatrix::getPivot() const
{
    vector<double> pivot;
    size_t j(0);
    for (size_t i(0); i<m_N; i++) {
        j = 0;
        while (m_matrix[i][j] == 0) {
            j++;
        }
        pivot.push_back(m_matrix[i][j]);
    }
    return pivot;
}

bool SquareMatrix::SPD() const
{
    SquareMatrix In(true, m_N), P, Q, L, U;
    tie(P,Q,L,U) = this->totalLUFactorization();
    if ( (not this->isEqual(this->transpose())) or (not P.isEqual(In)) or (not Q.isEqual(In)) ) {
        return false;
    }
    else {
        vector<double> v = U.getPivot();
        for (size_t i(0); i<v.size(); i++) {
            if (v[i]<1e-10) {
                return false;
            }
        }
    }
    return true;
}

SquareMatrix SquareMatrix::choleskyFactorization() const
{
    SquareMatrix L(m_N);
    double sum(0);
    if (not this->SPD()) {
        return L;
    }
    else {
        for (size_t i(0); i<m_N; i++) {
            sum = 0;
            for (size_t k(0); k<i; k++) {
                sum += pow(L(i,k), 2);
            }
            L(i, i, sqrt(m_matrix[i][i] - sum));
            for (size_t j(0); j<m_N; j++) {
                sum = 0;
                for (size_t l(0); l<i; l++) {
                    sum += L(i,l)*L(j,l);
                }
                L(j, i, ((m_matrix[i][j]-sum)/L(i,i)));
            }
        }
        return L;
    }
}

void SquareMatrix::printCholeskyDecomposition() const
{
    this->print();
    cout << endl;
    
    SquareMatrix L(this->choleskyFactorization());
    cout << "L =" << endl;
    L.print();
}

SquareMatrix::~SquareMatrix()
{}
