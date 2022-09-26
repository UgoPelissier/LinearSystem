#include <iostream>
#include "Matrices.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    
    Matrix4f M {
        {1.1,0,6,2},
        {8,0,-2,-2},
        {2,9,1,3},
        {2,1,-3,10}
    };
    
    SquareMatrix N({
        {1,1,1,1},
        {1,5,5,5},
        {1,5,14,14},
        {1,5,14,15.5}
    });
    
    /*
    SquareMatrix P, Q, L, U;
    tie(P,Q,L,U) = N.totalLUFactorization();
    N.printTotalLUDecomposition();
    cout << endl;
    
    vector<double> v = U.getPivot();
    for (size_t i(0); i<v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl;
    
    cout << endl;
    cout << N.SPD() << endl;
     */
    
    N.printCholeskyDecomposition();
    

}
