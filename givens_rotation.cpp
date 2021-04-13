#include "givens_rotation.h"
#include<iostream>      
#include<Eigen/Core>

using namespace std;
using namespace Eigen;
using namespace QR_algorithm;

int main() {
    MatrixXd matr(4, 4);
    matr << -4.4529e-01, -1.8641e+00, -2.8109e+00,  7.2941e+00,
             8.0124e+00,  6.2898e+00,  1.2058e+01, -1.6088e+01,
             0.0000e+00,  4.0087e-01,  1.1545e+00, -3.3722e-01,
             0.0000e+00,  0.0000e+00, -1.5744e-01,  3.0010e+00;
    MatrixXd uni = MatrixXd::Identity(4, 4);
    auto ini = matr;
    //RowVectorXf vec = matr.row(0);
    //RowVectorXf vec2(4);
    //vec2 << 5, 9, 0.1, -13;
    //vec << 9, 0, 1.5, -7.89;
    //matr.row(2) = (3.5 * vec) + (5.6 * vec2); 
    //cout << matr << endl; 
    //cout << uni << endl; 

    //MatrixXf uni = MatrixXf::Ones(4, 4);
    //cout << matr << endl;
    /*auto kek = find_rotates(matr);
    for (auto x : kek) {
        cout << x << endl;
    }*/
    Givens_rotation<float> kek = {0, 1, 0.6, -0.8};
    //cout << left_multiply(uni, kek) << endl;

    //hessenberg_QR<double>(uni, matr, 1000);
    cout << uni * matr * uni.transpose() - ini << endl;
    cout << matr << endl;
}