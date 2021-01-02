#include "hessenberg.h"
#include "householder_reflections.h"
#include "shifts.h"

#include<iostream>      
#include<Eigen/Core>

using namespace std;
using namespace Eigen;
using namespace Hessenberg;
using namespace Householder_reflections;
using namespace Shifts;

int main() {
    MatrixXd matr(4, 4);
    matr << -4.4529e-01, -1.8641e+00, -2.8109e+00,  7.2941e+00,
             8.0124e+00,  6.2898e+00,  1.2058e+01, -1.6088e+01,
             0.0000e+00,  4.0087e-01,  1.1545e+00, -3.3722e-01,
             0.0000e+00,  0.0000e+00, -1.5744e-01,  3.0010e+00;


    auto matr_conserve = matr;

    //House_refl<double> kek = {vec};

    //cout << matr << endl << endl;
    //cout << matr * kek << endl << endl;
    //cout << kek * matr << endl << endl;

    MatrixXd uni = MatrixXd::Identity(4, 4);

    /*cout << uni << endl;

    cout << uni * matr * uni.adjoint() - matr << endl;*/


//    cout << matr << endl << endl;

 //   cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl;
    
    bool kk = Wilkinson_double_shift<double>(uni, matr, 1e-10, 1'000'000);

    cout << kk << endl << endl;

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl;
    cout << matr << endl;


}