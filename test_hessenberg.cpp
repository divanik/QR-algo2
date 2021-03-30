#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "Eigen/Core"

using namespace Eigen;
using namespace QR_algorithm;

int main() {
    MatrixXd matr(4, 4);
    matr << -4.4529e-01, -1.8641e+00, -2.8109e+00,  7.2941e+00,
             8.0124e+00,  6.2898e+00,  1.2058e+01, -1.6088e+01,
             9.1334e-01,  4.0087e-01,  1.1545e+00, -3.3722e-01,
             3.0210e+00,  1.9283e+01, -1.5744e-01,  3.0010e+00;
    auto matr_conserve = matr;

    MatrixXd uni = MatrixXd::Identity(4, 4);

    QR_algorithm::make_hessenberg_form<double>(&uni, &matr, HOUSEHOLDER_REFLECTION);

    cout << matr << endl << endl;

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;
}

