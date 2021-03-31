#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"
#include "Eigen/Core"

using namespace Eigen;
using namespace QR_algorithm;

int main() {
    int size;
    cin >> size;
    MatrixXd matr = MatrixXd::Random(size, size);

    auto matr_conserve = matr;

    MatrixXd uni = MatrixXd::Identity(size, size);

    make_hessenberg_form<double>(GIVENS_ROTATION, &uni, &matr);

    /*double err = 0;
    for (int i = 0; i < size - 2; i++) {
        for (int j = i + 2; j < size; j++) {
            err += abs(matr(j, i)) * abs(matr(j, i));
        } 
    }

    err = sqrt(err);
    cout << err << endl << endl;*/

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;

    cout << matr << endl << endl;

    for (int i = 0; i < 200; i++) {
        given_step(&uni, &matr, false);
    }

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;

    cout << matr << endl << endl;
}

