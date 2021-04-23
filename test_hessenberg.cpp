#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "Eigen/Core"

using namespace Eigen;
using namespace QR_algorithm;

int main() {
    int size;
    cin >> size;
    MatrixXd matr = MatrixXd::Random(size, size);

    auto matr_conserve = matr;

    MatrixXd uni = MatrixXd::Identity(size, size);

    make_hessenberg_form<double>(HOUSEHOLDER_REFLECTION, &uni, &matr);

    double err = 0;
    for (int i = 0; i < size - 2; i++) {
        for (int j = i + 2; j < size; j++) {
            err += abs(matr(j, i)) * abs(matr(j, i));
        } 
    }

    err = sqrt(err);
    cout << err << endl << endl;

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;
}

