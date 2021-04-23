#include "givens_rotation.h"
#include "householder_reflection.h"
#include "qr_decomposition.h"
#include "Eigen/Core"

using namespace Eigen;
using namespace QR_algorithm;

int main() {
    int size1, size2;
    cin >> size1 >> size2;
    MatrixXd matr = MatrixXd::Random(size1, size2);

    auto matr_conserve = matr;

    MatrixXd uni = MatrixXd::Identity(size2, size1);

    cout << uni << endl;

    find_full_qr_decomposition<double>(HOUSEHOLDER_REFLECTION, &uni, &matr);
    cout << "ok" << endl;
    matr = matr.block(0, 0, size2, size2);
    double err = 0;
    for (int i = 0; i < size2; i++) {
        for (int j = i + 1; j < size2; j++) {
            err += abs(matr(j, i)) * abs(matr(j, i));
        } 
    }

    err = sqrt(err);
    cout << err << endl << endl;
    cout << "kek" << endl;
    auto uni1 = uni.adjoint();
    cout << "kek" << endl;
    cout << uni1 << endl << endl;
    cout << matr << endl << endl;
    cout << "lol " << (uni1 * matr - matr_conserve).norm() << endl << endl;

    cout << (uni1.adjoint() * uni1 - MatrixXd::Identity(size2, size2)).norm() << endl;

    //cout << uni << endl << endl;
    //cout << matr << endl << endl;
}

