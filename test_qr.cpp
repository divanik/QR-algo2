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

    MatrixXd uni = MatrixXd::Identity(size1, size1);

    //cout << uni << endl;

    //cout << "ok" << endl;
    find_full_qr_decomposition<double>(HOUSEHOLDER_REFLECTION, &uni, &matr);
    //cout << "ok" << endl;
    double err = 0;
    for (int i = 0; i < size2; i++) {
        for (int j = i + 1; j < size1; j++) {
            err += abs(matr(j, i)) * abs(matr(j, i));
        } 
    }

    err = sqrt(err);
    cout << endl;
    cout << err << endl;
    cout << (uni * matr - matr_conserve).norm() << endl;
    size_t k = min(size1, size2);
    cout << (uni.adjoint() * uni - MatrixXd::Identity(size1, size1)).norm() << endl;

    //cout << uni << endl << endl;
    //cout << matr << endl << endl;
}

