#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"
#include "algorithm_iterations.h"
#include "Eigen/Core"

#include <complex>
using namespace Eigen;
using namespace QR_algorithm;

using comp = complex<double>;

int main() {
    int size;
    cin >> size;
    MatrixX<double> matr0 = MatrixX<double>::Random(size, size);
    MatrixX<double> matr = MatrixX<double>::Zero(size, size);

    /*for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            matr0(i, j) = (matr0(i, j) + matr0(j, i)) / 2;
            matr0(j, i) = matr0(i, j);
        }
    }*/
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matr(i, j) = matr0(i, j);
        }
    }

    //cout << matr << endl << endl;

    auto matr_conserve = matr;

    MatrixX<double> uni = MatrixX<double>::Identity(size, size);

    make_hessenberg_form<double>(GIVENS_ROTATION, &uni, &matr);

    size_t iter;
    cin >> iter;
    SHIFT shift;
    int k;
    cin >> k;
    if (k == 1) {
        shift = RAYLEIGH;
    } else if (k == 2) {
        shift = WILKINSON;
    } else if (k == 3) {
        shift = IMPLICIT_WILKINSON;
    }

    CALCULATION_MODE cm = WITH_UNIT;

    shift_iterations<double>(iter, 1e-5, false, cm, shift, true, &uni, &matr);

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;

    cout << matr << endl << endl;
    double err = 0;
    for (int i = 1; i < size; i++) {
        err += abs(matr(i, i - 1)) * abs(matr(i, i - 1));
    }
    //cout << matr << endl;
    err = sqrt(err);
    //cout << matr << endl << endl;

    cout << err << endl;
}

