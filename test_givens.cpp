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
    MatrixX<comp> matr = MatrixX<comp>::Random(size, size);

    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            comp p = (matr(i, j) + matr(j, i)) / static_cast<comp>(2);
            matr(i, j) = p;
            matr(j, i) = conj(p);
        }
    }

    auto matr_conserve = matr;

    MatrixX<comp> uni = MatrixX<comp>::Identity(size, size);

    make_hessenberg_form<comp>(GIVENS_ROTATION, &uni, &matr);

    size_t iter;
    cin >> iter;
    //simple_shift_iterations<comp>(10, 1e-13, false, &uni, &matr);
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

    shift_iterations<comp>(iter, 1e-11, false, shift, nullptr, &matr);

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;
    double err = 0;
    for (int i = 1; i < size; i++) {
        err += abs(matr(i, i - 1)) * abs(matr(i, i - 1));
    }
    err = sqrt(err);
    cout << matr << endl << endl;

    cout << err << endl;
}

