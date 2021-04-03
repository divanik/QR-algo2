#include "given_rotation.h"
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

    auto matr_conserve = matr;

    MatrixX<comp> uni = MatrixX<comp>::Identity(size, size);

    make_hessenberg_form<comp>(GIVENS_ROTATION, &uni, &matr);

    /*comp err = 0;
    for (int i = 0; i < size - 2; i++) {
        for (int j = i + 2; j < size; j++) {
            err += abs(matr(j, i)) * abs(matr(j, i));
        } 
    }

    err = sqrt(err);
    cout << err << endl << endl;*/

    //cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;

    //cout << matr << endl << endl;

    size_t iter;
    cin >> iter;
    //simple_shift_iterations<comp>(10, 1e-13, false, &uni, &matr);
    SHIFT shift;
    int k;
    cin >> k;
    if (k) {
        shift = RAYLEIGH;
    } else {
        shift = WILKINSON;
    }

    simple_shift_iterations<comp>(iter, 1e-11, false, shift, &uni, &matr);

    cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << endl << endl;
    double err = 0;
    for (int i = 1; i < size; i++) {
        err += abs(matr(i, i - 1)) * abs(matr(i, i - 1));
    }
    err = sqrt(err);
    cout << matr << endl << endl;

    cout << err << endl;
}

