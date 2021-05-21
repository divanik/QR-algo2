#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"
#include "algorithm_iterations.h"
#include "Eigen/Core"

#include <complex>

using type = std::complex<double>;

int main() {
    using std::cout;
    using std::cin;
    using std::endl;

    int size;
    std::cin >> size;
    Eigen::MatrixX<type> matr0 = Eigen::MatrixX<type>::Random(size, size);
    Eigen::MatrixX<type> matr = Eigen::MatrixX<type>::Zero(size, size);

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

    Eigen::MatrixX<type> uni = Eigen::MatrixX<type>::Identity(size, size);

    QR_algorithm::make_hessenberg_form<type>(QR_algorithm::HT_GIVENS_ROTATION, &uni, &matr);

    size_t iter;
    std::cin >> iter;
    QR_algorithm::SHIFT shift;
    int k;
    std::cin >> k;
    if (k == 1) {
        shift = QR_algorithm::RAYLEIGH;
    } else if (k == 2) {
        shift = QR_algorithm::WILKINSON;
    } else if (k == 3) {
        shift = QR_algorithm::FRANCIS;
    }

    QR_algorithm::CALCULATION_MODE cm = QR_algorithm::WITH_UNIT;

    QR_algorithm::shift_iterations<type>(iter, 1e-5, false, cm, shift, false, &uni, &matr);



    std::cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << std::endl << std::endl;

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

