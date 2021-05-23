#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"
#include "algorithm_iterations.h"
#include "Eigen/Core"

#include <complex>

using type = double;

int main() {
    using std::cout;
    using std::cin;
    using std::endl;
    using namespace QR_algorithm;

    int size;
    std::cin >> size;
    Eigen::MatrixX<type> matr0 = Eigen::MatrixX<type>::Random(size, size);
    Eigen::MatrixX<type> matr = Eigen::MatrixX<type>::Zero(size, size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matr(i, j) = (matr0(i, j) + matr0(j, i)) / type(2);
        }
    }

    auto matr_conserve = matr;

    Eigen::MatrixX<type> uni = Eigen::MatrixX<type>::Identity(size, size);
    //cout << matr << endl;

    QR_algorithm::make_hessenberg_form<type>(QR_algorithm::HT_GIVENS_ROTATION, &uni, &matr);

    //cout << matr << endl;

    //cout << (matr_conserve - uni * matr * uni.adjoint()).norm() << endl;

    size_t iter;
    std::cin >> iter;   

    symmetrical_iterations<type>(iter, 1e-9, false, WITH_UNIT, &uni, &matr);

    std::cout << (uni * matr * uni.adjoint() - matr_conserve).norm() << std::endl << std::endl;

    //cout << matr << endl << endl;
    double err = 0;
    for (int i = 1; i < size; i++) {
        err += abs(matr(i, i - 1)) * abs(matr(i, i - 1));
    }

    cout << sqrt(err) << endl;
}

