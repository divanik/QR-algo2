#include "manager.h"
#include "test_shur.h"

using type = std::complex<double>;

using namespace QR_algorithm;
using namespace Eigen;

bool test_asymmetric_shur() {
    using std::cout;
    using std::cin;
    using std::endl;
    using Eigen::MatrixX;

    Manager<type> man;

    man.set_accurance(1e-6);
    man.set_shift("wilkinson");
    man.set_maximum_iterations(10000);

    std::vector<int> sizes = {3, 10, 20, 50, 100, 130, 170, 250};

    for (auto size : sizes) {
        //std::cout << size << std::endl;
        MatrixX<type> matr = MatrixX<type>::Random(size, size);
        MatrixX<type> uni = MatrixX<type>::Zero(size, size);
        MatrixX<type> center = MatrixX<type>::Zero(size, size);

        man.shur_decomposition(matr, &center, &uni);



        double comp_err = (uni * center * uni.adjoint() - matr).norm();

        double unit_err = (uni * uni.adjoint() - Eigen::MatrixX<type>::Identity(size, size)).norm();

        double res_err = 0;
        for (int i = 1; i < size; i++) {
            for (int j = 0; j < i; j++) {
                res_err += abs(center(i, j)) * abs(center(i, j));
            }
        }

        if (comp_err >= 1e-2) {
            return false;
        }
        if (unit_err >= 1e-2) {
            return false;
        }
        if (res_err >= 1e-2) {
            return false;
        }
    }
    return true;
}

bool test_symmetric_shur() {
    Manager<type> man;

    man.set_accurance(1e-5);
    man.set_shift("wilkinson");
    man.set_maximum_iterations(10000);
    man.set_symmetry_mode(true);

    std::vector<int> sizes = {3, 10, 20, 50, 100, 130, 170, 250};

    for (auto size : sizes) {
        //std::cout << size << std::endl;
        MatrixX<type> matr = MatrixX<type>::Random(size, size);
        MatrixX<type> matr2 = (matr + matr.adjoint()) / type(2);
        MatrixX<type> uni2 = MatrixX<type>::Zero(size, size);
        MatrixX<type> center2 = MatrixX<type>::Zero(size, size);

        man.shur_decomposition(matr2, &center2, &uni2);

        double comp_err = (uni2 * center2 * uni2.adjoint() - matr2).norm();

        double unit_err = (uni2 * uni2.adjoint() - Eigen::MatrixX<type>::Identity(size, size)).norm();

        double res_err = 0;

        for (int i = 1; i < size; i++) {
            for (int j = 0; j < i; j++) {
                res_err += abs(center2(i, j)) * abs(center2(i, j)) + abs(center2(j, i)) * abs(center2(j, i));
            }
        }
        if (comp_err >= 1e-3) {
            return false;
        }
        if (unit_err >= 1e-3) {
            return false;
        }
        if (res_err >= 1e-3) {
            return false;
        }
    }
    return true;
}