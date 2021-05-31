#include "manager.h"
#include "test_qr.h"

using type = std::complex<double>;

using namespace QR_algorithm;
using namespace Eigen;

bool test_qr() {
    using std::cout;
    using std::cin;
    using std::endl;
    using Eigen::MatrixX;
    using namespace QR_algorithm;

    Manager<type> man;

    //man.set_accurance(1e-6);
    //man.set_shift("wilkinson");
    //man.set_maximum_iterations(10000);

    std::vector<std::pair<int, int>> sizes = {{3, 3}, {3, 7}, {7, 3}, {50, 50}, {10, 20}, {20, 10}, {100, 150}, {200, 100}, {150, 200}};

    //cout << "ok" << endl;    
    for (auto size : sizes) {
        size_t size1 = size.first;
        size_t size2 = size.second;
        //cout << size1 << " " << size2 << endl;
        MatrixX<type> matr = MatrixX<type>::Random(size1, size2);
        MatrixX<type> uni = MatrixX<type>::Zero(size1, size1);
        MatrixX<type> center = MatrixX<type>::Zero(size1, size2);

        man.qr_decomposition(matr, &center, &uni);

        double comp_err = (uni * center - matr).norm();

        double unit_err = (uni * uni.adjoint() - MatrixX<type>::Identity(size1, size1)).norm();

        double res_err = 0;

        for (size_t i = 1; i < size1; i++) {
            for (size_t j = 0; j < std::min(i, size2); j++) {
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