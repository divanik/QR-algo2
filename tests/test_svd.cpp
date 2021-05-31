#include "manager.h"
#include "test_svd.h"

using type = double;

bool test_svd() {
    using namespace QR_algorithm;
    using Eigen::MatrixX;
    using std::cin;
    using std::cout;
    using std::endl;

    Manager<type> man;

    man.set_accurance(1e-5);
    man.set_shift("wilkinson");
    man.set_maximum_iterations(10000);


    std::vector<std::pair<int, int>> sizes = {{3, 3}, {3, 7}, {7, 3}, {50, 50}, {10, 20}, {20, 10}, {100, 150}, {200, 100}, {150, 200}};
    

    for (auto size : sizes) {
        size_t size1 = size.first;
        size_t size2 = size.second;
        //cout << size1 << " " << size2 << endl;
        size_t sz = std::min(size1, size2);
        MatrixX<type> matr0 = MatrixX<type>::Random(size1, size2);
        MatrixX<type> U = MatrixX<type>::Zero(size1, sz);
        MatrixX<type> Vh = MatrixX<type>::Zero(sz, size2);
        std::vector<type> sing_values(sz);

        auto matr_conserve = matr0;

        man.svd_decomposition(matr_conserve, &U, &sing_values, &Vh);
        //size_t sz = min(size1, size2);
        double unit_U_err = (U.adjoint() * U - Eigen::MatrixX<type>::Identity(sz, sz)).norm(); // Checking that U is unit
        double unit_V_err = (Vh * Vh.adjoint() - Eigen::MatrixX<type>::Identity(sz, sz)).norm(); // Checking that V is unit

        for (size_t i = 0; i < sz; i++) {
            Vh.row(i) *= sing_values[i];
        }

        double comp_err = (matr_conserve - U * Vh).norm(); // Checking that SVD has been counted correctly

        if (comp_err >= 1e-2) {
            return false;
        }
        if (unit_U_err >= 1e-2) {
            return false;
        }
        if (unit_V_err >= 1e-2) {
            return false;
        }
    }
    return true;
}
