#include "manager.h"

using type = std::complex<double>;

int main() {
    using namespace QR_algorithm;
    using Eigen::MatrixX;
    using std::cin;
    using std::cout;
    using std::endl;

    Manager<type> man;

    man.set_accurance(1e-8);
    man.set_shift("wilkinson");

    int size1;
    int size2;
    cin >> size1 >> size2;
    size_t sz = std::min(size1, size2);
    MatrixX<type> matr0 = MatrixX<type>::Random(size1, size2);
    MatrixX<type> U = MatrixX<type>::Zero(size1, sz);
    MatrixX<type> Vh = MatrixX<type>::Zero(sz, size2);
    std::vector<type> sing_values(sz);

    auto matr_conserve = matr0;

    man.svd_decomposition(matr_conserve, &U, &sing_values, &Vh);
    //size_t sz = min(size1, size2);
    cout << (U.adjoint() * U - Eigen::MatrixX<type>::Identity(sz, sz)).norm() << endl << endl;
    cout << (Vh * Vh.adjoint() - Eigen::MatrixX<type>::Identity(sz, sz)).norm() << endl << endl;

    for (size_t i = 0; i < sz; i++) {
        Vh.row(i) *= sing_values[i];
    }

    cout << (matr_conserve - U * Vh).norm() << endl;

    
}


    //cout << matr0 << endl << endl << uni0 << endl << endl;