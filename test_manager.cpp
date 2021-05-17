#include "manager.h"

using type = std::complex<double>;

using namespace QR_algorithm;
using namespace Eigen;

int main() {
    Manager<type> man;

    man.set_accurance(1e-8);
    man.set_shift("wilkinson");

    int size;
    cin >> size;
    MatrixX<type> matr0 = MatrixX<double>::Random(size, size);
    MatrixX<type> uni0 = MatrixX<double>::Zero(size, size);

    auto matr_conserve = matr0;

    man.shur_decomposition_inplace(&matr0, &uni0);

    cout << (uni0 * matr0 * uni0.adjoint() - matr_conserve).norm() << endl << endl;

    cout << uni0 * uni0.adjoint() << endl << endl;

    cout << matr0 << endl << endl << uni0 << endl << endl;

    
}