#include "manager.h"

using type = std::complex<double>;

using namespace QR_algorithm;
using namespace Eigen;

int main() {
    using std::cout;
    using std::cin;
    using std::endl;
    using Eigen::MatrixX;

    Manager<type> man;

    man.set_accurance(1e-4);
    man.set_shift("wilkinson");
    man.set_maximum_iterations(10000);

    int size;
    cin >> size;
    MatrixX<type> matr = MatrixX<double>::Random(size, size);
    MatrixX<type> uni = MatrixX<double>::Zero(size, size);
    MatrixX<type> center = MatrixX<double>::Zero(size, size);

    man.shur_decomposition(matr, &center, &uni);



    cout << (uni * center * uni.adjoint() - matr).norm() << endl << endl;

    cout << (uni * uni.adjoint() - Eigen::MatrixX<double>::Identity(size, size)).norm() << endl << endl;

    double err = 0;
    for (int i = 1; i < size; i++) {
        err += abs(center(i, i - 1)) * abs(center(i, i - 1));
    }

    cout << err << endl;



    //cout << matr0 << endl << endl << uni0 << endl << endl;

    
}