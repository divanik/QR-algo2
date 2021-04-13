#include "givens_rotation.h"
#include "householder_reflection.h"

#include<iostream>
#include<utility>

using namespace std;


namespace QR_algorithm {


template<typename T>
bool Wilkinson_double_shift (Eigen::MatrixX<T>& unit,
            Eigen::MatrixX<T>& center, double eps, size_t max_num) {
    size_t size = center.rows();
    for (size_t i = 0; i < max_num; i++) {
     //cout << i << endl;
        if (i % 100000 == 99999) {
            cout << i << endl << endl;
        }
        if ( (abs(center(size - 1, size - 2)) < eps) || (abs(center(size - 2, size - 3)) < eps) ) {
            return true;
        }
        Wilkinson_double_shift_step(unit, center);
    }
    return false;   
}

}
 