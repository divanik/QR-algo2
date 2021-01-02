#include "hessenberg.h"
#include "householder_reflections.h"

#include<iostream>
#include<utility>

using namespace std;


namespace Shifts {

template<typename T>
std::pair<bool, T> Rayleigh_shift (Eigen::MatrixX<T>& unit, 
            Eigen::MatrixX<T>& center, double eps, size_t max_num) {
    size_t size = center.rows();
    if (size == 0) {
        return {false, T(0)};
    }

    if (size == 1) {
        return {true, center(0, 0)};
    }

    for (size_t step = 0; step < max_num; ++step) {
        T shift = center(size - 1, size - 1);
        center -= shift * Eigen::MatrixX<T>::Identity(size, size);
        Hessenberg::single_step_hessenberg_QR(unit, center);
        center += shift * Eigen::MatrixX<T>::Identity(size, size);
        if (abs(center(size - 1, size - 2)) < eps) {
            cout << center(size - 1, size - 1) << endl;
            return {true, center(size - 1, size - 1)};
        }
    }
    return {false, T(0)}; 
}



template<typename T>
void Wilkinson_double_shift_step (Eigen::MatrixX<T>& unit,
            Eigen::MatrixX<T>& center) {

    size_t size = center.rows();

    if (size <= 2) {
        return;
    }

    T trace = center(size - 2, size - 2) + center(size - 1, size - 1);
    T det =   center(size - 2, size - 2) * center(size - 1, size - 1) -
              center(size - 2, size - 1) * center(size - 1, size - 2);

    Eigen::VectorX<T> hvec = center.col(0).head(3);
    Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(3);

    //cout << "ok" << endl;

    e1(0) = 1;
    Eigen::VectorX<T> to_get = center.block(0, 0, 3, 3) * hvec - 2 * real(trace) * hvec.head(3) + abs(det) * abs(det) * e1;

    //cout << "ok" << endl;

    Householder_reflections::House_refl<T> p0 = {Householder_reflections::find_reflector(to_get), 0};

    Householder_reflections::left_multiply(center, p0);

    Householder_reflections::right_multiply(center, p0);

    Householder_reflections::right_multiply(unit, p0);

    //cout << "ok" << endl;

    for (size_t i = 0; i < size - 2; i++) {
        //cout << "ok" << endl;
        size_t block_size = min(static_cast<size_t>(3), size - i - 1);
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, block_size, 1);

        Householder_reflections::House_refl<T> p = {Householder_reflections::find_reflector(current_vec), i + 1};
        Householder_reflections::left_multiply(center, p);

        Householder_reflections::right_multiply(center, p);

        Householder_reflections::right_multiply(unit, p);
    }
}

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
 