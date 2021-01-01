#include "hessenberg.h"

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
        //cout << step << " " << size << endl;
        T shift = center(size - 1, size - 1);
        //cout << step << endl;
        center -= shift * Eigen::MatrixX<T>::Identity(size, size);
        //cout << step << endl;
        //cout << step << endl;
        Hessenberg::single_step_hessenberg_QR(unit, center);
        //cout << step << endl;
        center += shift * Eigen::MatrixX<T>::Identity(size, size);
        if (abs(center(size - 1, size - 2)) < eps) {
            cout << center(size - 1, size - 1) << endl;
            return {true, center(size - 1, size - 1)};
        }

        //cout << center << endl;
        //cout << step << endl;
    }
    return {false, T(0)}; 
}


}