#pragma once

#include "given_rotation.h"
#include "householder_reflection.h"
#include "Eigen/Core"

enum HESSENBERG_TRANSFORM {
    HOUSEHOLDER_REFLECTION,
    GIVENS_ROTATION
};

using namespace std;

namespace QR_algorithm {

template<typename T>
Householder_reflection<T> find_householder_reflector(const Eigen::VectorX<T>& object, 
                            size_t beginning) {
    Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(object.size());
    e1(0) = T(1);
    auto x1 = object(0);
    T sign;
    if (abs(x1) == T(0))  {                       // (abs(x1) < eps), but seems to be stable
        sign = 1;
    } else {
        sign = x1 / abs(x1);
    }
    Eigen::VectorX<T> num = object - object.norm() * sign * e1;
    if (num.norm() != T(0)) {                     // (abs(x1) < eps)
        num /= num.norm();
    }
    Householder_reflection<T> cur_refl = {num, beginning};
    return {num, beginning};
}

template<typename T>
std::vector<Given_rotation<T>> find_givens_rotations(Eigen::VectorX<T> vect, size_t bottom) {
    if (bottom == 0) {
        bottom = 1;
    }
    size_t sz = vect.size();
    std::vector<Given_rotation<T>> rotates;
    rotates.reserve(sz - 1);
    size_t p = bottom - 1;
    for (size_t i = bottom; i < sz; i++) {

        // cout << i << endl;
        T a = vect(p);
        T c = vect(i);
        T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
        // cout << i << endl;
        Given_rotation<T> rotate  = {p, i, T(1), T(0)};
        if (len != T(0)) {
            rotate = {
                p, i,
                a / len,
                c / len
            };
        }
        rotates.push_back(rotate);
        left_multiply(rotate, &vect);
        // cout << vect << endl;
        // cout << i << endl;
        //std::cout << matr << std::endl;
    }
    return rotates;
}

template<typename T>
void fill_hessenberg_zeros(Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    for (size_t i = 0; i < sz; i++) {
        for (size_t j = i + 2; j < sz; j++) {
            center0(i, j) = 0;
        }
    } 
    return;
}

template<typename T>
void fill_sym_hessenberg_zeros(Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    for (size_t i = 0; i < sz; i++) {
        for (size_t j = i + 2; j < sz; j++) {
            center0(i, j) = 0;
            center0(j, i) = 0;
        }
    } 
    return;
}


template<typename T>
void make_hessenberg_form(HESSENBERG_TRANSFORM ht, Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    Eigen::MatrixX<T>& center0 = *center;
    size_t size = center0.rows();
    if (ht == HOUSEHOLDER_REFLECTION) {
        for (size_t i = 0; i < size - 2; i++) {

            Eigen::VectorX<T> current_vec = center0.block(i + 1, i, size - i - 1, 1);

            Householder_reflection cur_refl = find_householder_reflector(current_vec, i + 1);

            left_multiply(cur_refl, center);

            right_multiply(cur_refl, center);

            right_multiply(cur_refl, unit); 

        }
    } else if (ht == GIVENS_ROTATION) {
        for (size_t i = 0; i < size - 2; i++) {

            //std::cout << i << std::endl;

            Eigen::VectorX<T> current_vec = center0.block(i + 1, i, size - i - 1, 1);

            std::vector<Given_rotation<T>> cur_rots = find_givens_rotations(current_vec, 1);

            //cout << "ok" << endl;

            for (auto& x : cur_rots) {
                x.fir_ind += i + 1;
                x.sec_ind += i + 1;
                left_multiply(x, center);

                right_multiply(x.adjacent(), center);

                right_multiply(x.adjacent(), unit);
            }
            
            //cout << center0 << endl << endl;
        }     
    }
}


}


/*https://patents.google.com/patent/US8473539*/