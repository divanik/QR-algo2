#pragma once

#include "givens_rotation.h"
#include "householder_reflection.h"
#include "Eigen/Core"

#include <vector>

namespace QR_algorithm {

enum HESSENBERG_TRANSFORM {
    HT_HOUSEHOLDER_REFLECTION,
    HT_GIVENS_ROTATION
};

template<typename T>
void fill_hessenberg_zeros(Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    for (size_t j = 0; j < sz; j++) {
        for (size_t i = j + 2; i < sz; i++) {
            center0(i, j) = 0;
        }
    } 
    return;
}

template<typename T>
void fill_hessenberg_zeros_strings(Eigen::MatrixX<T>* center, size_t lef, size_t rig) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    for (size_t j = 0; j < sz; j++) {
        for (size_t i = max(lef, j + 2); i < min(rig + 1, sz); i++) {
            center0(i, j) = 0;
        }
    } 
    return;
}

template<typename T>
void fill_hessenberg_zeros_columns(Eigen::MatrixX<T>* center, size_t lef, size_t rig) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    for (size_t j = lef; j < max(rig + 1, sz); j++) {
        for (size_t i = j + 2; i < sz; i++) {
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
    if (ht == HT_HOUSEHOLDER_REFLECTION) {
        for (size_t i = 0; i < size - 2; i++) {
            Eigen::VectorX<T> current_vec = center0.block(i + 1, i, size - i - 1, 1);
            Householder_reflection<T> cur_refl = find_householder_reflector(current_vec, i + 1);

            left_multiply(cur_refl, center);
            right_multiply(cur_refl, center);
            right_multiply(cur_refl, unit); 
        }
    } else if (ht == HT_GIVENS_ROTATION) {
        for (size_t i = 0; i < size - 2; i++) {

            //std::cout << i << std::endl;

            Eigen::VectorX<T> current_vec = center0.block(i + 1, i, size - i - 1, 1);

            std::vector<Givens_rotation<T>> cur_rots = find_givens_rotations(current_vec, 1);

            //cout << "ok" << endl;

            for (auto& x : cur_rots) {
                x.make_shift(i + 1);

                left_multiply(x, center);
                right_multiply(x.adjacent(), center);
                right_multiply(x.adjacent(), unit);
            }
            
            //cout << center0 << endl << endl;
        }     
    }
}


}
