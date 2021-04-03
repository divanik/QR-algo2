#pragma once

#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"

namespace QR_algorithm {

template<typename T>
void sub_diag(T shift, size_t lef, size_t rig, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    for (int i = lef; i <= rig; i++) {
        center0(i, i) -= shift;
    }
    return;
}

template<typename T>
void given_step( bool make_each_step_zeros, size_t lef, size_t rig, 
                Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    typename std::vector<Given_rotation<T>> rotates;
    size_t sz = center0.rows();
    rotates.reserve(sz - 1);

    for (size_t i = lef; i <= rig; i++) {
        T a = center0(i,i);
        T c = center0(i + 1, i);
        T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
        Given_rotation<T> rotate = {
            i, i + 1,
            a / len,
            c / len
        };
        left_multiply(rotate, center);
        rotates.push_back(rotate.adjacent());
        //std::cout << matr << std::endl;
    }

    for (auto x : rotates) {
        right_multiply(x, center);
        right_multiply(x, unit);
    }
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(center);
    }
}


template<typename T>
void rayleigh_step (bool make_each_step_zeros, size_t lef, size_t rig, 
                        Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    T shift = center0(rig, rig);
    sub_diag(shift, lef, rig, center);
    given_step<T>(make_each_step_zeros, lef, rig, unit, center);
    sub_diag(-shift, lef, rig, center);
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(center);
    }   
}

/*
template<typename T>
void simple_wilkinson_step (bool make_each_step_zeros, size_t lef, size_t rig,
                                        Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = Center->rows();
    T prev   = center0(sz - 2, sz - 2);
    T corner = center0(sz - 1, sz - 1);
    T prod   = center0(sz - 2, sz - 1) * center0(sz - 1, sz - 2);
    T disc = (prev - corner) * (prev - corner) + 4 * prod;
    T x1 = (prev + corner + sqrt(disc)) / 2; 
    T x2 = (prev + corner - sqrt(disc)) / 2;
    T shift = (abs(x1 - corner) < abs(x2 - corner)) : x1 ? x2; 
    center0 -= shift * Eigen::MatrixX<T>::Identity(sz, sz);
    givens_step(Unit, Center);
    center0 += shift * Eigen::MatrixX<T>::Identity(sz, sz);    
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(center);
    }
}

//to debug
template<typename T>
void double_wilkinson_step (Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center,
                bool make_each_step_zeros = true, size_t lef = 0, size_t rig = -1) {

    auto& center0 = *center;
    size_t size = center0.rows();

    if (size <= 2) {
        return;
    }

    T trace = center0(size - 2, size - 2) + center0(size - 1, size - 1);
    T det = center0(size - 2, size - 2) * center0(size - 1, size - 1) 
                - center0(size - 2, size - 1) * center0(size - 1, size - 2);

    Eigen::VectorX<T> hvec = center0.col(0).head(3);
    Eigen::VectorX<T> hvec2 = hvec;
    hvec2(2, 0) = 0;
    Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(3);

    e1(0) = 1;
    Eigen::VectorX<T> to_get = center0.block(0, 0, 3, 3) * hvec - trace * hvec2 + det * e1;

    Householder_reflection<T> p0 = {find_householder_reflector(to_get), 0};

    left_multiply(p0, center);

    right_multiply(p0, center);

    right_multiply(p0, unit);

    for (size_t i = 0; i < size - 2; i++) {
        size_t block_size = min(static_cast<size_t>(3), size - i - 1);
        Eigen::VectorX<T> current_vec = center0.block(i + 1, i, block_size, 1);

        Householder_reflection<T> p = {find_reflector(current_vec), i + 1};
        left_multiply(p, center);

        right_multiply(p, center);

        right_multiply(p, center);
    }
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(center);
    }
}

template<typename T>
void implicit_symmetrical_step (Eigen::MatrixX<T>& unit, Eigen::MatrixX<T>& center, 
                bool make_each_step_zeros = true, size_t lef = 0, size_t rig = -1) {

    auto& center0 = *center;
    size_t size = center0.rows();

    if (size <= 1) {
        return;
    }

    T shift = (center0(size - 2, size - 2) - center0(size - 1, size - 1);
    T b = center0(size - 1, size - 2);
    if (shift == T(0)) {
        shift -= b;
    } else {
        T sign = shift / abs(shift);
        shift -= (b * b) / (shift + sign * sqrt(shift * shift + b * b));
    }

    Eigen::VectorX<T> hvec = center0.col(0).head(2);

    std::vector< Given_rotation<T> > p0 = {find_given_rotations(to_get), 1};

    left_multiply(p0[0], center);

    right_multiply(p0[0], center);

    right_multiply(p0[0], unit);

    for (size_t i = 0; i < size - 2; i++) {
        size_t block_size = 2;
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, block_size, 1);

        Householder_reflection<T> p = {find_given_rotations(current_vec), i + 1};

        left_multiply(p0[0], center);

        right_multiply(p0[0], center);

        right_multiply(p0[0], unit);
    }

    if (make_each_step_zeros) {
        fill_sym_hessenberg_zeros(Center);
    }
}
*/

}