#pragma once

#include "conjugate.h"

#include<iostream>
#include"Eigen/Core"
#include<vector>

namespace QR_algorithm {

template<typename T>
class Given_rotation {

public:

    Given_rotation(fir_ind_, sec_ind_, T cos_, T sin_) :    
        fir_ind(fir_ind_), sec_ind(sec_ind_), cos(cos_), sin(sin_) {
        ;if (fir_ind_ > sec_ind_) {
            ASSERT("Wrong arguments order in givens rotation!");
        }
    }

private:
    size_t fir_ind, sec_ind; //row and column with rotation matrix. fir_ind < sec_ind
    T cos, sin;

    Given_rotation adjacent () {
        return {fir_ind, sec_ind, conj(cos), T(-1) * sin};
    }

};

template<typename T>
std::ostream& operator<< (std::ostream& os, const Given_rotation<T>& gr) {
    os << "giv_rot: [" <<
        gr.fir_ind << ", " <<
        gr.sec_ind << ", " <<
        gr.cos << ", " <<
        gr.sin << "]";
    return os;

}

template<typename T>
Eigen::MatrixX<T> operator*(const Given_rotation<T>& giv_rot, const Eigen::MatrixX<T>& matr) {
    auto answer = matr;
    return (left_multiply(&answer, giv_rot));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const Given_rotation<T>& giv_rot) {
    auto answer = matr;
    return (right_multiply(&answer, giv_rot));
}


template<typename T>
void left_multiply( const Given_rotation<T>& giv_rot, Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    const Eigen::RowVectorX<T> covec1 = matr->row(giv_rot.fir_ind);
    const Eigen::RowVectorX<T> covec2 = matr->row(giv_rot.sec_ind);
    matr->row(giv_rot.fir_ind) = ((conj(giv_rot.cos) * covec1) + (conj(giv_rot.sin) * covec2));
    matr->row(giv_rot.sec_ind) = ((giv_rot.cos * covec2) - (giv_rot.sin * covec1));
    return;
}

template <typename T>
void left_multiply( const Given_rotation<T>& giv_rot, Eigen::VectorX<T>* vect) {
    if (vect == nullptr) {
        return;
    }
    const T c1 = (*vect)(giv_rot.fir_ind);
    const T c2 = (*vect)(giv_rot.sec_ind);
    (*vect)(giv_rot.fir_ind) = ((conj(giv_rot.cos) * c1) + (conj(giv_rot.sin) * c2));
    (*vect)(giv_rot.sec_ind) = ((giv_rot.cos * c2) - (giv_rot.sin * c1));
    return;
}



template<typename T>
void right_multiply(const Given_rotation<T>& giv_rot, Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    const Eigen::VectorX<T>  vec1 = matr->col(giv_rot.fir_ind);
    const Eigen::VectorX<T>  vec2 = matr->col(giv_rot.sec_ind);
    matr->col(giv_rot.fir_ind) = (conj(giv_rot.cos) * vec1) - (giv_rot.sin * vec2);
    matr->col(giv_rot.sec_ind) = (giv_rot.cos * vec2) + (conj(giv_rot.sin) * vec1);
    return;
}


template <typename T>
void left_multiply(const Given_rotation<T>& giv_rot, size_t lef, size_t rig, 
                                        Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    auto& matr0 = *matr;
    size_t fi = giv_rot.fir_ind;
    size_t se = giv_rot.sec_ind;
    for (int i = lef; i <= rig; i++) {
        T c1 = ((conj(giv_rot.cos) * matr0(fi, i)) + (conj(giv_rot.sin) * matr0(se, i)));
        T c2 = ((giv_rot.cos * matr0(se, i)) - (giv_rot.sin * matr0(fi, i)));
        matr0(fi, i) = c1;
        matr0(se, i) = c2;
    }
    return;
}



template<typename T>
void right_multiply(const Given_rotation<T>& giv_rot, size_t lef, size_t rig, 
                                        Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    auto& matr0 = *matr;
    size_t fi = giv_rot.fir_ind;
    size_t se = giv_rot.sec_ind;
    for (int i = lef; i <= rig; i++) {
        T c1 = ((conj(giv_rot.cos) * matr0(i, fi)) - (giv_rot.sin * matr0(i, se)));
        T c2 = ((giv_rot.cos * matr0(i, se)) + (conj(giv_rot.sin) * matr0(i, fi)));
        matr0(i, fi) = c1;
        matr0(i, se) = c2;
    }
    return;
}

}