#pragma once

#include<iostream>
#include<Eigen/Core>
#include<vector>

namespace QR_algorithm {

template<typename T>
class Given_rotation {

public:
    size_t fir_ind, sec_ind; //row and column with rotation matrix. fir_ind < sec_ind
    T cos, sin;

    Given_rotation adjacent () {
        return {fir_ind, sec_ind, conjugate(cos), T(-1) * sin};
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
    return (left_multiply(answer, giv_rot));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const Given_rotation<T>& giv_rot) {
    auto answer = matr;
    return (right_multiply(answer, giv_rot));
}


template<typename T>
Eigen::MatrixX<T>& left_multiply(Eigen::MatrixX<T>* matr, const Given_rotation<T>& giv_rot) {
    const Eigen::RowVectorX<T> covec1 = matr->row(giv_rot.fir_ind);
    const Eigen::RowVectorX<T> covec2 = matr->row(giv_rot.sec_ind);
    matr->row(giv_rot.fir_ind) = ((conjugate(giv_rot.cos) * covec1) + (conjugate(giv_rot.sin) * covec2));
    matr->row(giv_rot.sec_ind) = ((giv_rot.cos * covec2) - (giv_rot.sin * covec1));
    return matr;
}

template<typename T>
Eigen::MatrixX<T>& right_multiply(Eigen::MatrixX<T>* matr, const Given_rotation<T>& giv_rot) {
    const Eigen::VectorX<T>  vec1 = matr->col(giv_rot.fir_ind);
    const Eigen::VectorX<T>  vec2 = matr->col(giv_rot.sec_ind);
    matr->col(giv_rot.fir_ind) = (conjugate(giv_rot.cos) * vec1) - (giv_rot.sin * vec2);
    matr->col(giv_rot.sec_ind) = (giv_rot.cos * vec2) + (conjugate(giv_rot.sin) * vec1);
    return matr;
}

}