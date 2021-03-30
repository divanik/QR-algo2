#pragma once


#include"Eigen/Core"

#include<iostream>
#include<vector>

namespace  QR_algorithm {

/*template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;*/

using namespace std;

template<typename T>
class Householder_reflection {

public:
    Eigen::VectorX<T> reflect_vector;
    size_t beg;
    //size_t size;    The size of matrix, been multiplied by our Householder_reflector. Seeems to be useless.
};

template<typename T>
std::ostream& operator<< (std::ostream& os, const Householder_reflection<T>& hr) {
    size_t sz = hr.reflect_vector.size();
    Eigen::MatrixX<T> matr = Eigen::MatrixX<T>::Identity(sz, sz);
    matr -= 2 * hr.reflect_vector * hr.reflect_vector.adjoint();
    os << "reflect vector: " <<
        hr.reflect_vector << std::endl 
        << "matrix: " << std::endl <<
        matr;
    return os;

}

template<typename T>
Eigen::MatrixX<T> operator*(const Householder_reflection<T>& hou_refl, const Eigen::MatrixX<T>& matr) {
    auto answer = matr;
    return (left_multiply(&answer, hou_refl));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const Householder_reflection<T>& hou_refl) {
    auto answer = matr;
    return (right_multiply(&answer, hou_refl));
}


template<typename T>
void left_multiply(Eigen::MatrixX<T>* matr, const Householder_reflection<T>& hou_refl) {
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    Eigen::MatrixX<T> subrows = matr->block(beg, 0, size, matr->cols());
    const Eigen::VectorX<T>& u = hou_refl.reflect_vector;
    subrows = u.adjoint() * subrows;
    subrows = u * subrows;
    for (auto i = 0; i < size; i++) {
        matr->row(beg + i) -= 2 * subrows.row(i);
    }
    return;
}

template<typename T>
void right_multiply(Eigen::MatrixX<T>* matr, const Householder_reflection<T>& hou_refl) {
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    Eigen::MatrixX<T> subcols = matr->block(0, beg, matr->rows(), size);
    const Eigen::VectorX<T>& u = hou_refl.reflect_vector;
    subcols = subcols * u;
    subcols = subcols * u.adjoint();
    for (auto i = 0; i < size; i++) {
        matr->col(beg + i) -= 2 * subcols.col(i);
    }
    return;
}

}