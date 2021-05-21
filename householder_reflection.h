#pragma once


#include"Eigen/Core"

#include<iostream>
#include<vector>

namespace  QR_algorithm {

/*template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;*/

template<typename T>
class Householder_reflection {

private:
    Eigen::VectorX<T> reflect_vector;
    size_t beg;

public:
    Householder_reflection(Eigen::VectorX<T> reflect_vector_, size_t beg_) : 
        reflect_vector(reflect_vector_), beg(beg_) {
        T norm = reflect_vector.norm();
        if (norm != static_cast<T>(0)) {
            reflect_vector /= norm;
        }
    }
    
    void make_shift(int sh) {
        int begin = (int) beg;
        begin += (int)sh;
        if (begin < 0) {
            this->beg = 0;
        } else {
            this->beg = begin;
        }
    }

    template<typename U>
    friend Eigen::MatrixX<U> operator*(const Householder_reflection<T>& hou_refl, const Eigen::MatrixX<U>& matr);

    template<typename U>
    friend Eigen::MatrixX<U> operator*(const Eigen::MatrixX<U>& matr, const Householder_reflection<U>& hou_refl);

    template<typename U>
    friend void left_multiply(const Householder_reflection<U>& hou_refl, Eigen::MatrixX<U>* matr);


    template<typename U>
    friend void right_multiply(const Householder_reflection<U>& hou_refl, Eigen::MatrixX<U>* matr);


    template<typename U>
    friend Householder_reflection<U> find_householder_reflector(const Eigen::VectorX<U>& object, 
                            size_t beginning);

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
void left_multiply(const Householder_reflection<T>& hou_refl, Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    assert(beg + size <= matr->rows());
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
void right_multiply(const Householder_reflection<T>& hou_refl, Eigen::MatrixX<T>* matr) {
    if (matr == nullptr) {
        return;
    }
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    assert(beg + size <= matr->cols());
    Eigen::MatrixX<T> subcols = matr->block(0, beg, matr->rows(), size);
    const Eigen::VectorX<T>& u = hou_refl.reflect_vector;
    subcols = subcols * u;
    subcols = subcols * u.adjoint();
    for (auto i = 0; i < size; i++) {
        matr->col(beg + i) -= 2 * subcols.col(i);
    }
    return;
}

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

}