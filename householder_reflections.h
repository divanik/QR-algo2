#pragma once

#include<iostream>
#include<Eigen/Core>
#include<vector>

namespace  Householder_reflections {

/*template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;*/

using namespace std;

template<typename T>
class House_refl {

public:
    Eigen::VectorX<T> reflect_vector;
    size_t beg;
    //size_t size;    The size of matrix, been multiplied by our Householder_reflector. Seeems to be useless.
};

template<typename T>
std::ostream& operator<< (std::ostream& os, const House_refl<T>& hr) {
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
Eigen::MatrixX<T> operator*(const House_refl<T>& hou_refl, const Eigen::MatrixX<T>& matr) {
    auto answer = matr;
    return (left_multiply(answer, hou_refl));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const House_refl<T>& hou_refl) {
    auto answer = matr;
    return (right_multiply(answer, hou_refl));
}


template<typename T>
Eigen::MatrixX<T>& left_multiply(Eigen::MatrixX<T>& matr, const House_refl<T>& hou_refl) {
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    Eigen::MatrixX<T> subrows = matr.block(beg, 0, size, matr.cols());
    const Eigen::VectorX<T>& u = hou_refl.reflect_vector;
    subrows = u.adjoint() * subrows;
    subrows = u * subrows;
    for (auto i = 0; i < size; i++) {
        matr.row(beg + i) -= 2 * subrows.row(i);
    }
    return matr;
}

template<typename T>
Eigen::MatrixX<T>& right_multiply(Eigen::MatrixX<T>& matr, const House_refl<T>& hou_refl) {
    auto beg = hou_refl.beg;
    auto size = hou_refl.reflect_vector.size();
    Eigen::MatrixX<T> subcols = matr.block(0, beg, matr.rows(), size);
    const Eigen::VectorX<T>& u = hou_refl.reflect_vector;
    subcols = subcols * u;
    subcols = subcols * u.adjoint();
    for (auto i = 0; i < size; i++) {
        matr.col(beg + i) -= 2 * subcols.col(i);
    }
    return matr;
}

template<typename T>
Eigen::VectorX<T> find_reflector(Eigen::VectorX<T>& object) {
    Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(object.size());
    e1(0) = T(1);
    auto x1 = object(0);
    T sign;
    if (abs(x1) == T(0))  {                       // (abs(x1) < eps)
        sign = 1;
    } else {
        sign = x1 / abs(x1);
    }
    Eigen::VectorX<T> num = object - object.norm() * sign * e1;
    if (num.norm() != T(0)) {                     // (abs(x1) < eps)
        num /= num.norm();
    }
    return num;
}

template<typename T>
void make_hessenberg_form(Eigen::MatrixX<T>& unit, Eigen::MatrixX<T>& center) {
    size_t size = center.rows();

    for (size_t i = 0; i < size - 2; i++) {
        cout << "ok" << endl;
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, size - i - 1, 1);

        House_refl<T> cur_refl = {find_reflector(current_vec), i + 1};
        left_multiply(center, cur_refl);

        right_multiply(center, cur_refl);

        right_multiply(unit, cur_refl);
    }
    /*for (size_t i = 0; i < k - 2; i++) {
        for (size_t j = i + 2; j < k; j++) {
            center(j, i) = 0;
        }
    }*/
}





}