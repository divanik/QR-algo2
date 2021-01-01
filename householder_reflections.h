#pragma once

#include<iostream>
#include<Eigen/Core>
#include<vector>

namespace  Householder_reflections {

/*template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;*/

template<typename T>
class House_refl {

public:
    Eigen::VectorX<T> reflect_vector;

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
    auto subtrahend = 2 * (hou_refl.reflect_vector * 
                (hou_refl.reflect_vector.adjoint() * matr));
    return (matr -= subtrahend);
}

template<typename T>
Eigen::MatrixX<T>& right_multiply(Eigen::MatrixX<T>& matr, const House_refl<T>& hou_refl) {
    auto subtrahend = 2 * ((matr * hou_refl.reflect_vector) *
                    hou_refl.reflect_vector.adjoint());
    return (matr -= subtrahend);
}

template<typename T>
void make_hessenberg_form(Eigen::MatrixX<T>& unit, Eigen::MatrixX<T>& center) {
    size_t k = center.rows();

    for (size_t i = 0; i < k - 2; i++) {
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, k - i - 1, 1);
        //std::cout << current_vec << std::endl << std::endl;
        Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(k - i - 1);
        e1(0) = 1;
        auto x1 = current_vec(0);
        T sign;
        if (abs(x1) == 0)  {                       // (abs(x1) < eps)
            sign = 1;
        } else {
            sign = x1 / abs(x1);
        }
        Eigen::VectorX<T> num = current_vec - current_vec.norm() * sign * e1;
        num /= num.norm();

        Eigen::VectorX<T> num_cont = Eigen::VectorX<T>::Zero(k);
        for (size_t j = 0; j < k - i - 1; j++) {
            num_cont[j + i + 1] = num[j];
        }

        //std::cout << num << std::endl << std::endl;

        House_refl<T> cur_refl = {num_cont};
        left_multiply(center, cur_refl);
        right_multiply(center, cur_refl);

        //std::cout << center << std::endl;

        right_multiply(unit, cur_refl);
    }
    for (auto i = 0; i < k - 2; i++) {
        for (auto j = i + 2; j < k; j++) {
            center(j, i) = 0;
        }
    }
}





}