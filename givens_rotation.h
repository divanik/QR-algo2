#pragma once

#include "extra_math.h"

#include<iostream>
#include"Eigen/Core"
#include<vector>

namespace QR_algorithm {

template<typename T>
class Givens_rotation {

public:

    Givens_rotation(size_t fir_ind_, size_t sec_ind_, T cos_, T sin_) :    
        fir_ind(fir_ind_), sec_ind(sec_ind_), cos(cos_), sin(sin_) {
        /*if (fir_ind_ > sec_ind_) {
            //ASSERT("Wrong arguments order in givens rotation!");
        }*/
    }

    Givens_rotation adjacent () {
        return {fir_ind, sec_ind, conj(cos), T(-1) * sin};
    }

    void make_shift(size_t p) {
        fir_ind += p;
        sec_ind += p;
    }

private:
    size_t fir_ind, sec_ind; //row and column with rotation matrix. fir_ind < sec_ind
    T cos, sin;

    template<typename U>
    friend Eigen::MatrixX<U> operator*(const Givens_rotation<U>& giv_rot, const Eigen::MatrixX<U>& matr);

    template<typename U>
    friend Eigen::MatrixX<U> operator*(const Eigen::MatrixX<U>& matr, const Givens_rotation<U>& giv_rot);

    template<typename U>
    friend void left_multiply( const Givens_rotation<U>& giv_rot, Eigen::MatrixX<U>* matr);

    template<typename U>
    friend void left_multiply( const Givens_rotation<U>& giv_rot, Eigen::VectorX<U>* vect);

    template<typename U>
    friend void right_multiply(const Givens_rotation<U>& giv_rot, Eigen::MatrixX<U>* matr);

    template<typename U>
    friend void left_multiply(const Givens_rotation<U>& giv_rot, size_t lef, size_t rig, 
                                        Eigen::MatrixX<U>* matr);

    template<typename U>
    friend void right_multiply(const Givens_rotation<U>& giv_rot, size_t lef, size_t rig, 
                                        Eigen::MatrixX<U>* matr);

    template<typename U>
    friend std::vector<Givens_rotation<T>> find_givens_rotations(Eigen::VectorX<U> object, size_t beginning);
    
};

template<typename T>
std::ostream& operator<< (std::ostream& os, const Givens_rotation<T>& gr) {
    os << "giv_rot: [" <<
        gr.fir_ind << ", " <<
        gr.sec_ind << ", " <<
        gr.cos << ", " <<
        gr.sin << "]";
    return os;

}

template<typename T>
Eigen::MatrixX<T> operator*(const Givens_rotation<T>& giv_rot, const Eigen::MatrixX<T>& matr) {
    auto answer = matr;
    return (left_multiply(&answer, giv_rot));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const Givens_rotation<T>& giv_rot) {
    auto answer = matr;
    return (right_multiply(&answer, giv_rot));
}


template<typename T>
void left_multiply( const Givens_rotation<T>& giv_rot, Eigen::MatrixX<T>* matr) {
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
void left_multiply( const Givens_rotation<T>& giv_rot, Eigen::VectorX<T>* vect) {
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
void right_multiply(const Givens_rotation<T>& giv_rot, Eigen::MatrixX<T>* matr) {
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
void left_multiply(const Givens_rotation<T>& giv_rot, size_t lef, size_t rig, 
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
void right_multiply(const Givens_rotation<T>& giv_rot, size_t lef, size_t rig, 
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

template<typename T>
std::vector<Givens_rotation<T>> find_givens_rotations(Eigen::VectorX<T> object, size_t beginning) {
    if (beginning == 0) {
        beginning = 1;
    }
    size_t sz = object.size();
    std::vector<Givens_rotation<T>> rotates;
    rotates.reserve(sz - 1);
    size_t p = beginning - 1;
    for (size_t i = beginning; i < sz; i++) {

        // cout << i << endl;
        T a = object(p);
        T c = object(i);
        T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
        // cout << i << endl;
        if (len != T(0)) {
            a = a / len;
            c = c / len;
        } else {
            a = T(1);
            c = T(0);
        }
        Givens_rotation<T> rotate(p, i, a, c);
        rotates.push_back(rotate);
        left_multiply(rotate, &object);
        // cout << vect << endl;
        // cout << i << endl;
        //std::cout << matr << std::endl;
    }
    return rotates;
}

}