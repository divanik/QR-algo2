#pragma once

#include<iostream>
#include<Eigen/Core>
#include<vector>

namespace Hessenberg {

/*template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;*/

template<typename T>
class Given_Rot {

public:
    size_t fir_ind, sec_ind; //row and column with rotation matrix. fir_ind < sec_ind
    T cos, sin;

    //Given_Rot(size_t fir_ind_, size_t sec_ind_, T cos_, T sin)

    Given_Rot adjacent () {
        return {fir_ind, sec_ind, cos, T(-1) * sin};
    }

    /*template<typename T>
    friend std::ostream& operator<< (std::ostream& os, const Given_Rot<T>& gr);

    template<typename T>
    friend Eigen::MatrixX<T> operator*(const Given_Rot<T>& giv_rot, 
                                const Eigen::MatrixX<T>& matr);

    template<typename T>
    friend Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr,   
                                const Given_Rot<T>& giv_rot);

    template<typename T>
    friend Eigen::MatrixX<T>& left_multiply(Eigen::MatrixX<T>& matr, 
                                const Given_Rot<T>& giv_rot);

    template<typename T>
    friend Eigen::MatrixX<T>& right_multiply(Eigen::MatrixX<T>& matr, 
                                const Given_Rot<T>& giv_rot);
    template<typename T>
    friend std::vector<Given_Rot<T>> find_rotates(Eigen::MatrixX<T>& matr);

    template<typename T>
    friend void single_step_hessenberg_QR(Eigen::MatrixX<T>& Unit, 
                                Eigen::MatrixX<T>& Center);

    template<typename T>
    friend void hessenberg_QR(Eigen::MatrixX<T>& Unit, Eigen::MatrixX<T>& Center, 
                                const size_t steps_number);*/

};

template<typename T>
std::ostream& operator<< (std::ostream& os, const Given_Rot<T>& gr) {
    os << "giv_rot: [" <<
        gr.fir_ind << ", " <<
        gr.sec_ind << ", " <<
        gr.cos << ", " <<
        gr.sin << "]";
    return os;

}

template<typename T>
Eigen::MatrixX<T> operator*(const Given_Rot<T>& giv_rot, const Eigen::MatrixX<T>& matr) {
    auto answer = matr;
    return (left_multiply(answer, giv_rot));
}

template<typename T>
Eigen::MatrixX<T> operator*(const Eigen::MatrixX<T>& matr, const Given_Rot<T>& giv_rot) {
    auto answer = matr;
    return (right_multiply(answer, giv_rot));
}


template<typename T>
Eigen::MatrixX<T>& left_multiply(Eigen::MatrixX<T>& matr, const Given_Rot<T>& giv_rot) {
    const Eigen::RowVectorX<T> covec1 = matr.row(giv_rot.fir_ind);
    const Eigen::RowVectorX<T> covec2 = matr.row(giv_rot.sec_ind);
    matr.row(giv_rot.fir_ind) = ((giv_rot.cos * covec1) + (giv_rot.sin * covec2));
    matr.row(giv_rot.sec_ind) = ((giv_rot.cos * covec2) - (giv_rot.sin * covec1));
    return matr;
}

template<typename T>
Eigen::MatrixX<T>& right_multiply(Eigen::MatrixX<T>& matr, const Given_Rot<T>& giv_rot) {
    const Eigen::VectorX<T>  vec1 = matr.col(giv_rot.fir_ind);
    const Eigen::VectorX<T>  vec2 = matr.col(giv_rot.sec_ind);
    matr.col(giv_rot.fir_ind) = (giv_rot.cos * vec1) - (giv_rot.sin * vec2);
    matr.col(giv_rot.sec_ind) = (giv_rot.cos * vec2) + (giv_rot.sin * vec1);
    return matr;
}

template<typename T>
std::vector<Given_Rot<T>> find_rotates(Eigen::MatrixX<T>& matr) {
    typename std::vector<Given_Rot<T>> rotates;
    size_t k = matr.rows();
    rotates.reserve(k - 1);
    //std::cout << matr << std::endl;
    for (size_t i = 0; i < k - 1; i++) {
        T a = matr(i,i);
        T c = matr(i + 1, i);
        Given_Rot<T> rotate = {
            i, i + 1,
            a / sqrt(a * a + c * c),
            c / sqrt(a * a + c * c)
        };
        left_multiply(matr, rotate);
        rotates.push_back(rotate.adjacent());
        //std::cout << matr << std::endl;
    }
    return rotates;
}


template<typename T>
void single_step_hessenberg_QR(Eigen::MatrixX<T>& Unit, Eigen::MatrixX<T>& Center) {
    auto rotates = find_rotates(Center);
    for (auto x : rotates) {
        right_multiply(Center, x);
        right_multiply(Unit, x);
    }
}


template<typename T>
void hessenberg_QR(Eigen::MatrixX<T>& Unit, Eigen::MatrixX<T>& Center, const size_t steps_number) {
    typename Eigen::MatrixX<T> answer = Center;
    for (size_t step = 0; step < steps_number; step++) {
        single_step_hessenberg_QR(Unit, Center);
    }
}




/*

template<typename T>
pair<MatrixX<T>, MatrixX<T>> hessenberg_QR(const MatrixX<T>& Unit, const Matrix<T>& Center, const size_t steps_number) {
    Matrix<T> answer = matrix;
    for (size_t step = 0; step < steps_number; step++) {
        matrix = 
    }
} */

}