#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"




namespace QR_algorithm {

template<typename T>
void given_rotations(Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center, const size_t steps_number, T eps) {
    typename Eigen::MatrixX<T> answer = Center;
    size_t sz = Center.rows();
    for (size_t step = 0; step < steps_number; step++) {
        single_step_hessenberg_QR(Unit, Center);
        T err = 0;
        for (int i = 0; i < sz - 1; i++) {
            T k =  abs((*Center)(i + 1, i));
            err += k * k;
        }
        err = sqrt(err);
        if (err < eps) [
            return;
        ]
    }
}

template<typename T>
void given_rotations(Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center, const size_t steps_number, T eps) {
    typename Eigen::MatrixX<T> answer = Center;
    for (size_t step = 0; step < steps_number; step++) {
        single_step_hessenberg_QR(Unit, Center);
    }
}


}
