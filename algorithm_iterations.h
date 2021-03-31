#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"




namespace QR_algorithm {

enum SHIFT{
    NONE, 
    RAYLEIGH,
    WILKINSON,
    IMPLICIT_WILKINSON
};


template<typename T>
void given_iterations(Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center, const size_t steps_number, T eps,
                                            bool make_each_step_zeros) {
    size_t sz = center.rows();
    for (size_t step = 0; step < steps_number; step++) {
        given_step(unit, center, make_each_step_zeros);
        T err = 0;
        for (int i = 0; i < sz - 1; i++) {
            T k =  abs((*center)(i + 1, i));
            err += k * k;
        }
        err = sqrt(err);
        if (err < eps) {
            return;
        }
    }
}

/*template<typename T>
void simple_shift_iterations(Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center, const size_t steps_number, T eps,
                                            bool make_each_step_zeros, SHIFT shift) {
    typename Eigen::MatrixX<T> answer = Center;
    size_t sz = Center.rows();
    size_t cur = sz;
    for (size_t step = 0; step < steps_number; step++) {
        if (shift == WILKINSON) {
            simple_wilkinson_step(Unit, Center, make_each_step_zeros);
        } else if (shift == RAYLEIGH) {
            rayleigh_step(Unit, Center, make_each_step_zeros);
        }
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
void implicit_iterations(Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center, const size_t steps_number, T eps,
                                            bool make_each_step_zeros) {
    typename Eigen::MatrixX<T> answer = Center;
    size_t sz = Center.rows();
    size_t cur = sz;
    for (size_t step = 0; step < steps_number; step++) {
        given_step(Unit, Center, make_each_step_zeros);
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
}*/

}
