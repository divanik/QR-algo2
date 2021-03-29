#include "given_rotation.h"
#include "householder_reflection.h"

enum HESSENBERG_TRANSFORM {
    HOUSEHOLDER_REFLECTION,
    GIVENS_ROTATION
};

namespace QR_algorithm {



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
    Householder_reflection<T> cur_refl = {current_vec, beginning};
    return {num, beginning};
}

template<typename T>
std::vector<Given_rotation<T>> find_givens_rotations(const Eigen::VectorX<T>& vect, size_t bottom) {
    if (bottom = 0) {
        bottom = 1;
    }
    size_t sz = object.size();
    rotates.reserve(sz - 1);
    size_t p = bottom - 1;
    for (size_t i = bottom; i < sz; i++) {
        T a = vect(p);
        T c = vect(i);
        T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
        Given_rotation<T> rotate  = {p, i, T(1), T(0)};
        if (len != T(0)) {
            rotate = {
                p, i,
                a / sqrt(abs(a) * abs(a) + abs(c) * abs(c)),
                c / sqrt(abs(a) * abs(a) + abs(c) * abs(c))
            };
        }
        rotates.push_back(rotate);
        //std::cout << matr << std::endl;
    }
    return rotates;
}

template<typename T>
void make_hessenberg_form(Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center, HESSENBERG_TRANSFORM ht) {
    size_t size = center.rows();
    if (ht == HOUSEHOLDER_REFLECTION) {
        for (size_t i = 0; i < size - 2; i++) {

            Eigen::VectorX<T> current_vec = center.block(i + 1, i, size - i - 1, 1);

            cur_refl = find_householder_reflector(current_vec, i + 1);

            left_multiply(center, cur_refl);

            right_multiply(center, cur_refl);

            right_multiply(unit, cur_refl);
        }
    } else if (ht == GIVENS_ROTATION) {
        for (size_t i = 0; i < size - 2; i++) {

            Eigen::VectorX<T> current_vec = center.block(i + 1, i, size - i - 1, 1);

            std::vector<Given_rotation<T>> cur_rots = find_givens_rotations(current_vec, 1);

            for (auto& x : cur_rots) {
                left_multiply(center, x);

                right_multiply(center, x.adjacent());

                right_multiply(unit, x.adjacent());
            }

        }     
    }
}


}


/*https://patents.google.com/patent/US8473539*/