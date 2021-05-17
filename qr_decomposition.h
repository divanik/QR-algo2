#include "givens_rotation.h"
#include "householder_reflection.h"

#pragma once

#include "givens_rotation.h"
#include "householder_reflection.h"
#include "Eigen/Core"

#include <vector>

namespace QR_algorithm {

enum QR_TRANSFORM {
    QR_HOUSEHOLDER_REFLECTION,
    QR_GIVENS_ROTATION
};

template<typename T>
void find_full_qr_decomposition(QR_TRANSFORM qr_tr, Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    Eigen::MatrixX<T>& center0 = *center;
    size_t rows = center0.rows();
    size_t cols = center0.cols();
    if (qr_tr == QR_HOUSEHOLDER_REFLECTION) {
        for (size_t i = 0; i < min(rows, cols); i++) {
            Eigen::VectorX<T> current_vec = center0.block(i, i, rows - i, 1);
            Householder_reflection<T> cur_refl = find_householder_reflector(current_vec, 0);
            cur_refl.make_shift(i);
            left_multiply(cur_refl, center);
            right_multiply(cur_refl, unit); 
        }
    } else if (qr_tr == QR_GIVENS_ROTATION) {
        for (size_t i = 0; i < min(rows, cols); i++) {

            Eigen::VectorX<T> current_vec = center0.block(i, i, rows - i, 1);

            std::vector<Givens_rotation<T>> cur_rots = find_givens_rotations(current_vec, 1);

            for (auto& x : cur_rots) {
                x.make_shift(i);

                left_multiply(x, center);
                right_multiply(x.adjacent(), unit);
            }
        }     
    }
}

/*template<typename T>
void find_thin_qr_decomposition(QR_TRANSFORM qr_tr, Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    Eigen::MatrixX<T>& center0 = *center;
    if (rows <= cols) {
        find_full_qr(qr_tr, unit, center);
    }
}*/

}
