#pragma once


#include "algorithm_iterations.h"
#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "qr_decomposition.h"
#include "steps.h"

#include <algorithm>
#include <string>

namespace QR_algorithm {

template<typename T>
class Manager {
public:
    void set_symmetry_mode(bool mode) {
        symmetry_mode = mode;
    }

    bool set_calculation_mode(const std::string& mode) {
        if (mode == "eigenvalues_only") {
            calculation_mode = EIGENVALUES_ONLY;
            return true;
        } else if (mode == "without_unit") {
            calculation_mode = WITHOUT_UNIT;
            return true;
        } else if (mode == "with_unit") {
            calculation_mode = WITH_UNIT;
            return true;
        } else {
            return false;
        }
    }

    bool set_hessenberg_transform(const std::string& transform) {
        if (transform == "householder") {
            hessenberg_transform = HT_HOUSEHOLDER_REFLECTION;
            return true;
        } else if (transform == "givens") {
            hessenberg_transform = HT_GIVENS_ROTATION;
            return true;
        } else {
            return false;
        }
    }

    bool set_qr_transform(const std::string& transform) {
        if (transform == "householder") {
            qr_transform = QR_HOUSEHOLDER_REFLECTION;
            return true;
        } else if (transform == "givens") {
            qr_transform = QR_GIVENS_ROTATION;
            return true;
        } else {
            return false;
        }
    }

    bool set_shift(const std::string& shift) {
        if (shift == "none") {
            shift_mode = NONE;
            return true;
        } else if (shift == "rayleigh"){
            shift_mode = RAYLEIGH;
            return true;
        } else if (shift == "wilkinson") {
            shift_mode = WILKINSON;
            return true;
        } else if (shift == "francis") {
            shift_mode = FRANCIS;
            return true;
        } else {
            return false;
        }
    }

    void set_maximum_iterations(size_t max_iterations_) {
        max_iterations = max_iterations_;
    }

    void set_each_step_zeros_mode(bool make_each_step_zeros_) {
        make_each_step_zeros = make_each_step_zeros_;
    }

    void set_accurance (double accurance_) {
        if (accurance_ >= 0) {
            accurance = accurance_;
        } else {
            accurance = 0;
        }
    }

    void set_pseudo_shur_mode (bool pseudo_shur_) {
        pseudo_shur = pseudo_shur_;
    }

    void shur_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void shur_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void qr_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void qr_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    size_t svd_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* U, std::vector<T>* singular_values, Eigen::MatrixX<T>* Vh);

private:
    bool symmetry_mode = false;
    CALCULATION_MODE calculation_mode = WITH_UNIT;
    SHIFT shift_mode = WILKINSON;
    double accurance = 1e-4;
    size_t max_iterations = 1000;
    bool make_each_step_zeros = false;
    HESSENBERG_TRANSFORM hessenberg_transform = HT_HOUSEHOLDER_REFLECTION;
    QR_TRANSFORM qr_transform = QR_HOUSEHOLDER_REFLECTION;
    bool pseudo_shur = false;
    double eps_ignore_in_svd = 1e-8;
};

template<typename T>
void Manager<T>::shur_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    size_t check_size = center->rows();
    assert(check_size != 0);
    assert(check_size == center->cols());

    if (calculation_mode == WITH_UNIT) {
        assert(check_size == unit->rows());
        assert(check_size == unit->cols());
        *unit = Eigen::MatrixX<T>::Identity(check_size, check_size);
    }

    if (calculation_mode != WITH_UNIT) {
        make_hessenberg_form<T>(hessenberg_transform, nullptr, center);
    } else {
        make_hessenberg_form<T>(hessenberg_transform, unit, center);
    }

    if (symmetry_mode) {
        if (calculation_mode != WITH_UNIT) {
            symmetrical_iterations<T>(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                                unit, center);
        } else {
            symmetrical_iterations<T>(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                                nullptr, center);
        }
    } else {
        if (calculation_mode == WITH_UNIT) {
            shift_iterations<T>(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                        shift_mode, pseudo_shur, unit, center);
        } else {
            shift_iterations<T>(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                        shift_mode, pseudo_shur, nullptr, center);         
        }
    }
}

template<typename T>
void Manager<T>::qr_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    size_t check_size = center->rows();
    assert(check_size != 0);
    assert(check_size == unit->rows());
    assert(check_size == unit->cols());
    *unit = Eigen::MatrixX<T>::Identity(check_size, check_size);
    find_full_qr_decomposition(qr_transform, unit, center);
}

template<typename T>
void Manager<T>::shur_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    assert(matrix.rows() == center->rows());
    assert(matrix.cols() == center->cols());
    *center = matrix;
    shur_decomposition_inplace(center, unit);
}

template<typename T>
void Manager<T>::qr_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    assert(matrix.rows() == center->rows());
    assert(matrix.cols() == center->cols());
    *center = matrix;
    shur_decomposition_inplace(center, unit);
}

template<typename T>
size_t Manager<T>::svd_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* U, std::vector<T>* singular_values, Eigen::MatrixX<T>* Vh) {
    if (matrix.rows() <= matrix.cols()) {
        std::vector<T>& singular_values_ref = *singular_values;

        Eigen::MatrixX<T> center_matrix = matrix * matrix.adjoint();
        Eigen::MatrixX<T> center_matrix_conserve = center_matrix;
        size_t size = center_matrix.rows();
        Eigen::MatrixX<T> unit_matrix = Eigen::MatrixX<T>::Identity(size, size);

        set_symmetry_mode(false);
        shur_decomposition_inplace(&center_matrix, &unit_matrix);
        //std::cout << (center_matrix_conserve - unit_matrix * center_matrix * unit_matrix.adjoint()).norm() << std::endl;
        set_symmetry_mode(true);

        std::vector<std::pair<T, int>> sing_values_with_indeces(size);
        for (int i = 0; i < size; i++) {
            sing_values_with_indeces[i].first = center_matrix(i, i);
            sing_values_with_indeces[i].second = i;
        }
        std::sort(sing_values_with_indeces.begin(), sing_values_with_indeces.end(), [](std::pair<T, size_t> a, std::pair<T, size_t> b){
            return abs(a.first) > abs(b.first);
        });     

        for (int i = 0; i < size; i++) {
            singular_values_ref[i] = sqrt(sing_values_with_indeces[i].first);
            U->col(i) = unit_matrix.col(sing_values_with_indeces[i].second);
        }

        Eigen::MatrixX<T> res_center = Eigen::MatrixX<T>::Zero(size, size);

        for (int i = 0; i < size; i++) {
            res_center(i, i) = singular_values_ref[i];
        }

        *Vh = U->adjoint() * matrix;
        const Eigen::RowVectorX<T> zero_row = Eigen::MatrixX<T>::Zero(1, Vh->cols());
        const Eigen::VectorX<T> zero_col = Eigen::MatrixX<T>::Zero(U->rows(), 1);

        int answer = 0;
        for (int i = 0; i < size; i++) {
            if (abs(singular_values_ref[i]) < eps_ignore_in_svd) {
                Vh->row(i) = zero_row;
                U->col(i) = zero_col;
            } else {
                Vh->row(i) /= singular_values_ref[i];
                answer = i + 1;
            }
        } 
        return answer;
    } else {
        Eigen::MatrixX<T> Vhextra = Vh->adjoint();
        Eigen::MatrixX<T> Uextra  = U->adjoint();

        size_t answer = svd_decomposition(matrix.adjoint(), &Vhextra, singular_values, &Uextra);

        *Vh = Vhextra.adjoint();
        *U  = Uextra.adjoint();
        return answer;

    }
}

}