#pragma once


#include "algorithm_iterations.h"
#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "qr_decomposition.h"
#include "steps.h"

#include <string>

namespace QR_algorithm {

template<typename T>
class Manager {
public:
    void set_symmetry_mode(bool mode) {
        symmetry_mode = mode;
    }

    bool change_calculation_mode(const string& mode) {
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

    bool change_hessenberg_transform(const string& transform) {
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

    bool set_shift(const string& shift) {
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
        max_iterations_set = true;
    }

    void unset_maximum_iterations() {
        max_iterations = 0;
        max_iterations_set = false;
    }

    void set_each_step_zeros_mode(bool make_each_step_zeros_) {
        make_each_step_zeros = make_each_step_zeros_;
    }

    void set_accurance (double accurance_) {
        accurance = accurance_;
    }

    void set_pseudo_shur_mode (bool pseudo_shur_) {
        pseudo_shur = pseudo_shur_;
    }



    void shur_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void shur_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void qr_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

    void qr_decomposition_inplace(Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

private:
    bool symmetry_mode = false;
    CALCULATION_MODE calculation_mode = WITH_UNIT;
    SHIFT shift_mode = RAYLEIGH;
    double accurance = 1e-4;
    bool max_iterations_set = false;
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
void shur_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    assert(matrix.rows() == center->rows());
    assert(matrix.cols() == center->cols());
    *center = matrix;
    shur_decomposition_inplace(center, unit);
}

template<typename T>
void qr_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    assert(matrix.rows() == center->rows());
    assert(matrix.cols() == center->cols());
    *center = matrix;
    shur_decomposition_inplace(center, unit);
}

/*
template<typename T>
size_t svd_decomposition(const Eigen::MatrixX<T>& matrix, Eigen::MatrixX<T>* U, std::vector<T>* singular_values, Eigen::MatrixX<T>* V) {
    if (matrix.rows() <= matrix.cols()) {
        (*V) = matrix;

        Eigen::MatrixX<T> center_matrix = matrix * matrix.adjancent();
        size_t size = center_matrix.rows();
        Eigen::MatrixX<T> unit_matrix = Eigen::Identity(size, size);

        set_symmetry_mode(true);
        shur_decomposition_inplace(&center_matrix, &unit_matrix);
        set_symmetry_mode(false);
        std::vector<std::pair<T, int>> sing_values_with_indeces(size);
        for (int i = 0; i < size; i++) {
            sing_values_with_indeces[i].first = center_matrix(i, i);
        }
        sort(sing_values_with_indeces.begin(), sing_values_with_indeces.end(), [](std::pair<T, size_t> a, std::pair<T, size_t> b){
            return abs(a.first) > abs(b.first);
        });

        for (int i = 0; i < size; i++) {
            singular_values[i] = sing_values_with_indeces[i].first;
        }

        *V = (*U).adjancent() * (*V);

        const Eigen::RowVectorX<T> zero_row = Eigen::Zero(1, V->cols());

        for (int i = 0; i < size; i++) {
            if (abs(singular_values[i]) < eps_ignore_in_svd) {
                V->row(i) = zero_row;
            } else {
                V->row(i) /= singular_values[i];
        } 
    }
}
*/


}