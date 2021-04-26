#pragma once

#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"

namespace QR_algorithm {

template<typename T>
class Manager {
public:
    void set_symmetry_mode(bool mode) {
        symmetry_mode = mode;
    }

    bool change_calculation_mode(string& mode) {
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

    bool set_maximum_iterations(size_t max_iterations_) {
        max_iterations = max_iterations_;
        max_iterations_set = true;
        return true;
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

    void shur_decomposition(const Eigen::MatrixX<T>& matr, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit);

private:
    symmetry_mode = false;
    calculation_mode = WITH_UNIT;
    shift_mode = NONE;
    double accurance = 1e-4;
    bool max_iterations_set = false;
    size_t max_iterations;
    bool make_each_step_zeros = false;
    HESSENBERG_TRANSFORM ht = HOUSEHOLDER_REFLECTION;
    bool pseudo_shur = false;
};

template<typename T>
void Manager<T>::shur_decomposition(const Eigen::MatrixX<T>& matr, Eigen::MatrixX<T>* center, Eigen::MatrixX<T>* unit) {
    size_t check = matr.rows();
    assert(check == matr.cols());
    assert(check == center->rows());
    assert(check == center->cols());
    if (cm != WITH_UNIT) {
        assert(check == unit->rows());
        assert(check == unit->cols());
    }

    *center = matr;
    if (cm != WITH_UNIT) {
        make_hessenberg_form(ht, nullptr, &matr1);
    } else {
        size_t size = matr.rows();
        make_hessenberg_form(ht, unit, center);
    }
    if (symmetry_mode) {
        if (cm != WITH_UNIT) {
            symmetrical_iterations(max_iterations, accurance, make_each_step_zeros, cm,
                                unit, center);
        } else {
            symmetrical_iterations(max_iterations, accurance, make_each_step_zeros, cm,
                                nullptr, center);
        }
    } else {
        if (cm != WITH_UNIT) {
            shift_iterations(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                        shift, pseudo_shur, unit, center);
        } else {
            shift_iterations(max_iterations, accurance, make_each_step_zeros, calculation_mode,
                        shift, pseudo_shur, nullptr, center);         
        }
    }
}

}