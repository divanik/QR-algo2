#pragma once

#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"

namespace QR_algorithm {

template<typename T>
class Manager {
    symmetry_mode = false;
    enum {
        WITH_UNIT,
        WITHOUT_UNIT,
        EIGENVALUES_ONLY
    } calculation_mode;
    size_t max_iterations;
    bool make_each_step_zeros = false;

};



}