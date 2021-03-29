#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "steps.h"

namespace QR_algorithm {

template<typename T>
class Manager {
    enum SHIFT{
        NONE, 
        RAYLEIGH,
        WILKINSON,
        IMPLICIT_WILKINSON
    } shift;
    enum SYMMETRY {
        SYMMETRIC,
        GENERAL
    } symmetry;
    size_t max_iterations;
};



}