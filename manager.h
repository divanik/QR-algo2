#include<given_rotation.h>
#include<householder_reflection.h>
#include<shifts.h>

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
       
};