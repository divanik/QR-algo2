#pragma once

#include "givens_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
#include "shift_splitter.h"
#include "steps.h"

#include <set>
#include <vector>


namespace QR_algorithm {

enum SHIFT{
    NONE, 
    RAYLEIGH,
    WILKINSON,
    IMPLICIT_WILKINSON
};


template<typename T>
void givens_iterations(const size_t steps_number, double eps, bool make_each_step_zeros,
                            Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    size_t sz = center->rows();
    for (size_t step = 0; step < steps_number; step++) {
        givens_step<T>(make_each_step_zeros, 0, sz - 1, unit, center);
        double err = 0;
        for (int i = 0; i < sz - 1; i++) {
            double k =  abs((*center)(i + 1, i));
            err += k * k;
        }
        err = sqrt(err);
        if (err < eps) {
            return;
        }
    }
}

template<typename T>
void shift_iterations(const size_t steps_number, double eps, bool make_each_step_zeros, 
                            SHIFT shift, Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    Shift_splitter sh_sp(0, sz - 1);
    for (size_t step = 0; step < steps_number; step++) {
        for (auto& [lef, rig] : sh_sp) {

            if (shift == RAYLEIGH) {
                rayleigh_step(make_each_step_zeros, lef, rig, unit, center);
            } else if (shift == WILKINSON) {
                simple_wilkinson_step(make_each_step_zeros, lef, rig, unit, center);                
            } else if (shift == IMPLICIT_WILKINSON) {
                double_wilkinson_step(make_each_step_zeros, lef, rig, unit, center);
            }

            for (int i = lef; i < rig; i++) {
                if (abs(center0(i + 1, i)) < eps) {
                    sh_sp.fill_splitter(i);
                }
            }
            sh_sp.split_segs(lef, rig);
        }
        sh_sp.flush_buffer();
        /*for (auto& [lef, rig] : sh_sp) {
            cout << lef << " " << rig << endl << endl;
            cout << center0(lef, lef) << " " << center0(lef, rig) << endl;
            cout << center0(rig, lef) << " " << center0(rig, rig) << endl;
            cout << endl;
        }
        cout << endl;*/
    }
   /* for (auto& [lef, rig] : sh_sp) {
        cout << center0(lef, lef) << " " << center0(lef, rig) << endl;
        cout << center0(rig, lef) << " " << center0(rig, rig) << endl;
        cout << endl;
        
        cout << (center0(lef, lef) - center0(rig, rig)) * (center0(lef, lef) - center0(rig, rig)) + 
                4 * center0(rig, lef) *  center0(lef, rig) << endl << (center0(lef, lef) + center0(rig, rig)) << endl;
    }
    cout << endl;*/
}


template<typename T>
void symmetrical_iterations(const size_t steps_number, double eps, bool make_each_step_zeros, 
                                Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    Shift_splitter sh_sp(0, sz - 1);
    for (size_t step = 0; step < steps_number; step++) {
        for (auto& [lef, rig] : sh_sp) {
            symmetrical_step(make_each_step_zeros, lef, rig, unit, center);         
            for (int i = lef; i < rig; i++) {
                if (max(abs(center0(i + 1, i)), abs(center0(i, i + 1)))  < eps) {
                    sh_sp.fill_splitter(i);
                }
            }
            sh_sp.split_segs(lef, rig);
        }
        sh_sp.flush_buffer();
    }
}

}
