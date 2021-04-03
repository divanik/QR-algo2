#pragma once

#include "given_rotation.h"
#include "householder_reflection.h"
#include "hessenberg_form.h"
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
void given_iterations(const size_t steps_number, T eps, bool make_each_step_zeros,
                            Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    size_t sz = center->rows();
    for (size_t step = 0; step < steps_number; step++) {
        given_step<T>(unit, center, make_each_step_zeros, 0, sz - 1);
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

template<typename T>
void simple_shift_iterations(const size_t steps_number, T eps, bool make_each_step_zeros, 
                            Eigen::MatrixX<T>* unit, Eigen::MatrixX<T>* center) {
    auto& center0 = *center;
    size_t sz = center0.rows();
    std::set<pair<int, int>> segs;
    std::set<pair<int, int>> to_ins;
    std::set<pair<int, int>> to_del;
    segs.insert({0, sz - 1});
    std::vector<int> splitters;
    splitters.reserve(sz);
    for (size_t step = 0; step < steps_number; step++) {
        for (auto& [lef, rig] : segs) { 
            cout << lef << " " << rig << endl;
            rayleigh_step(make_each_step_zeros, lef, rig, unit, center);
            cout << "ok" << endl;
            T err = 0;
            for (int i = lef; i < rig; i++) {
                if (abs(center0(i + 1, i)) < eps) {
                    splitters.push_back(i);
                }
            }
            int cur_lef = 0;
            if (splitters.size()) {
                to_del.erase({lef, rig});
                int cur_lef = lef;
                for (auto& x : splitters) {
                    if (cur_lef != x) {
                        to_ins.insert({cur_lef, x});
                    }
                    cur_lef = x + 1;
                }
                to_ins.insert({cur_lef, rig});
                splitters.clear();
            }
            err = sqrt(err);
            if (err < eps) {
                return;
            }
        }

        for (auto x : to_del) {
            segs.erase(x);
        }
        to_del.clear();

        for (auto x : to_ins) {
            segs.insert(x);
        }
        to_ins.clear();
    }
}

/*

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
        if (err < eps) {
            return;
        }
    }
}*/

}
