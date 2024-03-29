#pragma once

#include "givens_rotation.h"
#include "hessenberg_form.h"
#include "householder_reflection.h"
#include "shift_splitter.h"
#include "steps.h"

#include <set>
#include <vector>

namespace QR_algorithm {

enum SHIFT { NONE, RAYLEIGH, WILKINSON, FRANCIS };

template <typename T>
void givens_iterations(const size_t steps_number, double eps,
                       bool make_each_step_zeros, CALCULATION_MODE cm,
                       Eigen::MatrixX<T> *unit, Eigen::MatrixX<T> *center) {
  size_t sz = center->rows();
  for (size_t step = 0; step < steps_number; step++) {
    givens_step<T>(make_each_step_zeros, 0, sz - 1, cm, unit, center);
    double err = 0;
    for (int i = 0; i < sz - 1; i++) {
      double k = abs((*center)(i + 1, i));
      err += k * k;
    }
    err = sqrt(err);
    if (err < eps) {
      return;
    }
  }
}

template <typename T>
void shift_iterations(size_t steps_number, double eps,
                      bool make_each_step_zeros, CALCULATION_MODE cm,
                      SHIFT shift, bool pseudo_shur, Eigen::MatrixX<T> *unit,
                      Eigen::MatrixX<T> *center) {
  auto &center0 = *center;
  size_t sz = center0.rows();
  Shift_splitter<T> sh_sp(0, sz - 1, center);
  for (size_t step = 0; step < steps_number; step++) {
    // cout << step << endl;
    if (sh_sp.empty()) {
      break;
    }
    for (auto &[lef, rig] : sh_sp) {
      if (shift == NONE) {
        givens_step(make_each_step_zeros, lef, rig, cm, unit, center);
      } else if (shift == RAYLEIGH) {
        rayleigh_step(make_each_step_zeros, lef, rig, cm, unit, center);
      } else if (shift == WILKINSON) {
        wilkinson_step(make_each_step_zeros, lef, rig, cm, unit, center);
      } else if (shift == FRANCIS) {
        francis_step(make_each_step_zeros, lef, rig, cm, unit, center);
      }

      for (int i = lef; i < rig; i++) {
        if (abs(center0(i + 1, i)) < eps) {
          sh_sp.fill_splitter(i);
        }
      }

      sh_sp.split_segs(lef, rig);
    }
    sh_sp.flush_buffer(pseudo_shur);
  }
}

template <typename T>
void symmetrical_iterations(const size_t steps_number, double eps,
                            bool make_each_step_zeros, CALCULATION_MODE cm,
                            Eigen::MatrixX<T> *unit,
                            Eigen::MatrixX<T> *center) {
  auto &center0 = *center;
  size_t sz = center0.rows();
  Shift_splitter<T> sh_sp(0, sz - 1, center);
  for (size_t step = 0; step < steps_number; step++) {
    if (sh_sp.empty()) {
      break;
    }
    for (auto &[lef, rig] : sh_sp) {
      symmetrical_step(make_each_step_zeros, lef, rig, cm, unit, center);
      for (int i = lef; i < rig; i++) {
        if (std::max(abs(center0(i + 1, i)), abs(center0(i, i + 1))) < eps) {
          sh_sp.fill_splitter(i);
        }
      }
      sh_sp.split_segs(lef, rig);
    }
    sh_sp.flush_buffer(false);
  }
}

} // namespace QR_algorithm
