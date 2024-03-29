#pragma once

#include "givens_rotation.h"
#include "hessenberg_form.h"
#include "householder_reflection.h"

namespace QR_algorithm {

enum CALCULATION_MODE { WITH_UNIT, WITHOUT_UNIT, EIGENVALUES_ONLY };

template <typename T>
void sub_diag(T shift, size_t lef, size_t rig, Eigen::MatrixX<T> *center) {
  auto &center0 = *center;
  for (int i = lef; i <= rig; i++) {
    center0(i, i) -= shift;
  }
  return;
}

template <typename T>
void givens_step(bool make_each_step_zeros, size_t lef, size_t rig,
                 CALCULATION_MODE cm, Eigen::MatrixX<T> *unit,
                 Eigen::MatrixX<T> *center) {
  if (rig == lef) {
    return;
  }

  auto &center0 = *center;
  typename std::vector<Givens_rotation<T>> rotates;
  size_t sz = center0.rows();
  rotates.reserve(sz - 1);

  for (size_t i = lef; i < rig; i++) {
    T a = center0(i, i);
    T c = center0(i + 1, i);
    T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
    Givens_rotation<T> rotate = {i, i + 1, a / len, c / len};
    left_multiply(rotate, center);
    rotates.push_back(rotate.adjacent());
  }

  for (auto x : rotates) {
    right_multiply(x, center);
    if (cm == WITH_UNIT) {
      right_multiply(x, unit);
    }
  }
  if (make_each_step_zeros) {
    // cout << "mesz" << endl;
    fill_hessenberg_zeros_strings(center, lef, rig);
  }
}

template <typename T>
void rayleigh_step(bool make_each_step_zeros, size_t lef, size_t rig,
                   CALCULATION_MODE cm, Eigen::MatrixX<T> *unit,
                   Eigen::MatrixX<T> *center) {
  if (rig == lef) {
    return;
  }

  auto &center0 = *center;
  T shift = center0(rig, rig);
  size_t sz = center0.rows();
  sub_diag(shift, lef, rig, center);
  givens_step<T>(make_each_step_zeros, lef, rig, cm, unit, center);
  sub_diag(-shift, lef, rig, center);
}

template <typename T>
void wilkinson_step(bool make_each_step_zeros, size_t lef, size_t rig,
                    CALCULATION_MODE cm, Eigen::MatrixX<T> *unit,
                    Eigen::MatrixX<T> *center) {
  if (rig == lef) {
    return;
  }

  auto &center0 = *center;
  size_t sz = center0.rows();
  T corner = center0(rig, rig);
  T prev = center0(rig - 1, rig - 1);
  T prod = center0(rig - 1, rig) * center0(rig, rig - 1);
  T disc = (prev - corner) * (prev - corner) + T(4) * prod;
  T x1 = (prev + corner + sqrt(disc)) / T(2);
  T x2 = (prev + corner - sqrt(disc)) / T(2);
  T shift = (abs(x1 - corner) < abs(x2 - corner)) ? x1 : x2;
  sub_diag(shift, lef, rig, center);
  givens_step(make_each_step_zeros, lef, rig, cm, unit, center);
  sub_diag(-shift, lef, rig, center);
}

template <typename T>
void francis_step(bool make_each_step_zeros, size_t lef, size_t rig,
                  CALCULATION_MODE cm, Eigen::MatrixX<T> *unit,
                  Eigen::MatrixX<T> *center) {

  auto &center0 = *center;

  if (rig == lef) {
    return;
  }

  T trace = center0(rig, rig) + center0(rig - 1, rig - 1);
  T det = center0(rig - 1, rig - 1) * center0(rig, rig) -
          center0(rig - 1, rig) * center0(rig, rig - 1);

  size_t block_size = std::min(static_cast<size_t>(3), rig - lef + 1);

  Eigen::VectorX<T> hvec = center0.block(lef, lef, block_size, 1);
  Eigen::VectorX<T> hvec2 = hvec;
  if (block_size == 3) {
    hvec2(2, 0) = 0;
  }
  Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(block_size);

  e1(0) = 1;
  Eigen::VectorX<T> to_get =
      center0.block(lef, lef, block_size, block_size) * hvec - trace * hvec2 +
      det * e1;

  Householder_reflection<T> p0 = find_householder_reflector(to_get, 0);

  p0.make_shift(lef);

  left_multiply(p0, center);

  right_multiply(p0, center);

  right_multiply(p0, unit);

  for (size_t i = lef; i < rig - 1; i++) {
    block_size = std::min(static_cast<size_t>(3), rig - i);
    Eigen::VectorX<T> current_vec = center0.block(i + 1, i, block_size, 1);
    Householder_reflection<T> p = find_householder_reflector(current_vec, 0);
    p.make_shift(i + 1);

    left_multiply(p, center);
    right_multiply(p, center);
    right_multiply(p, unit);
  }

  if (make_each_step_zeros) {
    fill_hessenberg_zeros(center);
  }
}

template <typename T>
void symmetrical_step(bool make_each_step_zeros, size_t lef, size_t rig,
                      CALCULATION_MODE cm, Eigen::MatrixX<T> *unit,
                      Eigen::MatrixX<T> *center) {

  auto &center0 = *center;

  if (lef == rig) {
    return;
  }

  T d = (center0(rig - 1, rig - 1) - center0(rig, rig)) / T(2);
  T b = center0(rig, rig - 1);
  T shift = center0(rig, rig);
  if (d == T(0)) {
    shift = -abs(b);
  } else {
    T sign = d / abs(d);
    shift -= (b * b) / (d + sign * sqrt(d * d + b * b));
  }

  Eigen::VectorX<T> hvec = center0.block(lef, lef, 2, 1);
  hvec(0) -= shift;
  Givens_rotation<T> p0 = find_givens_rotations(hvec, 1)[0];
  p0.make_shift(lef);

  left_multiply(p0, center);
  right_multiply(p0.adjacent(), center);
  right_multiply(p0.adjacent(), unit);

  for (int i = lef; i < rig - 1; i++) {
    using std::cout;
    using std::endl;
    size_t block_size = 2;
    Eigen::VectorX<T> current_vec = center0.block(i + 1, i, block_size, 1);

    Givens_rotation<T> p = find_givens_rotations(current_vec, 1)[0];
    p.make_shift(i + 1);
    size_t lef_bord =
        static_cast<size_t>(std::max(static_cast<int>(lef), i - 3));
    size_t rig_bord =
        static_cast<size_t>(std::min(static_cast<int>(rig), i + 3));
    left_multiply(p, lef_bord, rig_bord, center);
    right_multiply(p.adjacent(), lef_bord, rig_bord, center);
    right_multiply(p.adjacent(), unit);
  }

  if (make_each_step_zeros) {
    fill_sym_hessenberg_zeros(center);
  }
}

} // namespace QR_algorithm