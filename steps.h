#include<given_rotation.h>
#include<householder_reflection.h>

namespace QR_algorithm {
template<typename T>
void fill_hessenberg_zeros(Eigen::MatrixX<T>* Center) {
    size_t sz = Center->rows();
    for (size_t i = 0; i < sz; i++) {
        for (size_t j = i + 2; j < sz; j++) {
            (*Center)(i, j) = 0;
        }
    } 
    return;
}

template<typename T>
void givens_step(Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center, bool make_each_step_zeros) {
    typename std::vector<Given_rotation<T>> rotates;
    size_t sz = Center->rows();
    rotates.reserve(sz - 1);
    for (size_t i = 0; i < sz - 1; i++) {
        T a = (*Center)(i,i);
        T c = (*Center)(i + 1, i);
        T len = sqrt(abs(a) * abs(a) + abs(c) * abs(c));
        Given_rotation<T> rotate = {
            i, i + 1,
            a / len,
            c / len
        };
        left_multiply(matr, rotate);
        rotates.push_back(rotate.adjacent());
        //std::cout << matr << std::endl;
    }
    for (auto x : rotates) {
        right_multiply(Center, x);
        right_multiply(Unit, x);
    }
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(Center);
    }
    return rotates;
}

template<typename T>
void rayleigh_step (Eigen::MatrixX<T>* Unit, 
            Eigen::MatrixX<T>* Center) {
    T shift = (*Center)(size - 1, size - 1);
    Center -= shift * Eigen::MatrixX<T>::Identity(size, size);
    givens_step(Unit, Center);
    Center += shift * Eigen::MatrixX<T>::Identity(size, size); 
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(Center);
    }   
}

template<typename T>
void simple_wilkinson_step (Eigen::MatrixX<T>* Unit, 
            Eigen::MatrixX<T>* Center) {
    size_t sz = Center->rows();
    T prev   = (*Center)(sz - 2, sz - 2);
    T corner = (*Center)(sz - 1, sz - 1);
    T prod   = (*Center)(sz - 2, sz - 1) * (*Center)(sz - 1, sz - 2);
    T disc = (prev - corner) * (prev - corner) + 4 * prod;
    T x1 = (prev + corner + sqrt(disc)) / 2; 
    T x2 = (prev + corner - sqrt(disc)) / 2;
    T shift = (abs(x1 - corner) < abs(x2 - corner)) : x1 ? x2; 
    Center -= shift * Eigen::MatrixX<T>::Identity(sz, sz);
    givens_step(Unit, Center);
    Center += shift * Eigen::MatrixX<T>::Identity(sz, sz);    
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(Center);
    }
}

//to debug
template<typename T>
void double_wilkinson_step (Eigen::MatrixX<T>* Unit, Eigen::MatrixX<T>* Center) {

    size_t size = center.rows();

    if (size <= 2) {
        return;
    }

    T trace = Center(size - 2, size - 2) + Center(size - 1, size - 1);
    T det = Center(size - 2, size - 2) * Center(size - 1, size - 1) 
                - Center(size - 2, size - 1) * Center(size - 1, size - 2);

    Eigen::VectorX<T> hvec = Center.col(0).head(3);
    Eigen::VectorX<T> hvec2 = hvec;
    hvec2(2, 0) = 0;
    Eigen::VectorX<T> e1 = Eigen::VectorX<T>::Zero(3);

    e1(0) = 1;
    Eigen::VectorX<T> to_get = center.block(0, 0, 3, 3) * hvec - trace * hvec2 + det * e1;

    Householder_reflection<T> p0 = {find_householder_reflector(to_get), 0};

    left_multiply(center, p0);

    right_multiply(center, p0);

    right_multiply(unit, p0);

    for (size_t i = 0; i < size - 2; i++) {
        size_t block_size = min(static_cast<size_t>(3), size - i - 1);
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, block_size, 1);

        Householder_reflection<T> p = {find_reflector(current_vec), i + 1};
        left_multiply(center, p);

        right_multiply(center, p);

        right_multiply(unit, p);
    }
    if (make_each_step_zeros) {
        fill_hessenberg_zeros(Center);
    }
}

template<typename T>
void implicit_symmetrical_step (Eigen::MatrixX<T>& Unit, Eigen::MatrixX<T>& Center) {

    size_t size = center.rows();

    if (size <= 1) {
        return;
    }

    T shift = ((*Center)(size - 2, size - 2) - (*Center)(size - 1, size - 1);
    T b = (*Center)(size - 1, size - 2);
    if (shift == T(0)) {
        shift -= b;
    } else {
        T sign = shift / abs(shift);
        shift -= (b * b) / (shift + sign * sqrt(shift * shift + b * b));
    }

    Eigen::VectorX<T> hvec = Center.col(0).head(2);

    std::vector< Given_rotation<T> > p0 = {find_given_rotations(to_get), 1};

    left_multiply(center, p0[0]);

    right_multiply(center, p0[0]);

    right_multiply(unit, p0[0]);

    for (size_t i = 0; i < size - 2; i++) {
        size_t block_size = 2;
        Eigen::VectorX<T> current_vec = center.block(i + 1, i, block_size, 1);

        Householder_reflection<T> p = {find_given_rotations(current_vec), i + 1};

        left_multiply(center, p[0]);

        right_multiply(center, p[0]);

        right_multiply(unit, p[0]);
    }

    if (make_each_step_zeros) {
        fill_hessenberg_zeros(Center);
    }
}


}