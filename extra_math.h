#pragma once

#include<complex>

namespace QR_algorithm {

inline double conj (double a) {
    return a;
}

inline bool has_root_with_discriminant(double det) {
    return det >= 0;
}

inline bool has_root_with_discriminant (std::complex<double> det) {
    return true;
}

}