#include "extra_math.h"

namespace QR_algorithm {

double conj (double a) {
    return a;
}

bool has_root_with_discriminant(double det) {
    return det >= 0;
}

bool has_root_with_discriminant (std::complex<double> det) {
    return true;
}

}