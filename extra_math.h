#pragma once

#include<complex>

namespace QR_algorithm {

double conj (double a);

bool has_root_with_discriminant (double det);

bool has_root_with_discriminant (std::complex<double> det);

}