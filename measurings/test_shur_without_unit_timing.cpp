#include "manager.h"

#include <chrono>
#include <fstream>

using type = std::complex<double>;

using namespace QR_algorithm;
using namespace Eigen;

int main() {
    using std::cout;
    using std::cin;
    using std::endl;
    using Eigen::MatrixX;
    using micro_type = std::chrono::microseconds;
    using nano_type = std::chrono::nanoseconds;

    std::ofstream out("shur_timing.txt");

    Manager<type> man;

    man.set_accurance(1e-6);
    man.set_shift("wilkinson");
    man.set_maximum_iterations(10000);
    man.set_calculation_mode("without_unit");

    size_t tests = 60;
    size_t heated_tests = 15;

    std::vector<size_t> sizes;
    for (int i = 10; i <= 100; i += 10) {
        sizes.push_back(i);
    }

    cout << "ok" << endl;

    for (auto size : sizes) {
        double mini;
        double maxi;
        double sum = 0.0;
        double heated_sum = 0.0;
        for (size_t current_test = 0; current_test < tests; current_test++) {
            MatrixX<type> matr = MatrixX<type>::Random(size, size);
            MatrixX<type> uni = MatrixX<type>::Zero(size, size);
            MatrixX<type> center = MatrixX<type>::Zero(size, size);
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            man.shur_decomposition(matr, &center, &uni);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::chrono::duration<double> duration = end - begin;
            if (current_test) {
                mini = std::min(mini, duration.count());
                maxi = std::max(maxi, duration.count());
            } else {
                mini = duration.count();
                maxi = duration.count();
            }
            sum += duration.count();
            if (current_test > heated_tests) {
                heated_sum += duration.count();
            }
        }
        double aver = sum / tests;
        double heated_aver = heated_sum / (tests - heated_tests);
        out << size << " " << 
            mini << " " << 
            maxi << " " << 
            aver << " " << 
            heated_aver <<
            endl;
        cout << size << endl;
    }
    
}