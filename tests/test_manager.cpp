#include<iostream>

#include "test_svd.h"
#include "test_qr.h"
#include "test_shur.h"

int main() {

    using std::cout;
    using std::cin;
    using std::endl;

    if (!test_svd()) {
        cout << "svd failed" << endl;
    } else {
        cout << "svd succeded" << endl;  
    }

    if (!test_qr()) {
        cout << "test_qr failed" << endl;
    } else {
        cout << "test_qr succeded" << endl;
    }

    if (!test_asymmetric_shur()) {
        std::cout << "assymetric version failed" << std::endl;
    } else {
        std::cout << "assymetric version succeded" << std::endl;  
    }

    if (!test_symmetric_shur()) {
        std::cout << "symetric version failed" << std::endl;
    } else {
        std::cout << "symetric version succeded" << std::endl;  
    }
}


