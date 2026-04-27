#include <print>

#include "base.hh"

using namespace MATLIB;

int main() {
    Matrix<int> A(2, {3, 2}), B(3, {2, 3});
    auto C = A * B;
    std::println("{}", C);

    return 0;
}