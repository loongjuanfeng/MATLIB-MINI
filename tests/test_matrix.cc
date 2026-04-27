#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

#include "matlib.hh"

using namespace MATLIB;

TEST_CASE("Matrix Construction", "[Matrix][Construction]") {
    SECTION("Matrix(value, rols, cols)") {
        Matrix<int> M(5, 3, 2);
        CHECK(M.shape() == Matrix<int>::Shape{3, 2});
        CHECK(M(2, 1) == 5);
    }
}