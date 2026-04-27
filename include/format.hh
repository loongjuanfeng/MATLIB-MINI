#pragma once

#include <format>

#include "matlib.hh"

namespace std {

template <typename Scalar>
struct formatter<::MATLIB::Matrix<Scalar>> {
    constexpr auto parse(std::format_parse_context& context) {
        return context.begin();
    }

    auto format(const ::MATLIB::Matrix<Scalar>& matrix,
                std::format_context& context) const {
        auto out = context.out();
        for (std::size_t r = 0; r < matrix.shape().rows; ++r) {
            if (r == 0) {
                out = std::format_to(out, "[[");
            } else {
                out = std::format_to(out, "\n [");
            }
            for (std::size_t c = 0; c < matrix.shape().cols; ++c) {
                if (c > 0) out = std::format_to(out, ", ");
                out = std::format_to(out, "{}", matrix(r, c));
            }
            out = std::format_to(out, "]");
        }
        out = std::format_to(out, "]");
        return out;
    }
};

}  // namespace std