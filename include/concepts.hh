#pragma once

#include <cmath>
#include <concepts>

namespace MATLIB::concepts {

template <typename Type>
[[nodiscard]] Type absolute_value(Type value) {
    using std::abs;
    return abs(value);
}

template <typename Type>
[[nodiscard]] Type square_root(Type value) {
    using std::sqrt;
    return sqrt(value);
}

template <typename Type>
concept is_ring = std::regular<Type> && requires(Type value_1, Type value_2) {
    { value_1 + value_1 } -> std::same_as<Type>;
    { value_1 - value_1 } -> std::same_as<Type>;
    { value_1 * value_1 } -> std::same_as<Type>;
    { -value_1 } -> std::same_as<Type>;
    Type{0};
    Type{1};
};

template <typename Type>
concept is_field = is_ring<Type> && requires(Type value_1, Type value_2) {
    { value_1 / value_2 } -> std::same_as<Type>;
};

template <typename Type>
concept normed = requires(Type value) {
    { absolute_value(value) } -> std::same_as<Type>;
};

template <typename Type>
concept ordered = std::totally_ordered<Type>;

template <typename Type>
concept has_square_root = requires(Type value) {
    { square_root(value) } -> std::same_as<Type>;
};

}  // namespace MATLIB::concepts
