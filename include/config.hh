#pragma once

#include <cstddef>

namespace MATLIB::config {

static constexpr bool kStrict = true;
static constexpr bool kAlwaysParallel = false;

static constexpr std::size_t kMinParallel = 100'000;

}  // namespace MATLIB::config