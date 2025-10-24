#pragma once

#include <cstdint>
#include <algorithm>

template <typename T, std::uint32_t n>
struct anasir_states
{

    T x[n] = {};

    anasir_states() {};

    inline T &operator[](std::uint32_t const &i)
    {
        return x[i];
    }

    inline const T &operator[](std::uint32_t const &i) const
    {
        return x[i];
    }

    inline void fill(T const &val)
    {
        std::fill(x, x + n, val);
    }

    inline void zeros()
    {
        std::fill(x, x + n, 0.00);
    }
};
