#include <cmath>
#include <cstdint>
#include <numeric>

#include "anasir_types.h"
namespace anasir
{

    template <typename T, std::uint32_t n>
    void Normalize_Vector(T (&x)[n])
    {
        T sum = std::accumulate(x, x + n, 0.0, [](T accumulate, T element)
                                { return accumulate + element * element; });

        if (sum == 0.0)
        {
            // Skip normalization
            return;
        }

        T magnitude = std::sqrt(sum);

        for (std::uint32_t i = 0; i < n; i++)
        {
            x[i] = x[i] / magnitude;
        }
    }

    template <typename T>
    anasir_states<T, 3> quat_To_Euler(anasir_states<T, 4> const &q)
    {
        anasir_states<T, 3> states;

        // Roll
        states[0] = std::atan2(2.0 * (q[0] * q[1] + q[2] * q[3]), 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]));

        // test pitch
        T test = 2.0 * (q[0] * q[2] - q[3] * q[1]);

        if (std::abs(test) >= 1.0)
        {
            states[1] = -std::copysign(M_PI / 2.0, test);
        }
        else
        {
            states[1] = -std::asin(test);
        }

        // Yaw
        states[2] = -std::atan2(2.0 * (q[0] * q[3] + q[1] * q[2]), 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]));

        return states;
    }
} // namespace anasir
