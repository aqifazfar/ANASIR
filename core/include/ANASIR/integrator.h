#include <fstream>
#include <cmath>
#include <cstdint>
#include <functional>

template <std::uint32_t NStates>
void DOPRI(double (&x)[NStates], std::function<void(double (&)[NStates], double (&)[NStates])> odeFunc, double const &dt)
{

    double k[NStates][6] = {};
    double y[NStates] = {};
    double xtemp[NStates] = {};
    double x4[NStates] = {};
    double x5[NStates] = {};
    double error;

    while (true)
    {
        odeFunc(y, x);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][0] = dt * y[i];

            xtemp[i] = x[i] + (1.0 / 5.0 * k[i][0]);
        }

        odeFunc(y, xtemp);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][1] = dt * y[i];

            xtemp[i] = x[i] + (3.0 / 40.0 * k[i][0] + 9.0 / 40.0 * k[i][1]);
        }

        odeFunc(y, xtemp);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][2] = dt * y[i];

            xtemp[i] = x[i] + (44.0 / 45.0 * k[i][0] - 56.0 / 15.0 * k[i][1] + 32.0 / 9.0 * k[i][2]);
        }

        odeFunc(y, xtemp);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][3] = dt * y[i];

            xtemp[i] = x[i] + (19372.0 / 6561.0 * k[i][0] - 25360.0 / 2187.0 * k[i][1] + 64448.0 / 6561.0 * k[i][2] - 212.0 / 729.0 * k[i][3]);
        }

        odeFunc(y, xtemp);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][4] = dt * y[i];

            xtemp[i] = x[i] + (9017.0 / 3168.0 * k[i][0] - 355.0 / 33.0 * k[i][1] + 46732.0 / 5247.0 * k[i][2] - 49.0 / 176.0 * k[i][3] - 5103.0 / 18656.0 * k[i][4]);
        }

        odeFunc(y, xtemp);

        error = 0.0;
        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][5] = dt * y[i];

            x4[i] = x[i] + 5179.0 / 57600.0 * k[i][0] + 0.0 * k[i][1] + 7571.0 / 16695.0 * k[i][2] + 393.0 / 640.0 * k[i][3] - 92097.0 / 339200.0 * k[i][4] + 187.0 / 2100.0 * k[i][5];

            x5[i] = x[i] + 35.0 / 384.0 * k[i][0] + 500.0 / 1113.0 * k[i][2] + 125.0 / 192.0 * k[i][3] - 2187.0 / 6784.0 * k[i][4] + 11.0 / 84.0 * k[i][5];

            error = std::max(error, std::abs(x5[i] - x4[i]));
        }

        if (error <= 1e-6)
        {
            for (std::uint32_t i = 0; i < NStates; i++)
            {
                x[i] = x5[i];
            }
            break;
        }

        dt = dt * 0.9 * std::pow((1e-6 / error), 0.2);

        if (dt < 1e-6)
        {
            throw std::runtime_error("cannot proceed dt too small");
        }
    }

    dt = std::clamp(0.9 * dt * std::pow(1e-6 / error, 0.2), 1e-6, 1.0);
}

template <std::uint32_t NStates>
void RK4(double (&x)[NStates], std::function<void(double (&)[NStates], double (&)[NStates])> odeFunc, double const dt)
{

    double k[NStates][6] = {};
    double y[NStates] = {};
    double xtemp[NStates] = {};

    odeFunc(y, x);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][0] = dt * y[i];

        xtemp[i] = x[i] + (0.5 * k[i][0]);
    }

    odeFunc(y, xtemp);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][1] = dt * y[i];

        xtemp[i] = x[i] + (0.5 * k[i][1]);
    }

    odeFunc(y, xtemp);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][2] = dt * y[i];

        xtemp[i] = x[i] + (k[i][2]);
    }

    odeFunc(y, xtemp);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][3] = dt * y[i];

        x[i] = x[i] + 1.0 / 6.0 * (k[i][0] + 2.0 * k[i][1] + 2.0 * k[i][2] + k[i][3]);
    }
}

template <std::uint32_t NStates>
void Euler_Method(double (&x)[NStates], std::function<void(double (&)[NStates], double (&)[NStates])> odeFunc, double const dt)
{

    double xtemp[NStates] = {};

    odeFunc(xtemp, x);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        x[i] = x[i] + dt * xtemp[i];
    }
}
