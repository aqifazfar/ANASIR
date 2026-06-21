#include <cstdint>
#include <cmath>
#include <cmath>
#include <ctime>
#include <functional>

#define EARTH_SEMI_MAJOR 6378137.0 // earth semi major axis length, m
#define ECCENTRIC_2 0.006694379    // square of eccentric of earth

void Quaternion_Dynamics(double (&dq)[4], double (&q)[4], double const (&w)[3])
{
    dq[0] = 0.5 * (w[0] * q[0] - w[1] * q[1] - w[2] * q[2]);
    dq[1] = 0.5 * (w[0] * q[1] + w[1] * q[0] - w[2] * q[3]);
    dq[2] = 0.5 * (w[0] * q[2] + w[1] * q[3] + w[2] * q[0]);
    dq[3] = 0.5 * (w[0] * q[3] - w[1] * q[2] + w[2] * q[1]);

    double magnitude = sqrt(dq[0] * dq[0] + dq[1] * dq[1] + dq[2] * dq[2] + dq[3] * dq[3]);

    dq[0] /= magnitude;
    dq[1] /= magnitude;
    dq[2] /= magnitude;
    dq[3] /= magnitude;
}

template <std::uint32_t n>
void Normalize_Vector(double (&x)[n])
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += x[i] * x[i];
    }

    if (sum == 0.0)
    {
        // Skip normalization
        return;
    }

    double magnitude = std::sqrt(sum);

    for (std::uint32_t i = 0; i < n; i++)
    {
        x[i] = x[i] / magnitude;
    }
}

void Geodetic_to_ECEF(double (&ecef_coord)[3], double const &latitude, double const &longitude, double const &altitude)
{

    // e^2 = f(2-f)  = (1.0 / 298.257223563) * (2.0 - 1.0 / 298.257223563) = 6.6943799e-3 = 0.006694379
    double Rc = EARTH_SEMI_MAJOR / std::sqrt(1.0 - ECCENTRIC_2 * std::sin(longitude) * std::sin(longitude));

    double coslat = std::cos(latitude);

    ecef_coord[0] = (Rc + altitude) * coslat * std::cos(longitude);
    ecef_coord[1] = (Rc + altitude) * coslat * std::sin(longitude);

    // 1.0 - eccentric^2 = 0.99330562
    ecef_coord[2] = (Rc * 0.99330562 + altitude) * std::sin(latitude);
}

void ECEF_to_NED(double (&ned_coord)[3], double const (&ecef_coord_init)[3], double const &latitude_init, double const &longitude_init, double const (&ecef_coord)[3])
{
    double dx = ecef_coord[0] - ecef_coord_init[0];
    double dy = ecef_coord[1] - ecef_coord_init[1];
    double dz = ecef_coord[2] - ecef_coord_init[2];

    double sinlat = std::sin(latitude_init);
    double coslat = std::cos(latitude_init);
    double sinlong = std::sin(longitude_init);
    double coslong = std::cos(longitude_init);

    ned_coord[0] = -dx * sinlat * coslong - dy * sinlat * sinlong + dz * coslat;
    ned_coord[1] = -dx * sinlong + dy * coslong;
    ned_coord[2] = -dx * coslat * coslong - dy * coslat * sinlong - dz * sinlat;
}

void Geodetic_to_NED(double (&ned_coord)[3], double const &latitude_init, double const &longitude_init, double const &latitude, double const &longitude, double const &altitude)
{

    double temp1[3] = {0.0}; // init
    double temp2[3] = {0.0}; // next

    Geodetic_to_ECEF(temp1, latitude_init, longitude_init, altitude);
    Geodetic_to_ECEF(temp2, latitude, longitude, altitude);
    ECEF_to_NED(ned_coord, temp1, latitude_init, longitude_init, temp2);
}

//////////////////////////////////////////////////////////// INTEGRATORS ///////////////////////////////////////////////

enum Integrator_Type
{
    EULER = 1,
    RK4 = 2,
    RK5 = 3,
    DOPRI = 4,
};

// The format of the ODE Functions shall be:
// y = f(x,u)
// odefunc(y,x,u)

template <std::uint32_t NStates, std::uint32_t NVariables>
void Euler_Method(double (&x)[NStates], double const (&u)[NVariables], double const &dt, std::function<void(double (&)[NStates], double (&)[NStates], double const (&)[NVariables])> odeFunc)
{

    double y[NStates] = {0};

    odeFunc(y, x, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        x[i] += dt * y[i];
    }
}

template <std::uint32_t NStates, std::uint32_t NVariables>
void RK4_Method(double (&x)[NStates], double const (&u)[NVariables], double const &dt, std::function<void(double (&)[NStates], double (&)[NStates], double const (&)[NVariables])> odeFunc)
{

    double k[NStates][4] = {0};
    double y[NStates] = {0};
    double xtemp[NStates] = {0};

    odeFunc(y, x, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][0] = dt * y[i];

        xtemp[i] = x[i] + (0.5 * k[i][0]);
    }

    odeFunc(y, xtemp, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][1] = dt * y[i];

        xtemp[i] = x[i] + (0.5 * k[i][1]);
    }

    odeFunc(y, xtemp, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][2] = dt * y[i];

        xtemp[i] = x[i] + (k[i][2]);
    }

    odeFunc(y, xtemp, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][3] = dt * y[i];

        x[i] += (1.0 / 6.0) * (k[i][0] + 2.0 * k[i][1] + 2.0 * k[i][2] + k[i][3]);
    }
}

template <std::uint32_t NStates, std::uint32_t NVariables>
void RK5_Method(double (&x)[NStates], double const (&u)[NVariables], double const &dt, std::function<void(double (&)[NStates], double (&)[NStates], double const (&)[NVariables])> odeFunc)
{

    double k[NStates][6] = {0};
    double y[NStates] = {0};
    double xtemp[NStates] = {0};

    odeFunc(y, x, u);
    // k0
    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][0] = dt * y[i];

        xtemp[i] = x[i] + ((1.0 / 3.0) * k[i][0]);
    }

    odeFunc(y, xtemp, u);
    // k1
    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][1] = dt * y[i];

        xtemp[i] = x[i] + ((1.0 / 25.0) * (4.0 * k[i][0] + 6.0 * k[i][1]));
    }
    // k2
    odeFunc(y, xtemp, u);

    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][2] = dt * y[i];

        xtemp[i] = x[i] + ((1.0 / 4.0) * (k[i][0] - 12.0 * k[i][1] + 15.0 * k[i][2]));
    }

    odeFunc(y, xtemp, u);
    // k3
    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][3] = dt * y[i];

        xtemp[i] = x[i] + ((1.0 / 81.0) * (6.0 * k[i][0] + 90.0 * k[i][1] - 50.0 * k[i][2] + 8.0 * k[i][3]));
    }

    odeFunc(y, xtemp, u);
    // k4
    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][4] = dt * y[i];

        xtemp[i] = x[i] + ((1.0 / 75.0) * (6.0 * k[i][0] + 36.0 * k[i][1] + 10.0 * k[i][2] + 8.0 * k[i][3]));
    }

    odeFunc(y, xtemp, u);
    // k5 and update
    for (std::uint32_t i = 0; i < NStates; i++)
    {
        k[i][5] = dt * y[i];

        x[i] += (1.0 / 192.0) * (23.0 * k[i][0] + 125.0 * k[i][2] - 81.0 * k[i][4] + 125 * k[i][5]);
    }
}

template <std::uint32_t NStates, std::uint32_t NVariables>
void DOPRI_Method(double (&x)[NStates], double const (&u)[NVariables], double &dt, double const &tol, int max_iter, std::function<void(double (&)[NStates], double (&)[NStates], double const (&)[NVariables])> odeFunc)
{

    double k[NStates][6] = {0};
    double y[NStates] = {0};
    double xtemp[NStates] = {0};
    double x4[NStates] = {0};
    double x5[NStates] = {0};
    double error = 0.0;

    double _dt = dt;

    int count = 0;

    while (true)
    {
        count++;
        odeFunc(y, x, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][0] = _dt * y[i];

            xtemp[i] = x[i] + (1.0 / 5.0 * k[i][0]);
        }

        odeFunc(y, xtemp, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][1] = _dt * y[i];

            xtemp[i] = x[i] + (3.0 / 40.0 * k[i][0] + 9.0 / 40.0 * k[i][1]);
        }

        odeFunc(y, xtemp, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][2] = _dt * y[i];

            xtemp[i] = x[i] + (44.0 / 45.0 * k[i][0] - 56.0 / 15.0 * k[i][1] + 32.0 / 9.0 * k[i][2]);
        }

        odeFunc(y, xtemp, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][3] = _dt * y[i];

            xtemp[i] = x[i] + (19372.0 / 6561.0 * k[i][0] - 25360.0 / 2187.0 * k[i][1] + 64448.0 / 6561.0 * k[i][2] - 212.0 / 729.0 * k[i][3]);
        }

        odeFunc(y, xtemp, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][4] = _dt * y[i];

            xtemp[i] = x[i] + (9017.0 / 3168.0 * k[i][0] - 355.0 / 33.0 * k[i][1] + 46732.0 / 5247.0 * k[i][2] - 49.0 / 176.0 * k[i][3] - 5103.0 / 18656.0 * k[i][4]);
        }

        odeFunc(y, xtemp, u);

        for (std::uint32_t i = 0; i < NStates; i++)
        {
            k[i][5] = _dt * y[i];

            x4[i] = x[i] + 5179.0 / 57600.0 * k[i][0] + 0.0 * k[i][1] + 7571.0 / 16695.0 * k[i][2] + 393.0 / 640.0 * k[i][3] - 92097.0 / 339200.0 * k[i][4] + 187.0 / 2100.0 * k[i][5];

            x5[i] = x[i] + 35.0 / 384.0 * k[i][0] + 500.0 / 1113.0 * k[i][2] + 125.0 / 192.0 * k[i][3] - 2187.0 / 6784.0 * k[i][4] + 11.0 / 84.0 * k[i][5];

            error = std::max(error, std::abs(x5[i] - x4[i]) / (1e-6 + tol * std::max(std::abs(x4[i]), std::abs(x5[i]))));
        }

        if ((error <= tol) || (count == max_iter))
        {
            for (std::uint32_t i = 0; i < NStates; i++)
            {
                x[i] = x5[i];
            }
            break;
        }

        _dt = _dt * 0.9 * std::pow((tol / error), 0.2);

        if (_dt < 1e-6)
        {
            // Stop
            return;
        }
    }

    dt = std::clamp(_dt, 1e-6, 1.0);
}
