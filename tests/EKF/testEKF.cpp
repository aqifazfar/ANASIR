#include <iostream>
#include <vector>
#include <string>
#include "ANASIR/Extended_Kalman_Filter"
#include "ANASIR/helpers.h"
#include "../csv.h"
int main()
{

    io::CSVReader<13> in("uavdata.csv");
    std::cout << "Read\n";
    in.read_header(io::ignore_extra_column, "TimeFromStart (s)", "AccelX (m/s^2)", "AccelY (m/s^2)", "AccelZ (m/s^2)", "AngleRateX (rad/s)", "AngleRateY (rad/s)", "AngleRateZ (rad/s)", "GNSSPosAlt (m)", "GNSSPosLat (deg)", "GNSSPosLon (deg)", "MagFieldX (G)", "MagFieldY (G)", "MagFieldZ (G)");
    float accX = 0.0f, accY = 0.0f, accZ = 0.0f;
    float wX = 0.0f, wY = 0.0f, wZ = 0.0f;
    float alt = 0.0f, lat = 0.0f, lon = 0.0f;
    float hX = 0.0f, hY = 0.0f, hZ = 0.0f;
    float tfs = 0.0f;

    float dt = 0.02f;
    float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};

    EKF<float, 4, 3> ekf(q, 0.01f, 0.01f, 0.01f, dt);

    while (in.read_row(tfs, accX, accY, accZ, wX, wY, wZ, alt, lat, lon, hX, hY, hZ))
    {
        float fx[4] = {dt * 0.5f * (-wX * q[1] - wY * q[2] - wZ * q[3]), dt * 0.5f * (wX * q[0] + wZ * q[2] - wY * q[3]), dt * 0.5f * (wY * q[0] - wZ * q[1] + wX * q[3]), dt * 0.5f * (wZ * q[0] + wY * q[1] - wX * q[2])};
        float F[16] = {0.00f, -0.5f * dt * wX, -0.5f * dt * wY, -0.5f * dt * wZ, 0.5f * dt * wX, 0, 0.5f * dt * wZ, -0.5f * dt * wY, 0.5f * dt * wY, -0.5f * dt * wZ, 0, 0.5f * dt * wX, 0.5f * dt * wZ, 0.5f * dt * wY, -0.5f * dt * wX, 0};
        float hx[3] = {2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 1 - 2 * (q[1] * q[1] + q[2] * q[2])};
        float H[12] = {-2 * q[2], 2 * 2 * q[3], -2 * q[0], 2 * q[1], 2 * q[1], 2 * q[0], 2 * q[3], 2 * q[2], 0, -4 * q[1], -4 * q[2], 0};
        float z[3] = {accX, accY, accZ};

        anasir::Normalize_Vector<float, 4>(fx);
        anasir::Normalize_Vector<float, 3>(hx);

        anasir_states<float, 4> states = ekf.EKF_Process(fx, F, hx, H, z);

        anasir::Normalize_Vector<float, 4>(states.x);

        anasir_states<float, 3> angle = anasir::quat_To_Euler(states);

        std::cout << "roll = " << angle[0] * 180.0 / M_PI << " pitch = " << angle[1] * 180.0 / M_PI << " yaw = " << angle[2] * 180.0 / M_PI;

        std::cout << std::endl;
    }

    return 0;
}