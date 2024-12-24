#include <iostream>
#include <vector>
#include <string>
#include "ANASIR/Extended_Kalman_Filter.h"
#include "csv.h"
int main()
{

    io::CSVReader<13> in("uavdata.csv");
    std::cout << "Read\n";
    in.read_header(io::ignore_extra_column, "TimeFromStart (s)", "AccelX (m/s^2)", "AccelY (m/s^2)", "AccelZ (m/s^2)", "AngleRateX (rad/s)", "AngleRateY (rad/s)", "AngleRateZ (rad/s)", "GNSSPosAlt (m)", "GNSSPosLat (deg)", "GNSSPosLon (deg)", "MagFieldX (G)", "MagFieldY (G)", "MagFieldZ (G)");
    float accX, accY, accZ;
    float wX, wY, wZ;
    float alt, lat, lon;
    float hX, hY, hZ;
    float tfs;

    anasir::EKF<4, 6> ekfilter;
    float dt = 0.00f;
    float q[4] = {1, 0, 0, 0};
    float fx[4] = {dt * 0.5f * (-wX * q[1] - wY * q[2] - wZ * q[3]), dt * 0.5f * (wX * q[0] + wZ * q[2] - wY * q[3]), dt * 0.5f * (wY * q[0] - wZ * q[1] + wX * q[3]), dt * 0.5f * (wZ * q[0] + wY * q[1] - wX * q[2])};
    float F[16] = {0, -0.5f * dt * wX, -0.5f * dt * wY, -0.5f * dt * wZ, 0.5f * dt * wX, 0, 0.5f * dt * wZ, -0.5f * dt * wY, 0.5f * dt * wY, -0.5f * dt * wZ, 0, 0.5f * dt * wX, 0.5f * dt * wZ, 0.5f * dt * wY, -0.5f * dt * wX, 0};

    float hx[3] = {2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 1 - 2 * (q[1] * q[1] + q[2] * q[2])};
    float H[12] = {-2 * q[2], 2 * 2 * q[3], -2 * q[0], 2 * q[1], 2 * q[1], 2 * q[0], 2 * q[3], 2 * q[2], 0, -4 * q[1], -4 * q[2], 0};
    ekfilter.Set_Covariance(0.1, 0.1, 0.1);

    while (in.read_row(tfs, accX, accY, accZ, wX, wY, wZ, alt, lat, lon, hX, hY, hZ))
    {
        dt = tfs - dt;
    }

    return 0;
}