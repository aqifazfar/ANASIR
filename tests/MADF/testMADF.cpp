
#include <iostream>
#include <fstream>
#include "ANASIR/Madgwick_Filter"
#include "../csv.h"

int main()
{
    anasir::Madgwick_Filter<float> filter(0.01f, 0.05f, 0.01f, 0.02f, 0.1f);

    std::ofstream output("MADF_output.csv");

    output << "roll,pitch,yaw\n";

    io::CSVReader<9> in("imu_data.csv");
    std::cout << "Read\n";
    in.read_header(io::ignore_extra_column, "ax", "ay", "az", "gx", "gy", "gz", "hx", "hy", "hz");
    float accX = 0.0f, accY = 0.0f, accZ = 0.0f;
    float wX = 0.0f, wY = 0.0f, wZ = 0.0f;
    float hX = 0.0f, hY = 0.0f, hZ = 0.0f;

    while (in.read_row(accX, accY, accZ, wX, wY, wZ, hX, hY, hZ))
    {
        float acc[3] = {accX, accY, accZ};
        float mag[3]{hX, hY, hZ};
        float gyro[3] = {wX * M_PI / 180.0, wY * M_PI / 180.0, wZ * M_PI / 180.0};

        anasir_states<float, 3> states = anasir::quat_To_Euler((filter.Madgwick_Filter_AHRS(acc, gyro, mag)));

        output << states[0] * 180.0 / M_PI << "," << states[1] * 180.0 / M_PI << "," << states[2] * 180.0 / M_PI << "\n";
    }

    return 0;
}