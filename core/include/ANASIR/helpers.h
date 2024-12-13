#include <cmath>

// void anasir::Quat_To_Euler(double roll, double pitch, double yaw, double *q)
// {

//     double qwqy = q[0] * q[2];
//     double qxqz = q[1] * q[3];

//     double test = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
//     double unit = q[1] * q[2] + q[3] * q[0];

//     if (test > (0.499f * unit))
//     {

//         roll = 0.00f;
//         pitch = acos(0.00f);
//         yaw = 2 * std::atan2(q[1], q[0]);
//     }
//     else if (test < (0.499f * unit))
//     {
//         roll = 0.00f;
//         pitch = -acos(0.00f);
//         yaw = -2 * std::atan2(q[1], q[0]);
//     }
//     else
//     {
//         roll = std::atan2(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (q[1] * q[1] + q[2] * q[2]));
//         pitch = std::asin(2 * (qwqy - qxqz));
//         yaw = std::atan2(2 * (q[0] * q[3] + q[1] * q[2]), 1 - 2 * (q[2] * q[2] + q[3] * q[3]));
//     }
// }