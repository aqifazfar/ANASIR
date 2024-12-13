/*!*******************************************
 * @file MAHF.h
 * @brief This file contains the definitions of    * functions used in Mahony Filter
 * ANASIR's implementation
 *********************************************/

#include <cstdint>
#include <cmath>
#include <cstring>

namespace anasir
{
  template <typename T>
  class MAHF
  {
  private:
    T a[3];
    T g[3];
    T m[3];

    T q[4];
    T Kp;
    T Ki;
    T dt;

    static inline T
    sqr(T x)
    {

      return x * x;
    }

    static T Norm(T &x, T n)
    {

      T sum;

      sum = 0;

      for (std::uint8_t i = 0; i < n; i++)
      {

        sum += sqr(x[i]);
      }

      for (std::uint8_t i = 0; i < n; i++)
      {

        x[i] /= std::sqrt(sum);
      }
      return 0; // success
    }

  public:
    void Mahony_Filter_Init(T Kp, T Ki, T (&AccIn)[3], T (&GyroIn)[3], T (&MagIn)[3])
    {
      Kp = Kp;
      Ki = Ki;
      std::memset(m, MagIn, 3 * sizeof(T));
      std::memset(g, GyroIn, 3 * sizeof(T));
      std::memset(a, AccIn, 3 * sizeof(T));

      q[0] = 1.00f;
      q[1] = 0.00f;
      q[2] = 0.00f;
      q[3] = 0.00f;
    }

    T Mahony_Filter_IMU()
    {

      T vx, vy, vz, error[3];
      std::uint8_t i;

      // Compute only when accelerometer data are available
      if (a[0] == 0.0f && a[1] == 0.0f && a[2] == 0.0f)
      {

        return 1; // error
      }

      // normalized accelerometer
      Norm(a, 3);

      // Estimate direction of gravity, V
      vx = 2 * (q[1] * q[3] - q[0] * q[2]);

      vy = 2 * (q[0] * q[1] + q[2] * q[3]);

      vz = sqr(q[0]) - sqr(q[1]) - sqr(q[2]) + sqr(q[3]);

      // Calculate error using cross-product
      //  a x v
      error[0] = a[1] * vz - a[2] * vy;
      error[1] = a[2] * vx - a[0] * vz;
      error[2] = a[0] * vy - a[1] * vx;

      // Update Gyroscope measurements
      // w = w + Kp * w + Ki * w
      for (i = 0; i < 3; i++)
      {

        error[i] += error[i] * dt;

        g[i] += Kp * error[i] + Ki * error[i];
      }

      T qw, qx, qy, qz;

      qw = q[0];
      qx = q[1];
      qy = q[2];
      qz = q[3];

      q[0] += (dt * 0.5f) * (-qx * g[0] - qy * g[1] - qz * g[2]);

      q[1] += (dt * 0.5f) * (qw * g[0] + qy * g[2] - qz * g[1]);

      q[2] += (dt * 0.5f) * (qw * g[1] - qx * g[2] + qz * g[0]);

      q[3] += (dt * 0.5f) * (qw * g[2] + qx * g[1] - qy * g[0]);

      Norm(q, 4);

      return 0;
    }

    T Mahony_Filter_AHRS()
    {

      // Compute only when accelerometer data are available
      if (a[0] == 0.0f && a[1] == 0.0f && a[2] == 0.0f)
      {

        return 1; // error
      }

      // Use IMU algorithm if mx,my and mz is invalid
      if (m[0] == 0.0f && m[1] == 0.0f && m[2] == 0.0f)
      {

        Mahony_Filter_IMU();

        return 0; // success
      }

      Norm(a, 3);
      Norm(m, 3);

      T vx, vy, vz, hx, hy, bx, bz, wx, wy, wz, error[3];
      std::uint8_t i;

      // Compute reference direction of earth magnetic field, B

      hx = 2.0f * (m[0] * (0.5f - sqr(q[2]) - sqr(q[3])) + m[1] * (q[1] * q[2] - q[0] * q[3]) + m[2] * (q[1] * q[3] + q[0] * q[2]));

      hy = 2.0f * (m[0] * (q[1] * q[2]) + m[1] * (0.5f - sqr(q[1]) - sqr(q[3])) + m[2] * (q[2] * q[3]));

      bx = std::sqrt((sqr(hx) + sqr(hy)));

      bz = 2.0f * (m[0] * (q[1] * q[3] - q[0] * q[2]) + m[1] * (q[2] * q[3] + q[0] * q[1]) + m[2] * (0.5f - sqr(q[1]) - sqr(q[2])));

      // Estimate direction of gravity, V
      vx = 2.0f * (q[1] * q[3] - q[0] * q[2]);

      vy = 2.0f * (q[0] * q[1] + q[2] * q[3]);

      vz = sqr(q[0]) - sqr(q[1]) - sqr(q[2]) + sqr(q[3]);

      // Estimate direction of magnetic field , W

      wx = 2.0f * (bx * (0.5f - sqr(q[2]) - sqr(q[3])) + bz * (q[1] * q[3] - q[0] * q[2])) - m[0];

      wy = 2.0f * (bx * (q[1] * q[2] - q[0] * q[3]) + bz * (q[0] * q[1] + q[2] * q[3])) - m[1];

      wz = 2.0f * (bx * (q[0] * q[2] + q[1] * q[3]) + bz * (0.5f - sqr(q[1]) - sqr(q[2]))) - m[2];

      // Calculate error using cross-product
      //  ea = a x v and em = m  x w
      //  error = ea + em

      error[0] = (a[1] * vz - a[2] * vy) + (m[1] * wz - m[2] * wy);

      error[1] = (a[2] * vx - a[0] * vz) + (m[2] * wx - m[0] * wz);

      error[2] = (a[0] * vy - a[1] * vx) + (m[0] * wy - m[1] * wx);

      // Update Gyroscope measurements
      // w = w + Kp * w + Ki * w
      for (i = 0; i < 3; i++)
      {

        error[i] += error[i] * dt;

        g[i] += Kp * error[i] + Ki * error[i];
      }

      T qw, qx, qy, qz;

      qw = q[0];
      qx = q[1];
      qy = q[2];
      qz = q[3];

      q[0] += (dt * 0.5f) * (-qx * g[0] - qy * g[1] - qz * g[2]);

      q[1] += (dt * 0.5f) * (qw * g[0] + qy * g[2] - qz * g[1]);

      q[2] += (dt * 0.5f) * (qw * g[1] - qx * g[2] + qz * g[0]);

      q[3] += (dt * 0.5f) * (qw * g[2] + qx * g[1] - qy * g[0]);

      Norm(q, 4);

      return 0; // success
    }
  };
} // namespace anasir
