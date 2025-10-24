/*!*******************************************
 * @file MADF.h
 * @brief
 * This file contains the definitions of                * functions used in Madgwick Filter ANASIR's
 * implementation
 *********************************************/

#include <cstdint>
#include <cmath>
#include <cstring>

namespace anasir
{
  template <typename double>
  class MADF
  {
  private:
    double a[3];
    double g[3];
    double m[3];
    double q[4];
    double Beta;
    double dt;
    static inline double sqr(double x)
    {

      return x * x;
    }

    static double Norm(double *x, double n)
    {

      double sum;

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
    void Madgwick_Filter_Init(double Kp, double Ki, double Kd, double Beta, double dt, double (&AccIn)[3], double (&GyroIn)[3], double (&MagIn)[3])
    {
      Kp = Kp;
      Ki = Ki;
      std::memset(m, MagIn, 3 * sizeof(double));
      std::memset(g, GyroIn, 3 * sizeof(double));
      std::memset(a, AccIn, 3 * sizeof(double));

      dt = dt;

      q[0] = 1.00f;
      q[1] = 0.00f;
      q[2] = 0.00f;
      q[3] = 0.00f;
    }

    double Madgwick_Filter_IMU()
    {

      // Compute only when accelerometer data   are available
      if (a[0] == 0.0f && a[1] == 0.0f && a[2] == 0.0f)
      {

        return 1; // error
      }

      // Normalized accelerometer
      Norm(a, 3);

      // Orientation increment from gyroscope
      double qg0, qg1, qg2, qg3;

      qg0 = 0.5f * (-q[1] * g[0] - q[2] * g[1] - q[3] * g[2]);

      qg1 = 0.5f * (q[0] * g[0] + q[2] * g[2] - q[3] * g[1]);

      qg2 = 0.5f * (q[0] * g[1] - q[1] * g[2] + q[3] * g[0]);

      qg3 = 0.5f * (q[0] * g[2] + q[1] * g[1] - q[2] * g[0]);

      // Gradient Descent Algorithm
      double df[4], f[3];

      // Compute f(q , a)
      f[0] = 2.0f * (q[1] * q[3] - q[0] * q[2]) - a[0];

      f[1] = 2.0f * (q[0] * q[1] + q[2] * q[3]) - a[1];

      f[2] = 2.0f * (0.5f - sqr(q[1]) - sqr(q[2])) - a[3];

      // Compute delta f
      df[0] = (-2.0f * q[2] * f[0]) + (2.0f * q[1] * f[1]);

      df[1] = (2.0f * q[3] * f[0]) + (2.0f * q[0] * f[1]) + (-4.0f * q[1] * f[2]);

      df[2] = (-2.0f * q[0] * f[0]) + (2.0f * q[3] * f[1]) + (-4.0f * q[2] * f[2]);

      df[3] = (2.0f * q[1] * f[0]) + (2.0f * q[2] * f[1]);

      Norm(df, 4);

      // Scale with Beta
      df[0] *= Beta;
      df[1] *= Beta;
      df[2] *= Beta;
      df[3] *= Beta;
      

      // Integrate rate of change of quaternion
      q[0] += (qg0 - df[0]) * dt;
      q[1] += (qg1 - df[1]) * dt;
      q[2] += (qg2 - df[2]) * dt;
      q[3] += (qg3 - df[2]) * dt;

      // Normalized quaternion
      Norm(q, 4);

      return 0; // success
    }

    double Madgwick_Filter_AHRS()
    {

      // Compute only when accelerometer data   are available
      if (a[0] == 0.0f && a[1] == 0.0f && a[2] == 0.0f)
      {

        return 1; // error
      }

      // Use IMU algorithm if mx,my and mz is invalid
      if (m[0] == 0.0f && m[1] == 0.0f && m[2] == 0.0f)
      {

        Madgwick_Filter_IMU();

        return 0; // success
      }
      // Normalized accelerometer and magnetometer
      Norm(a, 3);
      Norm(m, 3);

      // Orientation increment from gyroscope
      double qg0, qg1, qg2, qg3;

      qg0 = 0.5f * (-q[1] * g[0] - q[2] * g[1] - q[3] * g[2]);

      qg1 = 0.5f * (q[0] * g[0] + q[2] * g[2] - q[3] * g[1]);

      qg2 = 0.5f * (q[0] * g[1] - q[1] * g[2] + q[3] * g[0]);

      qg3 = 0.5f * (q[0] * g[2] + q[1] * g[1] - q[2] * g[0]);

      // Compute reference direction of earth magnetic field, B
      double hx, hy, bx, bz;

      hx = 2.0f * (m[0] * (0.5f - sqr(q[2]) - sqr(q[3])) + m[1] * (q[1] * q[2] - q[0] * q[3]) + m[2] * (q[1] * q[3] + q[0] * q[2]));

      hy = 2.0f * (m[0] * (q[1] * q[2]) + m[1] * (0.5f - sqr(q[1]) - sqr(q[3])) + m[2] * (q[2] * q[3]));

      bx = std::sqrt((sqr(hx) + sqr(hy)));

      bz = 2.0f * (m[0] * (q[1] * q[3] - q[0] * q[2]) + m[1] * (q[2] * q[3] + q[0] * q[1]) + m[2] * (0.5f - sqr(q[1]) - sqr(q[2])));

      // Gradient Descent Algorithm
      double df[4], fqa[3], fqm[3];

      // Compute f(q,a)
      fqa[0] = 2.0f * (q[1] * q[3] - q[0] * q[2]) - a[0];

      fqa[1] = 2.0f * (q[0] * q[1] + q[2] * q[3]) - a[1];

      fqa[2] = 2.0f * (0.5f - sqr(q[1]) - sqr(q[2])) - a[2];

      // Compute f(q,m)
      fqm[0] = 2.0f * (bx * (0.5f - sqr(q[2]) - sqr(q[3])) + bz * (q[1] * q[3] - q[0] * q[2])) - m[0];

      fqm[1] = 2.0f * (bx * (q[1] * q[2] - q[0] * q[3]) + bz * (q[0] * q[1] + q[2] * q[3])) - m[1];

      fqm[2] = 2.0f * (bx * (q[0] * q[2] + q[1] * q[3]) + bz * (0.5f - sqr(q[1]) - sqr(q[2]))) - m[2];

      // Compute delta f
      df[0] = 2.0f * (-q[2] * fqa[0] + q[1] * fqa[1]) + 2.0f * (-bz * q[2] * fqm[0] + (-bx * q[3] + bz * q[1]) * fqm[1] + bx * q[2] * fqm[2]);

      df[1] = 2.0f * (q[3] * fqa[0] + q[0] * fqa[1] - 2.0f * q[1] * fqa[2]) + 2.0f * (bz * q[3] * fqm[0] + (bx * q[2] + bz * q[0]) * fqm[1] + (bx * q[3] - 2.0f * bz * q[1]) * fqm[2]);

      df[2] = 2.0f * (-q[0] * fqa[0] + q[3] * fqa[1] - 2.0f * q[2] * fqa[2]) + 2.0f * ((-2.0f * bx * q[2] - bz * q[0]) * fqm[0] + (bx * q[1] + bz * q[3]) * fqm[1] + (bx * q[0] - 2.0f * bz * q[2]) * fqm[2]);

      df[3] = 2.0f * (q[1] * fqa[0] + q[2] * fqa[1]) + 2.0f * ((-2.0f * bx * q[3] + bz * q[1]) * fqm[0] + (-bx * q[0] + bz * q[2]) * fqm[1] + bx * q[1] * fqm[2]);

      Norm(df, 4);

      // Scale with Beta
      df[0] *= Beta;
      df[1] *= Beta;
      df[2] *= Beta;
      df[3] *= Beta;

      // Integrate rate of change of quaternion
      q[0] += (qg0 - df[0]) * dt;
      q[1] += (qg1 - df[1]) * dt;
      q[2] += (qg2 - df[2]) * dt;
      q[3] += (qg3 - df[2]) * dt;

      // Normalized quaternion
      Norm(q, 4);

      return 0; // success
    }
  };
} // namespace anasir