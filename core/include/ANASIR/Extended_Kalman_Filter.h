/*!*******************************************
 * @file Extended_Kalman_Filter.h
 * @brief
 * This file contains the declarations of                * functions used in EKF ANASIR's
 * implementation
 *********************************************/
#include <cstdint>
#include <cstring>

#include "mat_math.h"
/**
 * @brief holds ekf parameters
 *
 */

namespace anasir
{
  template <std::uint32_t NStates, std::uint32_t NSensor>
  class EKF
  {
  private:
    double x[NStates];           /**<States Vector ; nx * 1>*/
    double fx[NStates];          /**<States-Transition Function ; nx * 1>*/
    double F[NStates * NStates]; /**<Jacobian of States-Transition Function; nx * nx*/

    double P[NStates * NStates]; /**<Noise Covariance Matrix ; nx * nx>*/
    double Q[NStates * NStates]; /**<Process Noise Matrix ; nx * nx>*/
    double R[NSensor * NSensor]; /**<Measurement Noise Matrix ; nz * nz>*/
    double z[NSensor];           /**<Measurement Vector ; nz * 1>*/

    double hx[NSensor];          /**<Measurement Function ; nz * 1>*/
    double H[NSensor * NStates]; /**<Jacobian of Measurement Function ; nz * nx */

  public:
    void Set_State_Model(double (&xIn)[NStates], double (&fxIn)[NStates], double (&FIn)[NStates * NStates])
    {
      x = xIn;
      fx = fxIn;
      F = FIn;
    }

    void Set_Measurement_Model(double (&zIn)[NSensor], double (&hxIn)[NSensor], double (&HIn)[NSensor * NStates])
    {
      z = zIn;
      hx = hxIn;
      H = HIn;
    }

    void Set_Covariance(double PInit, double QInit, double RInit)
    {
      Matrix::Diagonal<NStates>(P, PInit);
      Matrix::Diagonal<NStates>(Q, QInit);
      Matrix::Diagonal<NSensor>(R, RInit);
    }
    void
    Extended_Kalman_Filter()
    {
      double tmp0[NSensor * NSensor] = {};
      double tmp1[NSensor * NSensor] = {};
      double tmp2[NStates] = {};
      double tmp3[NStates * NStates] = {};
      double tmp4[NStates * NStates] = {};
      double tmp5[NStates * NSensor] = {};
      static double K[NStates * NSensor];

      // PREDICTION
      // States Prediction : Xn+1 = f(Xn)

      std::memcpy(x, fx, NStates * sizeof(double));

      // Covariance Matrix Prediction : Pn+1 = FPnFt + Q
      Matrix::Mul<NStates, NStates, NStates>(P, F, P);

      Matrix::Transpose<NStates, NStates>(F);

      Matrix::Mul<NStates, NStates, NStates>(P, P, F);

      Matrix::Add<NStates, NStates>(P, P, Q);

      // UPDATE
      // Computes Kalman Gain : K = PHt *(HPHt + R)^-1

      // H ->Ht
      Matrix::Transpose<NSensor, NStates>(H);

      // K = PHt
      Matrix::Mul<NStates, NStates, NSensor>(K, P, H);

      // (HPHt + R)^-1

      // Ht -> H
      Matrix::Transpose<NSensor, NStates>(H);

      // H * PHt = H * K
      Matrix::Mul<NSensor, NStates, NSensor>(tmp0, H, K);

      Matrix::Add<NSensor, NSensor>(tmp0, tmp0, R);

      Matrix::CholeskyDecomposition<NSensor>(tmp0, tmp0);

      Matrix::Inv<NSensor>(tmp0, tmp0); // Inverse of Upper/Lower triangle

      std::memcpy(tmp1, tmp0, sizeof(tmp0));

      Matrix::Transpose<NSensor, NSensor>(tmp0);

      Matrix::Mul<NSensor, NSensor, NSensor>(tmp1, tmp1, tmp0); // Full Inverse Matrix

      Matrix::Mul<NStates, NSensor, NSensor>(K, K, tmp1);

      // Update States : Xn = Xn-1 + K * (z - hx)

      Matrix::Sub<NSensor, 1>(z, z, hx);

      Matrix::Mul<NStates, NSensor, 1>(tmp2, K, z);

      Matrix::Add<NStates, 1>(x, x, tmp2);

      // Update Covariance Matrix : P = (I - KH) * P * (I - KH)t + KRKt

      Matrix::Diagonal<NStates>(tmp3, 1);

      Matrix::Mul<NStates, NSensor, NStates>(tmp4, K, H);

      Matrix::Sub<NStates, NStates>(tmp3, tmp3, tmp4);

      Matrix::Mul<NStates, NStates, NStates>(P, tmp3, P);

      // (I-KH) -> (I-KH)t
      Matrix::Transpose<NStates, NStates>(tmp3);

      Matrix::Mul<NStates, NStates, NStates>(P, P, tmp3);

      // KR
      Matrix::Mul<NStates, NSensor, NSensor>(tmp5, K, R);

      // K->Kt
      Matrix::Transpose<NStates, NSensor>(K);

      Matrix::Mul<NStates, NSensor, NStates>(tmp3, tmp5, K);

      // Kt->K
      Matrix::Transpose<NStates, NSensor>(K);

      // Update Covariance Matrix

      Matrix::Add<NStates, NStates>(P, P, tmp3);
    }
  };
} // namespace anasir