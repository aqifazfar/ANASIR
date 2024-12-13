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
  template <typename T, std::uint32_t NStates, std::uint32_t NSensor>
  class EKF
  {
  private:
    T x[NStates];           /**<States Vector ; nx * 1>*/
    T fx[NStates];          /**<States-Transition Function ; nx * 1>*/
    T F[NStates * NStates]; /**<Jacobian of States-Transition Function; nx * nx*/

    T P[NStates * NStates]; /**<Noise Covariance Matrix ; nx * nx>*/
    T Q[NStates * NStates]; /**<Process Noise Matrix ; nx * nx>*/
    T R[NSensor * NSensor]; /**<Measurement Noise Matrix ; nz * nz>*/
    T z[NSensor];           /**<Measurement Vector ; nz * 1>*/

    T hx[NSensor];          /**<Measurement Function ; nz * 1>*/
    T H[NSensor * NStates]; /**<Jacobian of Measurement Function ; nz * nx */

  public:
    void Set_State_Model(T (&xIn)[NStates], T (&fxIn)[NStates], T (&FIn)[NStates * NStates])
    {
      x = xIn;
      fx = fxIn;
      F = FIn;
    }

    void Set_Measurement_Model(T (&zIn)[NSensor], T (&hxIn)[NSensor], T (&HIn)[NSensor * NStates])
    {
      z = zIn;
      hx = hxIn;
      H = HIn;
    }

    void Set_Covariance(T PInit, T QInit, T RInit)
    {
      Matrix::Diagonal<T, NStates>(P, PInit);
      Matrix::Diagonal<T, NStates>(Q, QInit);
      Matrix::Diagonal<T, NSensor>(R, RInit);
    }
    void
    Extended_Kalman_Filter()
    {
      T tmp0[NSensor * NSensor] = {};
      T tmp1[NSensor * NSensor] = {};
      T tmp2[NStates] = {};
      T tmp3[NStates * NStates] = {};
      T tmp4[NStates * NStates] = {};
      T tmp5[NStates * NSensor] = {};
      static T K[NStates * NSensor];

      // PREDICTION
      // States Prediction : Xn+1 = f(Xn)

      std::memcpy(x, fx, NStates * sizeof(T));

      // Covariance Matrix Prediction : Pn+1 = FPnFt + Q
      Matrix::Mul<T, NStates, NStates, NStates>(P, F, P);

      Matrix::Transpose<T, NStates, NStates>(F);

      Matrix::Mul<T, NStates, NStates, NStates>(P, P, F);

      Matrix::Add<T, NStates, NStates>(P, P, Q);

      // UPDATE
      // Computes Kalman Gain : K = PHt *(HPHt + R)^-1

      // H ->Ht
      Matrix::Transpose<T, NSensor, NStates>(H);

      // K = PHt
      Matrix::Mul<T, NStates, NStates, NSensor>(K, P, H);

      // (HPHt + R)^-1

      // Ht -> H
      Matrix::Transpose<T, NSensor, NStates>(H);

      // H * PHt = H * K
      Matrix::Mul<T, NSensor, NStates, NSensor>(tmp0, H, K);

      Matrix::Add<T, NSensor, NSensor>(tmp0, tmp0, R);

      Matrix::CholeskyDecomposition<T, NSensor>(tmp0, tmp0);

      Matrix::Inv<T, NSensor>(tmp0, tmp0); // Inverse of Upper/Lower triangle

      std::memcpy(tmp1, tmp0, sizeof(tmp0));

      Matrix::Transpose<T, NSensor, NSensor>(tmp0);

      Matrix::Mul<T, NSensor, NSensor, NSensor>(tmp1, tmp1, tmp0); // Full Inverse Matrix

      Matrix::Mul<T, NStates, NSensor, NSensor>(K, K, tmp1);

      // Update States : Xn = Xn-1 + K * (z - hx)

      Matrix::Sub<T, NSensor, 1>(z, z, hx);

      Matrix::Mul<T, NStates, NSensor, 1>(tmp2, K, z);

      Matrix::Add<T, NStates, 1>(x, x, tmp2);

      // Update Covariance Matrix : P = (I - KH) * P * (I - KH)t + KRKt

      Matrix::Diagonal<T, NStates>(tmp3, 1);

      Matrix::Mul<T, NStates, NSensor, NStates>(tmp4, K, H);

      Matrix::Sub<T, NStates, NStates>(tmp3, tmp3, tmp4);

      Matrix::Mul<T, NStates, NStates, NStates>(P, tmp3, P);

      // (I-KH) -> (I-KH)t
      Matrix::Transpose<T, NStates, NStates>(tmp3);

      Matrix::Mul<T, NStates, NStates, NStates>(P, P, tmp3);

      // KR
      Matrix::Mul<T, NStates, NSensor, NSensor>(tmp5, K, R);

      // K->Kt
      Matrix::Transpose<T, NStates, NSensor>(K);

      Matrix::Mul<T, NStates, NSensor, NStates>(tmp3, tmp5, K);

      // Kt->K
      Matrix::Transpose<T, NStates, NSensor>(K);

      // Update Covariance Matrix

      Matrix::Add<T, NStates, NStates>(P, P, tmp3);
    }
  };
} // namespace anasir