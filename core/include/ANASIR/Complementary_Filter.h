/*!*******************************************
 * @file CF.h
 * @brief This file contains the definitions of    * functions used in Complementary Filter
 * ANASIR's implementation
 *********************************************/
#include <cstdint>
namespace anasir
{

  template <std::uint32_t n>
  class CF
  {
  private:
    T x[n]; /**<States Vector ; n * 1>*/
    T z[n]; /**<Measurement Vector ; n * 1>*/
    T ALPHA;

  public:
    void setCF(T (&x)[n], T (&z)[n], T ALPHA)
    {
      std::memset(x, x, sizeof(x));
      std::memset(z, z, sizeof(z));
      ALPHA = ALPHA;
    }
    void Complementary_Filter()
    {

      std::uint32_t i;

      for (i = 0; i < (n); i++)
      {
        x[i] = ALPHA * x[i] + (1 - ALPHA) * z[i];
      }
    }
  };
} // namespace anasir