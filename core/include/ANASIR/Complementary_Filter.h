/*!*******************************************
 * @file CF.h
 * @brief This file contains the definitions of    * functions used in Complementary Filter
 * ANASIR's implementation
 *********************************************/
#include <cstdint>
#include <cstring>
namespace anasir
{

  template <std::uint32_t n>
  class CF
  {
  private:
    double x[n]; /**<States Vector ; n * 1>*/
    double z[n]; /**<Measurement Vector ; n * 1>*/
    double ALPHA;

  public:
    void setCF(double (&x)[n], double (&z)[n], double ALPHA)
    {
      std::memset(this->x, x, sizeof(x));
      std::memset(this->z, z, sizeof(z));
      ALPHA = ALPHA;
    }
    void Complementary_Filter()
    {

      for (std::uint32_t i = 0; i < (n); i++)
      {
        x[i] = ALPHA * x[i] + (1 - ALPHA) * z[i];
      }
    }
  };
} // namespace anasir