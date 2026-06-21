#pragma once

#define EARTH_SEMI_MAJOR 6378.137     // earth semi major axis length, km
#define EARTH_SEMI_MINOR 6356.7523142 // earth semi minor axis length, km
#define RE 6371.2                     // geomagnetic reference radius (which is close to the mean Earth radius) km,a
#define MAX_ORDER 12                  // max order of Associate Legendre Functions,n
#define MAX_DEGREE 12                 // max degree of Associate Legendre Functions,m

#define A2 (EARTH_SEMI_MAJOR * EARTH_SEMI_MAJOR)
#define B2 (EARTH_SEMI_MINOR * EARTH_SEMI_MINOR)
#define C2 (A2 - B2)
#define A4 (A2 * A2)
#define B4 (B2 * B2)

class Geo_Mag
{
private:
    double horizontal_intensity = 0.0; // H in NOAA Doc
    double total_intensity = 0.0;      // F in NOAA Doc
    double inclination = 0.0;          // I in NOAA Doc
    double declination = 0.0;          // D in NOAA Doc
    double b[3] = {0.0};               // Magnetic vector in 3D NED

    double c[(MAX_ORDER + 1)][(MAX_ORDER + 1)] = {0.0};     //  Gauss Coefficients g adn h in upper triangular form
    double c_dot[(MAX_ORDER + 1)][(MAX_ORDER + 1)] = {0.0}; //  Gauss Coefficients g adn h in upper triangular form
    double tc[(MAX_ORDER + 1)][(MAX_ORDER + 1)] = {0.0};    // time adjust Gauss Coefficients
    double snorm[169] = {0.0};                              // schmdit n0rmalization factor
    double k[(MAX_ORDER + 1)][(MAX_ORDER + 1)] = {0.0};     // coefficients for legendre polynomials (upper triangular form)
    double fn[(MAX_ORDER + 1)] = {0.0};
    double fm[(MAX_ORDER + 1)] = {0.0};
    double epoch = 0.0;

public:
    Geo_Mag();
    ~Geo_Mag();
    void E0000(double const &latitude, double const &longitude, double const &altitude, double const &year);
    double *Get_Magnetic_Vector();
    double Get_Horizontal_Intensity();
    double Get_Total_Intensity();
    double Get_Inclination();
    double Get_Declination();
};
