#include <cstdio>
#include <cstdint>
#include <cmath>
#include <ctime>

#include "Geo_Mag.hpp"

Geo_Mag::Geo_Mag()
{

    FILE *fp = std::fopen("WMM.COF", "r");
    char temp[81];
    int m = 0, n = 1;
    char model[20];

    if (fp == NULL)
    {
        std::printf("ERROR: WMM.COF not found.\n\n");
        return;
    }

    std::fgets(temp, 80, fp);

    if (sscanf(temp, "%lf %19s", &epoch, model) < 2)
    {
        std::printf("Invalid header in model file WMM.COF\n\n");
        std::fclose(fp);
        return;
    }

    while (std::fgets(temp, 80, fp) != NULL)
    {
        double g_nm = 0.0, gdot_nm = 0.0, h_nm = 0.0, hdot_nm = 0.0;

        std::sscanf(temp, "%d %d %lf %lf %lf %lf", &n, &m, &g_nm, &h_nm, &gdot_nm, &hdot_nm);

        if (n > MAX_ORDER)
        {
            break;
        }

        c[m][n] = g_nm;
        c_dot[m][n] = gdot_nm;

        if (m != 0)
        {
            c[n][m - 1] = h_nm;
            c_dot[n][m - 1] = hdot_nm;
        }
    }

    std::fclose(fp);

    /*Unnormalized the Schmidt Normalized Gauss Coefficients*/
    snorm[0] = 1.0;
    fm[0] = 0.0;
    for (n = 1; n <= MAX_ORDER; n++)
    {
        snorm[n] = snorm[n - 1] * static_cast<double>(2.0 * n - 1.0) / static_cast<double>(n);
        double j = 2.0;
        for (m = 0; m <= n; m++)
        {
            k[m][n] = static_cast<double>((n - 1.0) * (n - 1.0) - m * m) / static_cast<double>((2.0 * n - 1.0) * (2.0 * n - 3.0));
            if (m > 0)
            {
                snorm[m * (MAX_ORDER + 1) + n] = snorm[(m - 1) * (MAX_ORDER + 1) + n] * std::sqrt(static_cast<double>((n - m + 1) * j) / static_cast<double>(n + m));
                j = 1.0;
                c[n][m - 1] *= snorm[m * (MAX_ORDER + 1) + n];
                c_dot[n][m - 1] *= snorm[m * (MAX_ORDER + 1) + n];
            }

            c[m][n] *= snorm[m * (MAX_ORDER + 1) + n];
            c_dot[m][n] *= snorm[m * (MAX_ORDER + 1) + n];
        }

        fn[n] = static_cast<double>(n + 1.0);
        fm[n] = static_cast<double>(n);
    }

    k[1][1] = 0.0;

#ifdef DEBUG

    std::printf("EPOCH: %lf \n\n", epoch);

    std::printf("\nsnorm[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE; j++)
        {
            std::printf("%lf ", snorm[i * 13 + j]);
        }
        std::printf("\n");
    }

    std::printf("\nc[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE + 1; j++)
        {
            std::printf("%lf ", c[i][j]);
        }
        std::printf("\n");
    }

    std::printf("\nc_dot[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE + 1; j++)
        {
            std::printf("%lf ", c_dot[i][j]);
        }
        std::printf("\n");
    }

    std::printf("\fm[]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {

        std::printf("%lf ", fm[i]);
    }

#endif
}
Geo_Mag::~Geo_Mag() {}

void Geo_Mag::E0000(double const &latitude, double const &longitude, double const &altitude, double const &year)
{

    double sp[(MAX_ORDER + 1)] = {0.0}; // sin(m*colat)
    double cp[(MAX_ORDER + 1)] = {0.0}; // cos(m*colat)

    double b_theta = 0.0; // b component in spherical coord (lon)
    double b_phi = 0.0;   // b component in spherical coord (lat)
    double b_r = 0.0;     // b component in spherical coord (r)

    double bpp = 0.0;
    double rlat = latitude * (M_PI / 180.0);
    double rlon = longitude * (M_PI / 180.0);
    double sin_lat = std::sin(rlat);
    double sin_lon = std::sin(rlon);
    double cos_lat = std::cos(rlat);
    double cos_lon = std::cos(rlon);

    double dt = year - epoch;

    /*Convert Geodetics Coordinates to Geodesics (Spherical Coord)*/
    double q = std::sqrt(A2 - C2 * sin_lat * sin_lat);
    double q1 = altitude * q;
    double q2 = (q1 + A2) / (q1 + B2);
    q2 = q2 * q2;
    double ct = sin_lat / std::sqrt(q2 * cos_lat * cos_lat + sin_lat * sin_lat);
    double st = std::sqrt(1.0 - ct * ct);
    double r = std::sqrt(altitude * altitude + 2.0 * q1 + (A4 - (A4 - B4) * sin_lat * sin_lat) / (q * q));

    double P_nm[169] = {0.0};                               // Associate Legendre Polynomials
    double dP_nm[(MAX_ORDER + 1)][(MAX_ORDER + 1)] = {0.0}; // Derivative of Assosciate Legendre Polynomial
    double pp[(MAX_ORDER + 1)] = {0.0};
    double aor = RE / r;   // a/r
    double ar = aor * aor; // (a/r)^2

    // pre-compute the sin(m *lon) and cos(m*lon) using recursive method
    // one can use che Chebyshev method or using trigo identities: sin(a+b), cos(a+b)
    // where a= lat, b = (m-1)*lat
    sp[1] = sin_lon;
    cp[0] = 1.0;
    cp[1] = cos_lon;
    for (int m = 2; m <= MAX_DEGREE; m++)
    {
        sp[m] = sp[1] * cp[m - 1] + cp[1] * sp[m - 1];
        cp[m] = cp[1] * cp[m - 1] - sp[1] * sp[m - 1];
    }

    pp[0] = 1.0;
    P_nm[0] = 1.0;
    for (int n = 1; n <= MAX_ORDER; n++)
    {
        ar = ar * aor; // (a/r)^(n+2)

        for (int m = 0; m <= n; m++)
        {
            /*Computing the unnormalized Associate Legendre Polynomials and its derivatives via recursion relations*/
            if (n == m)
            {
                P_nm[m * (MAX_ORDER + 1) + n] = st * P_nm[(m - 1) * (MAX_ORDER + 1) + n - 1];
                dP_nm[m][n] = st * dP_nm[m - 1][n - 1] + ct * P_nm[(m - 1) * (MAX_ORDER + 1) + n - 1];
            }

            else if ((n == 1) && (m == 0))
            {
                P_nm[m * (MAX_ORDER + 1) + n] = ct * P_nm[m * (MAX_ORDER + 1) + n - 1];
                dP_nm[m][n] = ct * dP_nm[m][n - 1] - st * P_nm[m * (MAX_ORDER + 1) + n - 1];
            }

            else
            {
                if (m > (n - 2))
                {
                    P_nm[m * (MAX_ORDER + 1) + n - 2] = 0.0;
                    dP_nm[m][n - 2] = 0.0;
                }

                P_nm[m * (MAX_ORDER + 1) + n] = ct * P_nm[m * (MAX_ORDER + 1) + n - 1] - k[m][n] * P_nm[m * (MAX_ORDER + 1) + n - 2];
                dP_nm[m][n] = ct * dP_nm[m][n - 1] - st * P_nm[m * (MAX_ORDER + 1) + n - 1] - k[m][n] * dP_nm[m][n - 2];
            }

            // Time adjust the gauss coefficients
            tc[m][n] = c[m][n] + dt * c_dot[m][n];

            if (m != 0)
            {
                tc[n][m - 1] = c[n][m - 1] + dt * c_dot[n][m - 1];
            }

            // accumulate terms of the spherical harmonic expansions
            double temp1 = 0.0;
            double temp2 = 0.0;
            if (m == 0)
            {
                temp1 = tc[m][n] * cp[m];
                temp2 = tc[m][n] * sp[m];
            }
            else
            {
                temp1 = tc[m][n] * cp[m] + tc[n][m - 1] * sp[m];
                temp2 = tc[m][n] * sp[m] - tc[n][m - 1] * cp[m];
            }

            b_theta -= ar * temp1 * dP_nm[m][n];
            b_phi += fm[m] * temp2 * ar * P_nm[m * (MAX_ORDER + 1) + n];
            b_r += fn[n] * temp1 * ar * P_nm[m * (MAX_ORDER + 1) + n];

            // special case: north/south poles
            if ((st == 0.0) && (m == 1))
            {
                if (n == 1)
                {
                    pp[n] = pp[n - 1];
                }
                else
                {
                    pp[n] = ct * pp[n - 1] - k[m][n] * pp[n - 2];
                }

                bpp += fm[m] * temp2 * ar * pp[n];
            }
        }
    }

    if (st == 0.0)
    {
        b_phi = bpp;
    }
    else
    {
        b_phi /= st;
    }

    /*Rotate magnetic vector component from spherical to geodetic coord*/
    double d = std::sqrt(A2 * cos_lat * cos_lat + B2 * sin_lat * sin_lat);
    double ca = (altitude + d) / r;
    double sa = C2 * cos_lat * sin_lat / (r * d);

    b[0] = -b_theta * ca - b_r * sa;
    b[1] = b_phi;
    b[2] = b_theta * sa - b_r * ca;

    // compute characteristics
    horizontal_intensity = std::sqrt(b[0] * b[0] + b[1] * b[1]);
    total_intensity = std::sqrt(horizontal_intensity * horizontal_intensity + b[2] * b[2]);
    inclination = std::atan2(b[2], horizontal_intensity);
    declination = std::atan2(b[1], b[0]);

#ifdef DEBUG
    std::printf("\ntc[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE + 1; j++)
        {
            std::printf("%lf ", tc[i][j]);
        }
        std::printf("\n");
    }

    std::printf("\nP_nm[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE + 1; j++)
        {
            std::printf("%lf ", P_nm[i * 13 + j]);
        }
        std::printf("\n");
    }

    std::printf("\ndP_nm[][]:\n");
    for (int i = 0; i < MAX_ORDER + 1; i++)
    {
        for (int j = 0; j < MAX_DEGREE + 1; j++)
        {
            std::printf("%lf ", dP_nm[i][j]);
        }
        std::printf("\n");
    }

    printf("b_theta: %lf b_phi: %lf b_r: %lf\n\n", b_theta, b_phi, b_r);

#endif
}

double *Geo_Mag::Get_Magnetic_Vector()
{
    return b;
}

double Geo_Mag::Get_Horizontal_Intensity()
{
    return horizontal_intensity;
}

double Geo_Mag::Get_Total_Intensity()
{
    return total_intensity;
}

double Geo_Mag::Get_Inclination()
{
    return inclination;
}

double Geo_Mag::Get_Declination()
{
    return declination;
}