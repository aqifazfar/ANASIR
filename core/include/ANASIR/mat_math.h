#include <iostream>
#include <cstdint>
#include <cmath>

namespace Matrix
{
	template <std::uint32_t m, std::uint32_t n>
	void Sub(double (&dst)[m * n], const double (&A)[m * n], const double (&B)[m * n])
	{
		for (std::uint32_t i = 0; i < (m * n); i++)
		{
			dst[i] = A[i] - B[i];
		}
	}

	template <std::uint32_t m, std::uint32_t n>
	void Add(double (&dst)[m * n], const double (&A)[m * n], const double (&B)[m * n])
	{
		for (std::uint32_t i = 0; i < (m * n); i++)
		{
			dst[i] = A[i] + B[i];
		}
	}

	template <std::uint32_t ARows, std::uint32_t AColumns, std::uint32_t BColumns>
	void Mul(double (&dst)[ARows * BColumns], const double (&A)[ARows * AColumns], const double (&B)[AColumns * BColumns])
	{
		std::uint32_t i, j, k;

		for (i = 0; i < ARows; i++)
		{
			for (j = 0; j < BColumns; j++)
			{
				dst[i * BColumns + j] = 0;

				for (k = 0; k < AColumns; k++)
				{
					dst[i * BColumns + j] += A[i * AColumns + k] * B[k * BColumns + j];
				}
			}
		}
	}
	template <std::uint32_t n>
	void CholeskyDecomposition(double (&dst)[n * n], const double (&A)[n * n])
	{
		std::uint32_t i, j, k;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j <= i; j++)
			{
				dst[i * n + j] = 0;

				double sum = 0;

				for (k = 0; k < j; k++)
				{
					sum += dst[i * n + k] * dst[j * n + k];
				}

				if (i == j)
				{
					dst[i * n + j] = sqrt(A[i * n + j] - sum);
				}
				else
				{
					dst[i * n + j] = (1 / dst[j * n + j]) * (A[i * n + j] - sum);
				}
			}
		}
	}

	template <std::uint32_t n>
	void Inv(double (&dst)[n * n], const double (&A)[n * n])
	{
		std::uint32_t i, j, k;
		double sum;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (i == j)
				{
					dst[i * n + j] = 1 / A[i * n + j];
				}
				else
				{
					sum = 0;

					for (k = j; k < i; k++)
					{
						sum += A[i * n + k] * dst[k * n + j];
					}
					dst[i * n + j] = -sum / A[i * n + i];
				}
			}
		}
	}

	template <std::uint32_t m, std::uint32_t n>
	void Scale(double (&dst)[m * n], const double (&A)[m * n], const double Scalar)
	{
		std::uint32_t i;
		for (i = 0; i < (m * n); i++)
		{
			dst[i] = A[i] * Scalar;
		}
	}

	template <std::uint32_t m, std::uint32_t n>
	void Transpose(double (&A)[m * n])
	{
		std::uint32_t i, j;
		double temp[m * n] = {};
		std::memcpy(temp, A, sizeof(A));

		for (i = 0; i < m; i++)
		{
			for (j = 0; j < n; j++)
			{
				A[i * m + j] = temp[j * m + i];
			}
		};
	}

	template <std::uint32_t n>
	void Diagonal(double (&dst)[n * n], double val)
	{
		std::uint32_t i, j;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (i == j)
				{
					dst[i * n + j] = val;
				}
				else
				{
					dst[i * n + j] = 0;
				}
			}
		}
	}
}