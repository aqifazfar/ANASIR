#include <cstdint>
#include <cassert>
#include <cstdio>

template <typename T, std::uint32_t M, std::uint32_t N>
class Matrix
{
private:
    T _data[M * N]{};

public:
    inline T &operator()(std::uint32_t i, std::uint32_t j)
    {
        return this->_data[i * N + j];
    }

    inline T const &operator()(std::uint32_t i, std::uint32_t j) const
    {
        return this->_data[i * N + j];
    }
    inline T &operator[](std::uint32_t i)
    {
        return this->_data[i];
    }

    inline Matrix<T, M, N> &operator=(Matrix<T, M, N> const &in)
    {

        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                this->_data[i * N + j] = in._data[i * N + j];
            }
        }

        return *this;
    }

    template <std::uint32_t P>
    inline Matrix<T, M, P> operator*(Matrix<T, N, P> const &other) const
    {
        Matrix<T, M, P> res;
        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < P; j++)
            {
                for (std::size_t k = 0; k < N; k++)
                {
                    res(i, j) += this->_data[i * N + k] * other.Data()[k * P + j];
                }
            }
        }
        return res;
    }

    inline Matrix<T, M, N> operator*(T const &scalar) const
    {
        Matrix<T, M, N> res;
        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                res(i, j) = scalar * this->_data[i * N + j];
            }
        }
        return res;
    }

    friend Matrix<T, M, N> operator*(T const &scalar, Matrix<T, M, N> const &mat)
    {
        return mat * scalar;
    }

    inline Matrix<T, M, N> operator+(Matrix<T, M, N> const &other) const
    {

        Matrix<T, M, N> res;
        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                res._data[i * N + j] = this->_data[i * N + j] + other._data[i * N + j];
            }
        }
        return res;
    }

    inline Matrix<T, M, N> operator-(Matrix<T, M, N> const &other) const
    {

        Matrix<T, M, N> res;
        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                res._data[i * N + j] = this->_data[i * N + j] - other._data[i * N + j];
            }
        }
        return res;
    }

    inline T *Data()
    {
        return this->_data;
    }

    inline T const *Data() const
    {
        return this->_data;
    }

    inline void Print()
    {

        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                std::printf("%.2f ", this->_data[i * N + j]);
            }
            std::printf("\n");
        }
    }

    inline std::uint32_t Rows()
    {
        return M;
    }

    inline std::uint32_t Cols()
    {
        return N;
    }

    inline void Fill(T val)
    {

        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                this->_data[i * N + j] = val;
            }
        }
    }

    inline void Eye()
    {
        this->Fill(0.0);
        for (std::size_t i = 0; i < M; i++)
        {

            this->_data[i * N + i] = 1;
        }
    }

    inline void Diagonal(T const &val)
    {
        this->Fill(0.00);

        for (std::size_t i = 0; i < N; i++)
        {
            this->_data[i * N + i] = val;
        }
    }

    template <std::uint32_t rows, std::uint32_t cols>
    inline void Copy_Block(Matrix<T, rows, cols> const &in, std::uint32_t stxRow, std::uint32_t stxCol)
    {
        assert(((stxRow + rows) <= M) && ((stxCol + cols) <= N));

        for (std::size_t i = 0; i < rows; i++)
        {
            std::size_t row_offset = (stxRow + i) * N;
            for (std::size_t j = 0; j < cols; j++)
            {
                this->_data[row_offset + stxCol + j] = in(i, j);
            }
        }
    }

    inline Matrix<T, N, M> Transpose()
    {
        Matrix<T, N, M> res;

        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                res(j, i) = this->_data[i * N + j];
            }
        }

        return res;
    }

    inline void Inverse()
    {
        Matrix<T, N, N> L;
        Matrix<T, N, N> U;
        Matrix<T, N, N> L_inverse;
        Matrix<T, N, N> U_inverse;
        Matrix<T, N, N> temp;

        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                temp(i, j) = this->_data[i * N + j];
            }
        }

        for (std::size_t i = 0; i < N; i++)
        {
            L(i, i) = 1.0;
            for (std::size_t j = 0; j < N; j++)
            {
                T sum = 0.00;

                if (i <= j)
                {
                    // Compute U
                    for (std::size_t k = 0; k < i; k++)
                    {
                        sum += L(i, k) * U(k, j);
                    }

                    U(i, j) = temp(i, j) - sum;
                }
                else
                {
                    // Compute L
                    for (std::size_t k = 0; k < j; k++)
                    {
                        sum += L(i, k) * U(k, j);
                    }

                    L(i, j) = (temp(i, j) - sum) / U(j, j);
                }
            }
        }

        // Compute inverse
        for (std::size_t i = 0; i < N; i++)
        {
            L_inverse(i, i) = 1.0 / L(i, i);
            U_inverse(i, i) = 1.0 / U(i, i);
        }

        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                // for lower
                if (i > j)
                {
                    T sum_l = 0.00;
                    for (std::size_t k = j; k < i; k++)
                    {
                        sum_l += L(i, k) * L_inverse(k, j);
                    }
                    L_inverse(i, j) = -sum_l / L(i, i);
                }

                // for upper
                if (i < j)
                {
                    T sum_u = 0.0;

                    for (std::size_t k = i + 1; k <= j; k++)
                    {
                        sum_u += U(i, k) * U_inverse(k, j);
                    }

                    U_inverse(i, j) = -sum_u / U(i, i);
                }
            }
        }

        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                this->_data[i * N + j] = 0.0;
                for (std::size_t k = 0; k < N; k++)
                {
                    this->_data[i * N + j] += L_inverse(i, k) * U_inverse(k, j);
                }
            }
        }
    }
};

template <typename T, std::uint32_t N>
class Vector : public Matrix<T, N, 1>
{
public:
    using Matrix<T, N, 1>::operator=;

    inline Matrix<T, N, N> hat()
    {
        assert((N == 3) && "N is not equal to 3");

        Matrix<T, N, N> res;

        res(0, 0) = 0.0;
        res(0, 1) = -this->Data()[2];
        res(0, 2) = this->Data()[1];

        res(1, 0) = this->Data()[2];
        res(1, 1) = 0.0;
        res(1, 2) = -this->Data()[0];

        res(2, 0) = -this->Data()[1];
        res(2, 1) = this->Data()[0];
        res(2, 2) = 0.0;

        return res;
    }

    inline void Normalize()
    {
        T sum = 0.0;
        for (std::size_t i = 0; i < N; i++)
        {
            sum += this->Data()[i] * this->Data()[i];
        }

        if (sum == 0.0)
        {
            // Skip normalization
            return;
        }

        T magnitude = std::sqrt(sum);

        for (std::uint32_t i = 0; i < N; i++)
        {
            this->_data[i] = this->_data[i] / magnitude;
        }
    }

    inline T Dot(Vector<T, N> const &other)
    {
        T res = 0.0;
        for (std::size_t i = 0; i < N; i++)
        {
            res += this->_data[i] * other[i];
        }

        return res;
    }

    inline Vector<T, N> Cross(Vector<T, N> const &other)
    {
        assert(N == 3);

        Vector<T, N> res;
        res[0] = this->_data[1] * other[2] - this->_data[2] * other[1];
        res[1] = this->_data[2] * other[0] - this->_data[0] * other[2];
        res[2] = this->_data[0] * other[1] - this->_data[1] * other[0];
        return res;
    }

    inline Matrix<T, N, N> To_Diagonal()
    {

        Matrix<T, N, N> res;

        for (std::size_t i = 0; i < N; i++)
        {
            res(i, i) = this->Data()[i];
        }

        return res;
    }
};

// #include <ctime>
// #include <iostream>
// std::size_t main()
// {
//     Matrix<T, 3, 3> a;

//     a(0, 0) = 2.0;
//     a(0, 1) = 0.0;
//     a(0, 2) = -1.0;
//     a(1, 0) = 5.0;
//     a(1, 1) = 1.0;
//     a(1, 2) = 0.0;
//     a(2, 0) = 0.0;
//     a(2, 1) = 1.0;
//     a(2, 2) = 3.0;

//     // a.inverse();
//     a.print();
// }