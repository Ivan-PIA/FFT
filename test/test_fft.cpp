// test_fft.cpp
#include <gtest/gtest.h>
#include "../fft.h"
#include <complex>
#include <vector>
#include <cmath>
#include <random>

using Complex = std::complex<double>;

bool are_equal(const Complex& a, const Complex& b, double epsilon = 1e-9) {
    return std::abs(a.real() - b.real()) < epsilon && std::abs(a.imag() - b.imag()) < epsilon;
}

std::vector<Complex> generate_random_qpsk(int N) {
    std::vector<Complex> data(N);
    std::mt19937 rng(1);
    std::uniform_int_distribution<int> bit_dist(0, 1);

    auto bits_to_qpsk = [](int b1, int b0) -> Complex {
        double norm = std::sqrt(0.5);
        double real = (b1 ? -1 : 1);
        double imag = (b0 ? -1 : 1);
        return Complex(real, imag) * norm;
    };

    for (int i = 0; i < N; ++i) {
        int b1 = bit_dist(rng);
        int b0 = bit_dist(rng);
        data[i] = bits_to_qpsk(b1, b0);
    }
    return data;
}

TEST(FFTTest, ForwardInverseMatchQPSK) {

    std::vector<int> sizes = {4, 8, 9, 25, 18, 123, 2048, 1024, 22, 55};

    for(int N : sizes){
        auto input = generate_random_qpsk(N);

        auto fft = FFT::compute(input, false);  // FFT
        auto ifff = FFT::compute(fft, true);  // IFFT

        for (int i = 0; i < N; ++i) {
            EXPECT_TRUE(are_equal(input[i], ifff[i]))
                << "Mismatch at index " << i << ": "
                << input[i] << " vs " << ifff[i];
        }
    }
}
