#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <iomanip>

#include "fft.h"

int main() {
    using Complex = std::complex<double>;

    int N = 21;  // длина, кратная 2 3 5
    std::vector<Complex> data(N);

    std::mt19937 rng(1);  // фиксируем сид
    std::uniform_int_distribution<int> bit_dist(0, 1);

    // QPSK 
    auto bits_to_qpsk = [](int b1, int b0) -> Complex {
        double norm = std::sqrt(0.5);
        double real = (b1 ? -1 : 1);
        double imag = (b0 ? -1 : 1);
        return Complex(real, imag) * norm;
    };

    
    for (size_t i = 0; i < N; ++i) {
        int b1 = bit_dist(rng);
        int b0 = bit_dist(rng);
        data[i] = bits_to_qpsk(b1, b0);
    }

    std::cout << "\nSize = " << N << "\n";

    std::cout << "\nQPSK сигнал:\n";
    for (const auto& x : data) {
        std::cout << "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    // FFT
    
    std::vector<Complex> data_ifft = FFT::compute(data, false);

    // Печать спектра (амплитуда)
    std::cout << "\nFFT:\n\n";
    for (const auto& x : data_ifft) {
        std::cout <<  "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    // IFFT
    std::vector<Complex> data_fft = FFT::compute(data_ifft, true);

    std::cout << "\nIFFT:\n\n";
    
    for (const auto& x : data_fft) {
        std::cout << "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    return 0;
}
