#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <iomanip>

#include "fft.h"

int main() {
    using Complex = std::complex<double>;

    const size_t N = 128;  // длина, кратная 2, 3, 5
    std::vector<Complex> data(N);

    std::mt19937 rng(123);  // фиксируем seed
    std::uniform_int_distribution<int> bit_dist(0, 1);

    // QPSK mapping: 2 бита -> 1 комплексная точка
    auto bits_to_qpsk = [](int b1, int b0) -> Complex {
        double norm = std::sqrt(0.5);
        double real = (b1 ? -1 : 1);
        double imag = (b0 ? -1 : 1);
        return Complex(real, imag) * norm;
    };

    // Генерация N QPSK символов
    for (size_t i = 0; i < N; ++i) {
        int b1 = bit_dist(rng);
        int b0 = bit_dist(rng);
        data[i] = bits_to_qpsk(b1, b0);
    }

    // Печать QPSK сигнала
    std::cout << "QPSK сигнал:\n";
    for (const auto& x : data) {
        std::cout << std::fixed << std::setprecision(4)
                  << "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    // FFT
    
    std::vector<Complex> data_fft = FFT::compute(data, false);

    // Печать спектра (амплитуда)
    std::cout << "\nСпектр (амплитуды):\n";
    for (const auto& x : data_fft) {
        std::cout << std::fixed << std::setprecision(4)
                  << "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    
    std::vector<Complex> data_ifft = FFT::compute(data_fft, true);

    std::cout << "demod QPSK сигнал:\n";
    
    for (const auto& x : data_ifft) {
        std::cout << std::fixed << std::setprecision(4)
                  << "(" << x.real() << ", " << x.imag() << "i)\n";
    }

    return 0;
}
