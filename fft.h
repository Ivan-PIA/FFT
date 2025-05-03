#pragma once

#include <complex>
#include <vector>

class FFT {
public:
    using Complex = std::complex<double>;
    static std::vector<Complex> compute(std::vector<Complex>& data, bool inverse);

private:
    static int revers(int num, int l);
    static bool is_valid_size(size_t N);
    static void compute_fft(std::vector<Complex>& data,std::vector<Complex>& data_fft, size_t N, bool inverse);
    static void fft(const std::vector<Complex>& data,std::vector<Complex>& data_fft,size_t N, bool inverse);
    static std::vector<FFT::Complex> fft2(const std::vector<Complex> &num, bool invert);
    static size_t smallest_radix(size_t N);
};
