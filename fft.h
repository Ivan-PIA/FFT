#pragma once

#include <complex>
#include <vector>

class FFT {
public:
    using Complex = std::complex<double>;
    static std::vector<Complex> compute(std::vector<Complex>& data, bool inverse);

private:
    static int revers(int num, int l);
    static void compute_fft(std::vector<Complex>& data,std::vector<Complex>& data_fft, size_t N, bool inverse);
    static std::vector<Complex> fft_mixed(const std::vector<Complex>& data, bool inverse);
    static std::vector<Complex> fft_2(const std::vector<Complex>& x, bool inverse);
    static std::vector<Complex> fft_3(const std::vector<Complex>& x, bool inverse);
    static std::vector<Complex> fft_5(const std::vector<Complex>& x, bool inverse);
    static std::vector<Complex> dft(const std::vector<Complex> &data, bool inverse);
    static void fft(const std::vector<Complex>& data,std::vector<Complex>& data_fft,size_t N, bool inverse);

};
