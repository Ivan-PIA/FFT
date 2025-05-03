#include "fft.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

const double PI = 3.141592653589793;


std::vector<FFT::Complex> FFT::compute(std::vector<Complex>& data, bool inverse) {
    size_t N = data.size();
    
    if (N == 0 || !is_valid_size(N)) {
        throw std::invalid_argument("FFT size must be a product of 2, 3, and/or 5");
    }

    std::vector<Complex> data_fft(N);
    fft(data, data_fft, N, inverse);
    // data_fft = fft2(data,inverse );
    return data_fft;
}


int FFT::revers(int idx, int log2_N){
    int revers_ind = 0;

    for(int i = 0; i < log2_N; i++){
        revers_ind = revers_ind << 1;
        revers_ind = revers_ind | (idx >> i) & 1;
    }
    
    return revers_ind;
}

bool FFT::is_valid_size(size_t N) {
    while (N % 2 == 0) N /= 2;
    while (N % 3 == 0) N /= 3;
    while (N % 5 == 0) N /= 5;
    return N == 1;
}

void FFT::fft(const std::vector<Complex>& data,std::vector<Complex>& data_fft, size_t N, bool inverse) { 
    data_fft = data;
    
    int log2 = std::log2(N);

    // for(int i = 0; i < N; i++){
    //     data_fft[revers(i,log2)] = data_fft[i];
    // }
    
    for (int i = 0; i < N; i++)
    {
        if (i < revers(i, log2))
            std::swap(data_fft[i], data_fft[revers(i, log2)]);
    }

    for(int l = 2; l <= N; l*=2){
        double ang = 2 * PI / l * (inverse ? 1 : -1);
        Complex exp = std::polar(1.0, ang);
        
        for(int i = 0; i < N; i +=l){
            Complex w(1);
            for(int j = 0; j < l/2; j++ ){
                Complex u = data_fft[i+j];
                Complex v = data_fft[i+j+l/2] * w;
                data_fft[i+j] = u + v;
                data_fft[i+j+l/2] = u - v;
                w *= exp;
            }
        }
    }

    if (inverse)
    {
        for (Complex & x : data_fft)
            x /= N;
    }

}


size_t FFT::smallest_radix(size_t N) {
    if (N % 2 == 0) return 2;
    if (N % 3 == 0) return 3;
    if (N % 5 == 0) return 5;
    return N;
}




