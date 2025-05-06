#include "fft.h"
#include <iostream>
#include <cmath>


const double PI = 3.141592653589793;


std::vector<FFT::Complex> FFT::compute(std::vector<Complex>& data, bool inverse) {
    int N = data.size();

    std::vector<Complex> data_fft(N);
    data_fft = fft_mixed(data,inverse);

    if(inverse){
        for(Complex &x:data_fft){
            x /= N; 
        }
    }

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


std::vector<FFT::Complex> FFT::fft_mixed(const std::vector<Complex>& data, bool inverse) {

    int N = data.size();
    // std::cout << "size = "<< N << std::endl;

    if (N % 2 == 0){ 
        // std::cout << "Using fft_2" << std::endl;
        return fft_2(data, inverse);
    }
    else if (N % 3 == 0){
        // std::cout << "Using fft_3" << std::endl;
        return fft_3(data, inverse);
    }
    else if (N % 5 == 0){
        // std::cout << "Using fft_5" << std::endl;
        return fft_5(data, inverse);
    }

    return dft(data, inverse);
}


std::vector<FFT::Complex> FFT::fft_2(const std::vector<Complex>& data, bool inverse) {
    int N = data.size();

    std::vector<Complex> part1(N / 2);
    std::vector<Complex> part2(N / 2);
    for (int i = 0; i < N / 2; i++) {
        part1[i] = data[2 * i];
        part2[i] = data[2 * i + 1];
    }

    std::vector<Complex> fft_part1 = fft_mixed(part1,inverse);
    std::vector<Complex> fft_part2 = fft_mixed(part2,inverse);

    std::vector<Complex> result(N);
    for (int k = 0; k < N / 2; k++) {
        double ang = (2 * PI * (inverse ? 1 : -1) * k) / N;
        Complex w = std::polar(1.0, ang);
        result[k]         = fft_part1[k] + w * fft_part2[k];
        result[k + N / 2] = fft_part1[k] - w * fft_part2[k];
    }

    return result;
}

std::vector<FFT::Complex> FFT::fft_3(const std::vector<Complex>& data, bool inverse) {
    
    int N = data.size();
    
    std::vector<Complex> part1(N / 3);
    std::vector<Complex> part2(N / 3);
    std::vector<Complex> part3(N / 3);

    for (int i = 0; i < N / 3; i++) {
        part1[i] = data[3 * i];
        part2[i] = data[3 * i + 1];
        part3[i] = data[3 * i + 2];
    }

    std::vector<Complex> fft_part1 = fft_mixed(part1,inverse);
    std::vector<Complex> fft_part2 = fft_mixed(part2,inverse);
    std::vector<Complex> fft_part3 = fft_mixed(part3,inverse);

    std::vector<Complex> result(N);

    for (int k = 0; k < N / 3; k++) {
        Complex ws[5][4]; 
        int sign = inverse ? 1 : -1;
    
        for (int i = 0; i < 3; i++) {
            int shift = k + i * N / 3;
            for (int j = 0; j < 2; j++) {
                double angle = 2.0 * PI * sign * (j+1) * shift / N;
                ws[i][j] = std::polar(1.0, angle);
            }
        }
    
        Complex a = fft_part1[k];
        Complex b = fft_part2[k];
        Complex c = fft_part3[k];
    
        for (int i = 0; i < 3; i++) {
            result[k + i * N / 3] = a + ws[i][0] * b + ws[i][1] * c;
        }
    }

    return result; 

}

std::vector<FFT::Complex> FFT::fft_5(const std::vector<Complex>& data, bool inverse) {

    int N = data.size();
    
    std::vector<Complex> part1(N / 5);
    std::vector<Complex> part2(N / 5);
    std::vector<Complex> part3(N / 5);
    std::vector<Complex> part4(N / 5);
    std::vector<Complex> part5(N / 5);

    for (int i = 0; i < N / 5; i++) {
        part1[i] = data[5 * i];
        part2[i] = data[5 * i + 1];
        part3[i] = data[5 * i + 2];
        part4[i] = data[5 * i + 3];
        part5[i] = data[5 * i + 4];
    }

    std::vector<Complex> fft_part1 = fft_mixed(part1,inverse);
    std::vector<Complex> fft_part2 = fft_mixed(part2,inverse);
    std::vector<Complex> fft_part3 = fft_mixed(part3,inverse);
    std::vector<Complex> fft_part4 = fft_mixed(part4,inverse);
    std::vector<Complex> fft_part5 = fft_mixed(part5,inverse);

    std::vector<Complex> result(N);

    for (int k = 0; k < N / 5; k++) {
        Complex ws[5][4]; // матрица поворотных коэффициентов
        int sign = inverse ? 1 : -1;
    
        for (int i = 0; i < 5; i++) {
            int shift = k + i * N / 5;
            for (int j = 1; j <= 4; j++) {
                double angle = 2.0 * PI * sign * j * shift / N;
                ws[i][j - 1] = std::polar(1.0, angle);
            }
        }
    
        Complex a = fft_part1[k];
        Complex b = fft_part2[k];
        Complex c = fft_part3[k];
        Complex d = fft_part4[k];
        Complex e = fft_part5[k];
    
        for (int i = 0; i < 5; ++i) {
            result[k + i * N / 5] = a
                + ws[i][0] * b
                + ws[i][1] * c
                + ws[i][2] * d
                + ws[i][3] * e;
        }
    }

    return result; 

}


std::vector<FFT::Complex> FFT::dft(const std::vector<Complex> &data, bool inverse){

    int N = data.size();
    std::vector<Complex> data_dft(N);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            double ang = (2 * PI * (inverse ? 1 : -1) * i * j) / N;
            Complex exp = std::polar(1.0, ang);
            data_dft[i] += data[j] * exp;
        }
    }

    return data_dft;
}



void FFT::fft(const std::vector<Complex>& data,std::vector<Complex>& data_fft, size_t N, bool inverse) { 
    
    data_fft = data;
    int log2 = std::log2(N);
    std::cout<< "log2(N) = " << log2 << std::endl;
    
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




