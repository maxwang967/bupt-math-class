//
// Created by 余京泽 on 2020/5/25.
//

#ifndef PROJECT_FFT_MY_FFT_H
#define PROJECT_FFT_MY_FFT_H

#include "Complex.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>

const double PI = 3.14159265358979;
const double PI_X2 = 2 * PI;

using namespace ComplexNumber;
extern void FFT(Complex* dst, Complex* src, int p);
extern void iFFT(Complex* dst, Complex* src, int p);
extern void fft(Complex* dst, Complex* src, int p);
extern void ifft(Complex* dst, Complex* src, int p);

typedef std::vector<unsigned char> Integer;
extern void multiply(Integer & rst, Integer const& a, Integer const& b);

#endif //PROJECT_FFT_MY_FFT_H
