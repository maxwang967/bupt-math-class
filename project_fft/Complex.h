//
// Created by 余京泽 on 2020/5/23.
//

#ifndef PROJECT_FFT_COMPLEX_H
#define PROJECT_FFT_COMPLEX_H

#pragma once
#include <iostream>
#include <cmath>
using namespace std;
namespace ComplexNumber {
    class Complex
    {
        friend ostream& operator<<(ostream&, const Complex&);
        friend istream& operator >> (istream&,Complex&);
        friend bool operator==(const Complex&, const Complex&);
        friend Complex operator+(double a, Complex z) { return Complex(a + z.Rez(), z.Imz()); }
        friend Complex operator*(double a, Complex z) { return Complex(a*z.Rez(), a*z.Imz()); }
    public:
        Complex();
        Complex(double, double);
        double Abs() { return sqrt(real*real + imag*imag); }
        double Rez() { return real; }
        double Imz() { return imag; }
        bool PureIm() { return real == 0 ? true : false; }

        Complex operator+(const Complex&);
        Complex operator-(const Complex&);
        Complex operator*(const Complex&);
        Complex operator/(const Complex&);
        Complex& operator=(const Complex&);
        Complex& operator+= (const Complex&);
        Complex& operator-= (const Complex&);
        Complex& operator*= (const Complex&);
        Complex& operator/= (const Complex&);
        Complex operator~ () { return Complex(real, -imag); }

    private:
        double real;
        double imag;
        Complex multiply(const Complex, const Complex);
        Complex divide(const Complex, const Complex);
    };
}

#endif //PROJECT_FFT_COMPLEX_H
