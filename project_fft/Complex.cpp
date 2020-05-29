//
// Created by 余京泽 on 2020/5/23.
//

#include "Complex.h"
using namespace std;
using namespace ComplexNumber;
Complex::Complex():real(0),imag(0)
{
}

Complex::Complex(double Re, double Im):real(Re),imag(Im)
{
}

Complex Complex::operator+(const Complex& cplx)
{
    Complex res;
    res.real = this->real + cplx.real;
    res.imag = this->imag + cplx.imag;
    return res;
}

Complex Complex::operator-(const Complex& cplx)
{
    Complex res;
    res.real = this->real - cplx.real;
    res.imag = this->imag - cplx.imag;
    return res;
}

Complex Complex::operator*(const Complex& cplx)
{
    return multiply(*this,cplx);
}
Complex Complex::operator/(const Complex &cplx)
{
    return divide(*this,cplx);
}
Complex& Complex::operator=(const Complex& cplx)
{
    real = cplx.real;
    imag = cplx.imag;
    return *this;
}

Complex& Complex::operator+= (const Complex& cplx)
{
    real += cplx.real;
    imag += cplx.imag;
    return *this;
}

Complex& Complex::operator-= (const Complex& cplx)
{
    real -= cplx.real;
    imag -= cplx.imag;
    return *this;
}

Complex& Complex::operator*= (const Complex& cplx)
{
    *this = multiply(*this, cplx);
    return *this;
}

Complex& Complex::operator/= (const Complex& cplx)
{
    *this = divide(*this, cplx);
    return *this;
}

bool ComplexNumber::operator==(const Complex & c1, const Complex & c2)
{
    if (c1.real == c2.real&&c1.imag == c2.imag)return true;
    else return false;
}

Complex Complex::multiply(const Complex c1, const Complex c2)
{
    Complex res;
    res.real = c1.real*c2.real - c1.imag*c2.imag;
    res.imag = c1.real*c2.imag + c1.imag*c2.real;
    return res;

}

Complex Complex::divide(const Complex c1, const Complex c2)
{
    Complex res;
    res.real = (((c1.real*c2.real) + (c1.imag*c2.imag)) / (c2.real*c2.real + c2.imag*c2.imag));
    res.imag = (((c1.imag*c2.real) - (c1.real*c2.imag)) / (c2.real*c2.real + c2.imag*c2.imag));
    return res;
}

ostream & ComplexNumber::operator<<(ostream& os, const Complex& cplx)
{
    if (cplx.imag== 0.0) {
        os << cplx.real;
    }
    else
    {
        if (cplx.real == 0.0) {
            os << cplx.imag << "i";
        }
        else {
            if (cplx.imag < 0) {
                os << cplx.real << cplx.imag << "i";
            }
            else {
                os << cplx.real << "+" << cplx.imag << "i";
            }

        }
    }
    return os;
}


istream& ComplexNumber::operator >> (istream& is, Complex& cplx)
{
    is >> cplx.real >> cplx.imag;
    return is;
}