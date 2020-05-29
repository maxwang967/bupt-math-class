//
// Created by 余京泽 on 2020/5/25.
//

#include "my_fft.h"

/* r=log2(N) */
void FFT(Complex *TD, Complex *FD, int r)
{
    const int count = 1 << r;

    int csz = sizeof(Complex)*count;
    auto* W = (Complex*)malloc(csz / 2);  //存基底
    auto* X1 = (Complex*)malloc(csz);
    auto* X2 = (Complex*)malloc(csz);
    Complex* X = nullptr;

    int  i, j, k;
    int  dist, p;
    double f = PI_X2 / count;
    double a = 0;
    for (i = 0; i < count / 2; ++i)
    {
        W[i] = Complex(cos(a), -sin(a));
        a += f;
    }

    for (i = 0; i < count; ++i)
    {
        X1[i] = TD[i];
    }

    for (k = 0; k < r; ++k)
    {
        for (j = 0; j < (1 << k); ++j)
        {
            dist = 1 << (r - k);
            for (i = 0; i < dist / 2; ++i)
            {
                p = j * dist;
                X2[i + p] = X1[i + p] + X1[i + p + dist / 2];
                X2[i + p + dist / 2] = (X1[i + p] - X1[i + p + dist / 2])* W[i * (1 << k)];
            }
        }
        X = X1;
        X1 = X2;
        X2 = X;
    }

    for (j = 0; j < count; ++j)
    {
        p = 0;
        for (i = 0; i < r; ++i)
        {
            if (j&(1 << i))
            {
                p += 1 << (r - i - 1);
            }
        }
        FD[j] = X1[p];
    }

    free(W);
    free(X1);
    free(X2);
}

/* r=log2(N) */
void iFFT(Complex *TD, Complex *FD, int r)
{
    const int count = 1 << r;
    int csz = sizeof(Complex)*count;
    auto* W = (Complex*)malloc(csz / 2);  //存基底
    auto* X1 = (Complex*)malloc(csz);
    auto* X2 = (Complex*)malloc(csz);
    Complex* X = nullptr;

    int  i, j, k;
    int  dist, p;
    double f = PI_X2 / count;
    double a = 0;
    for (i = 0; i < count / 2; ++i)
    {
        W[i] = Complex(cos(a), sin(a));
        a += f;
    }

    for (i = 0; i < count; ++i)
    {
        X1[i] = TD[i];
    }

    for (k = 0; k < r; ++k)
    {
        for (j = 0; j < (1 << k); ++j)
        {
            dist = 1 << (r - k);
            for (i = 0; i < dist / 2; ++i)
            {
                p = j * dist;
                X2[i + p] = X1[i + p] + X1[i + p + dist / 2];
                X2[i + p + dist / 2] = (X1[i + p] - X1[i + p + dist / 2])* W[i * (1 << k)];
            }
        }
        X = X1;
        X1 = X2;
        X2 = X;
    }

    for (j = 0; j < count; ++j)
    {
        p = 0;
        for (i = 0; i < r; ++i)
        {
            if (j&(1 << i))
            {
                p += 1 << (r - i - 1);
            }
        }
        FD[j] = X1[p]/Complex(count,0);
    }

    free(W);
    free(X1);
    free(X2);
}

void multiply(Integer & rst, Integer const& a, Integer const& b){
    std::vector<int> this_a,this_b;
    for (int i=0;i<a.capacity();i++){
        this_a.push_back(a[i]-'0');
//        this_b.push_back(pow(2.0,b[i]-'0'));
        this_b.push_back(b[i]-'0');
    }

    int lenStr1,lenStr2,sumN;
    lenStr1 = this_a.size();
    lenStr2 = this_b.size();;
    sumN = lenStr1+lenStr2;
    int N=1;
    for(int i = 1;i < sumN;i*=2) {
        N *=2;
    }
//    auto* num1 = new Complex[N];
//    auto* num2 = new Complex[N];
//    auto* num3 = new Complex[N];
//    auto* num1_f = new Complex[N];
//    auto* num2_f = new Complex[N];
//    auto* num3_f = new Complex[N];
    Complex num1[N],num2[N],num3[N],num1_f[N],num2_f[N],num3_f[N];
    for (int i = lenStr1-1;i >= 0 ;i--) {
        num1[lenStr1-i-1] = Complex((double)(this_a[i]),0);
    }
    for (int i = lenStr2-1;i >= 0 ;i--) {
        num2[lenStr2-i-1] = Complex((double)(this_b[i]),0);
    }
    FFT(num1,num1,log2(N));
    FFT(num2,num2,log2(N));
    for(int i = 0;i < N;i++) {
        num3[i] = num1[i]*num2[i];
    }
    iFFT(num3,num3,log2(N));
    int curr,roundup=0;
    for(int i = 0;i < sumN;i++){
        curr = round(num3[i].Rez()+roundup);
        rst.push_back((curr)%10);
        roundup = (curr)/10;
    }
//    vector<int>::reverse_iterator rit = rst.rbegin();
//    while ((*rit) == 0){
//        rst.pop_back();
//        rit = rst.rbegin();
//    }
//    stringstream rststream;
//    for (;rit != rst.rend(); rit++) {
//        rststream<<*rit;//+'0';
//    }
//    string result;
//    rststream>>result;
//    cout<<endl<<"两个乘数运算结果:  "<<result;

}

void fft(Complex* dist, Complex* src, int N) {
    int n = 0;
    for (int i = 1;i < N;i*=2) {//求N的二进制位数
        n++;
    }
    for (int i = 0;i <= N-1;i++) {//进行位反转置换
        int a = i;
        int b = 0;
        for (int j = 0;j < n;j++) {//生成a的反转b
            b = (b<<1)+(a&1);
            a >>= 1;
        }
        if(b > i) {//进行置换
            dist[b] = src[b]+src[i];
            dist[i] = dist[b]-src[i];
            dist[b] = dist[b]-dist[i];
        }
    }
    for (int s = 1, m = 1;s <= n;s++) {//进行迭代过程
        m *= 2;
        Complex temp,u,omiga,omigaM = Complex(cos(-2*M_PI/m),sin(-2*M_PI/m));
        for (int k = 0;k < N; k = k+m) {
            omiga = Complex(1,0);
            for (int j = 0;j <= m/2-1;j++) {//蝶形运算
                temp = omiga * dist[k + j + m/2];
                u = dist[k+j];
                dist[k+j] = u + temp;
                dist[k+j+m/2] = u - temp;
                omiga = omiga*omigaM;
            }
        }
    }
}

void ifft(Complex* dist, Complex* src, int N) {
    int n = 0;
    for (int i = 1;i < N;i*=2) {//求N的二进制位数
        n++;
    }
    for (int i = 0;i <= N-1;i++) {//进行位反转置换
        int a = i;
        int b = 0;
        for (int j = 0;j < n;j++) {//生成a的反转b
            b = (b<<1)+(a&1);
            a >>= 1;
        }
        if(b > i) {//进行置换
            dist[b] = src[b]+src[i];
            dist[i] = dist[b]-src[i];
            dist[b] = dist[b]-dist[i];
        }
    }
    for (int s = 1, m = 1;s <= n;s++) {//进行迭代过程
        m *= 2;
        Complex temp,u,omiga,omigaM = Complex(cos(2*M_PI/m),sin(2*M_PI/m));
        for (int k = 0;k < N; k = k+m) {
            omiga = Complex(1,0);
            for (int j = 0;j <= m/2-1;j++) {//蝶形运算
                temp = omiga * dist[k + j + m/2];
                u = dist[k+j];
                dist[k+j] = u + temp;
                dist[k+j+m/2] = u - temp;
                omiga = omiga*omigaM;
            }
        }
    }
    for(int i = 0;i < N;i++) {
        dist[i] = dist[i]/Complex(N,0);
    }
}