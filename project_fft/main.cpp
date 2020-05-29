#include <iostream>
#include <vector>
#include <sstream>
#include "Complex.h"
#include "my_fft.h"

int main() {
    using namespace ComplexNumber;
    static Complex j(0, 1);

//    Complex c = 1.4 + 3.1*i;
//    Complex test = 1 + 1.1*i;
//    Complex test1 = test * c;
//    Complex	conj = ~c;//求共轭复数
//    cout << conj << endl;
//    cout << conj*c << endl;

    int length_signal = 8;
//    auto* const signal = (Complex*)malloc(sizeof(Complex)*length_signal);
//    auto* const signal_out = (Complex*)malloc(sizeof(Complex)*length_signal);
    Complex signal1[length_signal], signal_out1[length_signal];
//    float s[16] = {0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1};
    float s[8] = {1,2,3,4,0,0,0,0};
    for (int i=0;i<length_signal;i++){
//        signal[i] = Complex(s[i],0);
        signal1[i] = Complex(s[i],0);
    }

    FFT(signal1, signal_out1, log2(length_signal));

    for(int i=0;i<length_signal;i++){
        printf("%f\n", signal_out1[i].Abs());
    }

    iFFT(signal_out1, signal1, log2(length_signal));

//    for(int i=0;i<length_signal;i++){
//        printf("%f+%f*j\n", signal1[i].Rez(), signal1[i].Imz());
//    }

//    Integer a{'1','2','3','4','5','6','7','8'};
//    Integer b{'1','2','3','4','5','6','7','8'};
    Integer a{'9','2','3','4'};
    Integer b{'9','2','3','4'};
    Integer res;

    multiply(res,a,b);

    char* numStr1 = "222222222222222222222233";
    char* numStr2 = "23333333333333333333333333331";
//    int lenStr1,lenStr2,sumN;
//    lenStr1 = strlen(numStr1);
//    lenStr2 = strlen(numStr2);
//    sumN = lenStr1+lenStr2;
//    int N=1;
//    for(int i = 1;i < sumN;i*=2) {
//        N *=2;
//    }
//    auto* num1 = new Complex[N];
//    auto* num2 = new Complex[N];
//    auto* num3 = new Complex[N];
//    for (int i = lenStr1-1;i >= 0 ;i--) {
//        num1[lenStr1-i-1] = Complex((double)(numStr1[i]-'0'),0);
//    }
//    for (int i = lenStr2-1;i >= 0 ;i--) {
//        num2[lenStr2-i-1] = Complex((double)(numStr2[i]-'0'),0);
//    }
//    FFT(num1,num1,N);
//    FFT(num2,num2,N);
//    for(int i = 0;i < N;i++) {
//        num3[i] = num1[i]*num2[i];
//    }
//    iFFT(num3,num3,N);
//    vector<int> rst;
//    int curr,roundup=0;
//    for(int i = 0;i < sumN;i++) {
//        curr = num3[i].Rez()+roundup;
//        rst.push_back(curr%10);
//        roundup = curr/10;
//    }
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
//    cout<<endl<<"两个乘数： "<<numStr1<<","<<numStr2<<" 运算结果:  "<<result;


    return 0;
}
