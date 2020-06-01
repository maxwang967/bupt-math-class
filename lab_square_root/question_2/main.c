#include <math.h>
#include <stdio.h>
#include <time.h>

unsigned get_s(int start, int end, unsigned c)
{
    int mid = (start + end) / 2;
    long long a = 1 << (2 * mid);

    if (a >= c)
    {
        if ((a >> 2) < c)
            return mid;
        return get_s(start, mid, c);
    }
    else
    {
        if (mid == start)
        {
            return end;
        }
        return get_s(mid, start, c);
    }
}

unsigned my_isqrt(unsigned c)
{
    if (c <= 1)
        return c;
    unsigned c1 = c - 1;
    int s = get_s(0, 31, c);
    int x0 = 1 << s;
    int x1 = (x0 + (c >> s)) >> 1;
    while (x1 < x0)
    {
        x0 = x1;
        x1 = (x0 + c / x0) >> 1;
    }
    return x0;
}

int isqrt2(unsigned x)
{
    unsigned a = 1;
    unsigned b = (x >> 5) + 8;
    if (b > 65535)
        b = 65535; //  a <= sqrt(x) <= b
    do
    {
        int m = (a + b) >> 1;
        if (m * m > x)
            b = m - 1;
        else
            a = m + 1;
    } while (b >= a);
    return a - 1;
}

int isqrt3(unsigned x)
{
    if (x <= 1)
        return x;
    int x1 = x - 1;
    int s = 1;
    if (x1 > 65535)
    {
        s += 8;
        x1 >>= 16;
    }
    if (x1 > 255)
    {
        s += 4;
        x1 >>= 8;
    }
    if (x1 > 15)
    {
        s += 2;
        x1 >>= 4;
    }
    if (x1 > 3)
    {
        s += 1;
    }
    int x0 = 1 << s;
    x1 = (x0 + (x >> s)) >> 1;
    while (x1 < x0)
    {
        x0 = x1;
        x1 = (x0 + x / x0) >> 1;
    }
    return x0;
}

unsigned int isqrt4(unsigned long M)
{
    unsigned int N, i;
    unsigned long tmp, ttp;
    if (M == 0)
        return 0;
    N = 0;
    tmp = (M >> 30);
    M <<= 2;
    if (tmp > 1)
    {
        N++;
        tmp -= N;
    }
    for (i = 15; i > 0; i--)
    {
        N <<= 1;
        tmp <<= 2;
        tmp += (M >> 30);
        ttp = N;
        ttp = (ttp << 1) + 1;
        M <<= 2;
        if (tmp >= ttp)
        {
            tmp -= ttp;
            N++;
        }
    }
    return N;
}

int main()
{
    clock_t start, end;
    double err = 0;
    double max_err = 0;
    double rel_err = 0;
    double max_rel_err = 0;
    start = clock();
    for (int i = 0; i <= __INT_MAX__; i++)
    {
        double r = my_isqrt(i);
        double e = fabs(i - );
        double r_e = e / i;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / num_of_sample;
    rel_err = rel_err / num_of_sample;
    end = clock();
    printf("time1=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC);
    printf("err1=%.20f\n", err);
    printf("rel_err1=%.20f\n", rel_err);
    printf("max_err1=%.20f\n", max_err);
    printf("max_rel_err1=%.20f\n", max_rel_err);
}