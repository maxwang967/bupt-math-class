#include <math.h>
#include <stdio.h>
#include <time.h>

double n_loop = 0;

unsigned get_s(int start, int end, unsigned c)
{
    int mid = (start + end) / 2;
    long long a = 1 << (2 * mid);

    if (a >= c)
    {
        if ((a >> 2) < c)
            return mid;
        return get_s(start, mid - 1, c);
    }
    else
    {
        if (mid == start)
        {
            return end;
        }
        return get_s(mid + 1, end, c);
    }
}

unsigned my_isqrt(unsigned c)
{
    if (c <= 1)
        return c;
    int s = get_s(0, 16, c);
    int x0 = 1 << s;
    int x1 = (x0 + (c >> s)) >> 1;
    while (x1 < x0)
    {
        x0 = x1;
        x1 = (x0 + c / x0) >> 1;
        n_loop++;
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
        n_loop++;
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
    unsigned int x1 = x - 1;
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
    unsigned int x0 = 1 << s;
    x1 = (x0 + (x >> s)) >> 1;
    while (x1 < x0)
    {
        n_loop++;
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
        n_loop++;
        N <<= 1;
        tmp <<= 2;
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

void test_unsigned()
{
    clock_t start, end;
    double err = 0;
    double max_err = 0;
    double rel_err = 0;
    double max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i != 0; i++)
    {
        unsigned r = sqrt(i);
        // if (i % (__UINT32_MAX__ / 100) == 0)
        // {
        //     printf("%d\n", i / (__UINT32_MAX__ / 100));
        // }
    }
    end = clock();
    double time1 = (double)(end - start) * 1000 / CLOCKS_PER_SEC;
    printf("time1=%.20f\n", time1);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i != 0; i++)
    {
        unsigned r = my_isqrt(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
        if (i % (__UINT32_MAX__ / 10) == 0)
        {
            printf("%d\n", i / (__UINT32_MAX__ / 10));
        }
    }
    err = err / __UINT32_MAX__;
    rel_err = rel_err / __UINT32_MAX__;
    end = clock();
    printf("time2=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err2=%.20f\n", err);
    printf("rel_err2=%.20f\n", rel_err);
    printf("max_err2=%.20f\n", max_err);
    printf("max_rel_err2=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / __UINT32_MAX__);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i != 0; i++)
    {
        unsigned r = isqrt2(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
        if (i % (__UINT32_MAX__ / 10) == 0)
        {
            printf("%d\n", i / (__UINT32_MAX__ / 10));
        }
    }
    err = err / __UINT32_MAX__;
    rel_err = rel_err / __UINT32_MAX__;
    end = clock();
    printf("time3=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err3=%.20f\n", err);
    printf("rel_err3=%.20f\n", rel_err);
    printf("max_err3=%.20f\n", max_err);
    printf("max_rel_err3=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / __UINT32_MAX__);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i != 0; i++)
    {
        unsigned r = isqrt3(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
        if (i % (__UINT32_MAX__ / 10) == 0)
        {
            printf("%d\n", i / (__UINT32_MAX__ / 10));
        }
    }
    err = err / __UINT32_MAX__;
    rel_err = rel_err / __UINT32_MAX__;
    end = clock();
    printf("time4=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err4=%.20f\n", err);
    printf("rel_err4=%.20f\n", rel_err);
    printf("max_err4=%.20f\n", max_err);
    printf("max_rel_err4=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / __UINT32_MAX__);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i != 0; i++)
    {
        unsigned r = isqrt4(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
        if (i % (__UINT32_MAX__ / 10) == 0)
        {
            printf("%d\n", i / (__UINT32_MAX__ / 10));
        }
    }
    err = err / __UINT32_MAX__;
    rel_err = rel_err / __UINT32_MAX__;
    end = clock();
    printf("time5=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err5=%.20f\n", err);
    printf("rel_err5=%.20f\n", rel_err);
    printf("max_err5=%.20f\n", max_err);
    printf("max_rel_err5=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / __UINT32_MAX__);
}

unsigned my_isqrt_op(unsigned c)
{
    if (c <= 1)
        return c;
    unsigned a = 1;
    unsigned s = 0;
    // 根据权重得出的二叉搜索树
    if (c < a << 4)
    {
        if (c < a << 2)
        {
            if (c < a << 1)
                s = 1;
            else
                s = 2;
        }
        else
        {
            if (c < a << 3)
                s = 3;
            else
                s = 4;
        }
    }
    else
    {
        if (c < a << 6)
        {
            if (c < a << 5)
                s = 5;
            else
                s = 6;
        }
        else
        {
            if (c < a << 7)
                s = 7;
            else
                s = 8;
            if (c >= a << 8)
            {
                if (c < a << 12)
                {
                    if (c < a << 10)
                    {
                        if (c < a << 9)
                            s = 9;
                        else
                            s = 10;
                    }
                    else
                    {
                        if (c < a << 11)
                            s = 11;
                        else
                            s = 12;
                    }
                }
                else
                {
                    if (c < a << 14)
                    {
                        if (c < a << 13)
                            s = 13;
                        else
                            s = 14;
                    }
                    else
                    {
                        if (c < a << 15)
                            s = 15;
                        else
                            s = 16;
                    }
                }
            }
        }
    }
    int x0 = 1 << s;
    int x1 = (x0 + (c >> s)) >> 1;
    while (x1 < x0)
    {
        x0 = x1;
        x1 = (x0 + c / x0) >> 1;
        n_loop++;
    }
    return x0;
}

void test_op_l()
{
    clock_t start, end;
    double err = 0;
    double max_err = 0;
    double rel_err = 0;
    double max_rel_err = 0;
    const unsigned n = ((1 << 8) + 1) * ((1 << 8) + 1) - 1;
    start = clock();
    for (unsigned i = 1; i <= n; i++)
    {
        unsigned r = sqrt(i);
        // if (i % (__UINT32_MAX__ / 100) == 0)
        // {
        //     printf("%d\n", i / (__UINT32_MAX__ / 100));
        // }
    }
    end = clock();
    double time1 = (double)(end - start) * 1000 / CLOCKS_PER_SEC;
    printf("time1=%.20f\n", time1);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i <= n; i++)
    {
        unsigned r = my_isqrt_op(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time2=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err2=%.20f\n", err);
    printf("rel_err2=%.20f\n", rel_err);
    printf("max_err2=%.20f\n", max_err);
    printf("max_rel_err2=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i <= n; i++)
    {
        unsigned r = isqrt2(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time3=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err3=%.20f\n", err);
    printf("rel_err3=%.20f\n", rel_err);
    printf("max_err3=%.20f\n", max_err);
    printf("max_rel_err3=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i <= n; i++)
    {
        unsigned r = isqrt3(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time4=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err4=%.20f\n", err);
    printf("rel_err4=%.20f\n", rel_err);
    printf("max_err4=%.20f\n", max_err);
    printf("max_rel_err4=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = 1; i <= n; i++)
    {
        unsigned r = isqrt4(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time5=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err5=%.20f\n", err);
    printf("rel_err5=%.20f\n", rel_err);
    printf("max_err5=%.20f\n", max_err);
    printf("max_rel_err5=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);
}

void test_op_g()
{
    clock_t start, end;
    double err = 0;
    double max_err = 0;
    double rel_err = 0;
    double max_rel_err = 0;
    unsigned n = __UINT32_MAX__ - ((1 << 8) + 1) * ((1 << 8) + 1);
    start = clock();
    for (unsigned i = ((1 << 8) + 1) * ((1 << 8) + 1); i != 0; i++)
    {
        unsigned r = sqrt(i);
        // if (i % (__UINT32_MAX__ / 100) == 0)
        // {
        //     printf("%d\n", i / (__UINT32_MAX__ / 100));
        // }
    }
    end = clock();
    double time1 = (double)(end - start) * 1000 / CLOCKS_PER_SEC;
    printf("time1=%.20f\n", time1);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = ((1 << 8) + 1) * ((1 << 8) + 1); i != 0; i++)
    {
        unsigned r = my_isqrt_op(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time2=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err2=%.20f\n", err);
    printf("rel_err2=%.20f\n", rel_err);
    printf("max_err2=%.20f\n", max_err);
    printf("max_rel_err2=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = ((1 << 8) + 1) * ((1 << 8) + 1); i != 0; i++)
    {
        unsigned r = isqrt2(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time3=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err3=%.20f\n", err);
    printf("rel_err3=%.20f\n", rel_err);
    printf("max_err3=%.20f\n", max_err);
    printf("max_rel_err3=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = ((1 << 8) + 1) * ((1 << 8) + 1); i != 0; i++)
    {
        unsigned r = isqrt3(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time4=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err4=%.20f\n", err);
    printf("rel_err4=%.20f\n", rel_err);
    printf("max_err4=%.20f\n", max_err);
    printf("max_rel_err4=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);

    n_loop = 0;
    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (unsigned i = ((1 << 8) + 1) * ((1 << 8) + 1); i != 0; i++)
    {
        unsigned r = isqrt4(i);
        double t = sqrt(i);
        double e = fabs(r - t);
        double r_e = e / t;
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / n;
    rel_err = rel_err / n;
    end = clock();
    printf("time5=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC - time1);
    printf("err5=%.20f\n", err);
    printf("rel_err5=%.20f\n", rel_err);
    printf("max_err5=%.20f\n", max_err);
    printf("max_rel_err5=%.20f\n", max_rel_err);
    printf("n_loop=%.20f\n", n_loop / n);
}

int main()
{
    // test_unsigned();
    test_op_l();
    test_op_g();
}