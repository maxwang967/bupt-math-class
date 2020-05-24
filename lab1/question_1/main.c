#include <math.h>
#include <stdio.h>
#include <time.h>

float my_sqrt(float c)
{
    if (c < 0)
        return -1;
    double t = c;
    double rel_err = 1;
    while (fabs(rel_err) > 1e-15)
    {
        t = (c / t + t) / 2.0;
        rel_err = 1 - c / (t * t);
    }
    return t;
}

float my_sqrt_d(float c)
{
    if (c < 0)
        return -1;
    double t = c;
    double rel_err = 1;
    while (fabs(rel_err) > 1e-15)
    {
        t = (c / t + t) / 2.0;
        rel_err = 1 - c / (t * t);
    }
    return t;
}

float q_sqrt(float c)
{
    float c_half = 0.5f * c;
    int i = *(int *)&c;
    i = 0x5f375a86 - (i >> 1);
    c = *(float *)&i;
    c = c * (1.5f - c_half * c * c);
    c = c * (1.5f - c_half * c * c);
    return 1 / c;
}

int test_float()
{
    const int num_of_sample = 1000000;
    float *numbers = malloc(sizeof(float) * num_of_sample);
    for (int i = 0; i < num_of_sample; i++)
    {
        numbers[i] = (float)rand() / RAND_MAX * 1000000;
    }

    clock_t start, end;
    float err = 0;
    float max_err = 0;
    float rel_err = 0;
    float max_rel_err = 0;
    start = clock();
    for (int i = 0; i < num_of_sample; i++)
    {
        float r = my_sqrt(numbers[i]);
        float e = fabs(numbers[i] - r * r);
        float r_e = e / numbers[i];
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

    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (int i = 0; i < num_of_sample; i++)
    {
        float r = q_sqrt(numbers[i]);
        float e = fabs(numbers[i] - r * r);
        float r_e = e / numbers[i];
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / num_of_sample;
    rel_err = rel_err / num_of_sample;
    end = clock();
    printf("time2=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC);
    printf("err2=%.20f\n", err);
    printf("rel_err2=%.20f\n", rel_err);
    printf("max_err2=%.20f\n", max_err);
    printf("max_rel_err2=%.20f\n", max_rel_err);

    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (int i = 0; i < num_of_sample; i++)
    {
        float r = sqrt(numbers[i]);
        float e = fabs(numbers[i] - r * r);
        float r_e = e / numbers[i];
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / num_of_sample;
    rel_err = rel_err / num_of_sample;
    end = clock();
    printf("time3=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC);
    printf("err3=%.20f\n", err);
    printf("rel_err3=%.20f\n", rel_err);
    printf("max_err3=%.20f\n", max_err);
    printf("max_rel_err3=%.20f\n", max_rel_err);

    free(numbers);
}

int test_double()
{
    const int num_of_sample = 1000000;
    double *numbers = malloc(sizeof(double) * num_of_sample);
    for (int i = 0; i < num_of_sample; i++)
    {
        numbers[i] = (double)rand() / RAND_MAX * 1000000;
    }

    clock_t start, end;
    double err = 0;
    double max_err = 0;
    double rel_err = 0;
    double max_rel_err = 0;
    start = clock();
    for (int i = 0; i < num_of_sample; i++)
    {
        double r = my_sqrt(numbers[i]);
        double e = fabs(numbers[i] - r * r);
        double r_e = e / numbers[i];
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

    err = 0;
    max_err = 0;
    rel_err = 0;
    max_rel_err = 0;
    start = clock();
    for (int i = 0; i < num_of_sample; i++)
    {
        double r = my_sqrt(numbers[i]);
        double e = fabs(numbers[i] - r * r);
        double r_e = e / numbers[i];
        max_err = fmax(max_err, e);
        max_rel_err = fmax(max_rel_err, r_e);
        err += e;
        rel_err += r_e;
    }
    err = err / num_of_sample;
    rel_err = rel_err / num_of_sample;
    end = clock();
    printf("time2=%.20f\n", (double)(end - start) * 1000 / CLOCKS_PER_SEC);
    printf("err2=%.20f\n", err);
    printf("rel_err2=%.20f\n", rel_err);
    printf("max_err2=%.20f\n", max_err);
    printf("max_rel_err2=%.20f\n", max_rel_err);
}

int main()
{
    test_double();
}