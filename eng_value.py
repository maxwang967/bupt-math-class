# -*- coding: utf-8 -*-
# @Time    : 2020/5/22 8:32 下午
# @Author  : morningstarwang
# @FileName: eng_value.py
# @Blog: wangchenxing.com
from math import sqrt, sin, pi

from numpy import linalg as la


def divide_by_number(a, number):
    """
    a divided by number
    :param a: matrix
    :param number: a constant number
    :return:
    """
    a = [
        [
            y / number for y in x
        ] for x in a
    ]
    return a


def dot(a, b):
    """
    matrix multiplication
    :param a: shape=(m,n)
    :param b: shape=(n,p)
    :return: shape=(m,p)
    """
    assert len(a[0]) == len(b)
    m = len(a)
    n = len(a[0])
    p = len(b[0])
    # m,n * n,p => m,p
    result = [
        [
            0 for y in range(p)
        ] for x in range(m)
    ]
    for i in range(m):
        for j in range(p):
            for k in range(n):
                result[i][j] += a[i][k] * b[k][j]
    return result


def power_eng(a, iter_num):
    """
    幂法求矩阵最大特征值
    :param a: 待求矩阵
    :param iter_num: 迭代次数
    :return: 最大特征值，特征向量
    """
    m = len(a)
    # 全1矩阵，shape=(m, 1)
    v0 = [
        [1] for _ in range(m)
    ]
    v1 = dot(a, v0)
    k = 0
    m = 0
    # while abs(max(v0)[0] - max(v1)[0]) > delta:
    while k <= iter_num:
        v0 = v1
        # v0_max = max(v0)[0]
        v1 = dot(a, v0)
        m = list(map(abs, max(v1)))[0]
        v1 = divide_by_number(v1, m)
        # v1 = dot(a, v0)
        k += 1
        # print(f"delta={abs(max(v0)[0] - max(v1)[0])}")
    # eig = max(v1[0])
    # v = divide_by_number(v1, eig)
    return m, v1


if __name__ == '__main__':
    A = [
        [611, 196, -192, 407, -8, -52, -49, 29],
        [196, 899, 113, -192, -71, -43, -8, -44],
        [-192, 113, 899, 196, 61, 49, 8, 52],
        [407, -192, 196, 611, 8, 44, 59, -23],
        [-8, -71, 61, 8, 411, -599, 208, 208],
        [-52, -43, 49, 44, -599, 411, 208, 208],
        [-49, -8, 8, 59, 208, 208, 99, -911],
        [29, -44, 52, -23, 208, 208, -911, 99]
    ]
    B = [
        [1, 1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10],
        [1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11],
        [1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12],
        [1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13],
        [1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14],
        [1 / 6, 1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14, 1 / 15],
        [1 / 7, 1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14, 1 / 15, 1 / 16],
        [1 / 8, 1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14, 1 / 15, 1 / 16, 1 / 17],
        [1 / 9, 1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14, 1 / 15, 1 / 16, 1 / 17, 1 / 18],
        [1 / 10, 1 / 11, 1 / 12, 1 / 13, 1 / 14, 1 / 15, 1 / 16, 1 / 17, 1 / 18, 1 / 19]
    ]
    C = [
        [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        [11, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        [10, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        [9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        [8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1],
        [7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1],
        [6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1],
        [5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, 1],
        [4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 1],
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1],
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    ]
    D = [
        [
            (sqrt(2 / 21) * sin((j * k * pi) / 21)) for k in range(1, 21)
        ] for j in range(1, 21)
    ]
    all_data = {
        "A": A,
        "B": B,
        "C": C,
        "D": D
    }
    for key in all_data.keys():
        print(f"Matrix {key}:")
        my_lambda1, my_v = power_eng(all_data[key], 1000)
        lambda1, v = la.eig(all_data[key])
        print(f"my_lambda1={my_lambda1}")
        print(f"my_v={[x[0] for x in my_v]}")
        print(f"ground_truth_lambda1={max(lambda1)}")
        print(f"differ={abs(my_lambda1 - max(lambda1))}")

