# -*- coding: utf-8 -*-
# @Time    : 2020/5/22 8:32 下午
# @Author  : morningstarwang
# @FileName: eng_value.py
# @Blog: wangchenxing.com
from copy import deepcopy
from datetime import datetime
from math import sqrt, sin, pi, atan, cos

from numpy import linalg as la


def multiply_with_number(a, number):
    """
    multiply_with_number
    :param a: matrix
    :param number: a constant
    :return:
    """
    a = [
        [
            y * number for y in x
        ] for x in a
    ]
    return a


def divide_by_number(a, number):
    """
    a divided by number
    :param a: matrix
    :param number: a constant
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
            0 for _ in range(p)
        ] for _ in range(m)
    ]
    for i in range(m):
        for j in range(p):
            for k in range(n):
                result[i][j] += a[i][k] * b[k][j]
    return result


def transpose(a):
    """
    transpose matrix
    :param a: matrix
    :return: transposed matrix
    """
    m = len(a)
    n = len(a[0])
    result = [
        [
            0 for _ in range(m)
        ] for _ in range(n)
    ]
    for i in range(n):
        for j in range(m):
            result[i][j] = a[j][i]
    return result


def power_eng(a, n):
    """
    幂法求矩阵最大特征值
    :param a: 待求矩阵
    :param n: 矩阵大小
    :return: 最大特征值，特征向量
    """
    iter_num = 1000
    # 全1矩阵，shape=(m, 1)
    v0 = [
        [1] for _ in range(n)
    ]
    env = dot(a, v0)
    k = 0
    pld = 0
    # while abs(max(v0)[0] - max(v1)[0]) > delta:
    while k <= iter_num:
        v0 = env
        # v0_max = max(v0)[0]
        env = dot(a, v0)
        pld = list(map(abs, max(env)))[0]
        env = divide_by_number(env, pld)
        # v1 = dot(a, v0)
        k += 1
        # print(f"delta={abs(max(v0)[0] - max(v1)[0])}")
    # eig = max(v1[0])
    # v = divide_by_number(v1, eig)
    return pld, env, True


def jacobi_eng(a, n):
    """
    Jacobi迭代法求矩阵特征值
    :param a: 待求矩阵
    :param n: 矩阵大小
    :return: 对角矩阵，对角线上为特征值
    """
    assert len(a) == len(a[0])
    iter_num = 0
    while iter_num < 1000:
        a_t = deepcopy(a)
        a_max = 0
        j = -1
        k = -1
        for i1 in range(n):
            for i2 in range(n):
                # 去除下三角和主对角线
                if i1 >= i2:
                    continue
                if abs(a[i1][i2]) > a_max:
                    a_max = abs(a[i1][i2])
                    j = i1
                    k = i2
        a_jk = a[j][k]
        if a_jk == 0:
            break
        if a[j][j] == a[k][k]:
            theta = pi / 4
        else:
            theta = atan(
                (2 * a[j][k]) / (a[j][j] - a[k][k])
            ) / 2
        s = sin(theta)
        s2 = sin(2 * theta)
        c = cos(theta)
        c2 = cos(theta * 2)
        # 旋转矩阵
        # 1 jj,kk
        a_t[j][j] = a[j][j] * c * c + a[k][k] * s * s + 2 * a[j][k] * s * c
        a_t[k][k] = a[j][j] * s * s + a[k][k] * c * c - 2 * a[j][k] * s * c
        # 2 jk,kj
        a_t[j][k] = 0.5 * (a[k][k] - a[j][j]) * s2 + a[j][k] * c2
        a_t[k][j] = 0.5 * (a[k][k] - a[j][j]) * s2 + a[j][k] * c2
        # 3 j*,*j
        # 4 k*,*k
        for i in range(len(a)):
            if i != j and i != k:
                a_t[j][i] = a[j][i] * c + a[k][i] * s
                a_t[i][j] = a[j][i] * c + a[k][i] * s
                a_t[k][i] = a[k][i] * c - a[j][i] * s
                a_t[i][k] = a[k][i] * c - a[j][i] * s
        # 5 **
        # DO NOTHING
        a = a_t
        del a_t
        iter_num += 1
    ev = [x[idx] for idx, x in enumerate(a)]
    return ev, True


def gauss_hessen(a, n):
    """
    将矩阵转换为上Hessenberg矩阵
    :param a: 待求矩阵
    :param n: 矩阵大小
    :return: 上Hessenberg矩阵
    """
    for m in range(n-1):
        x = 0.0
        i = m
        # 寻找主元
        for j in range(n):
            if abs(a[j][m-1]) > abs(x):
                x = a[j][m-1]
                i = j
        # 交换行列
        if i != m:
            for j in range(m-1, n):
                # SWAP(a[i][j], a[m][j])
                t = a[i][j]
                a[i][j] = a[m][j]
                a[m][j] = t
            for j in range(n):
                # SWAP(a[j][i], a[j][m])
                t = a[j][i]
                a[j][i] = a[j][m]
                a[j][m] = t
        # 执行消去法
        if x != 0.0:
            for i in range(m+1, n):
                y = a[i][m - 1]
                if y != 0:
                    y /= x
                    a[i][m-1] = y
                    for j in range(m, n):
                        a[i][j] -= y * a[m][j]
                    for j in range(0, n):
                        a[j][m] += y * a[j][i]
    # i > j + 1的元素输出为0
    for i in range(n):
        for j in range(n):
            if i > j + 1:
                a[i][j] = 0
    return a


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
    E = [
        [
            1 if i == j else (-1 if i > j else (1 if j == 49 else 0)) for j in range(50)
        ] for i in range(50)
    ]
    all_data = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E
    }

    # Power Method
    print("---------------------")
    print(f"Power Method:")
    for key in all_data.keys():
        if "E" == key:
            continue
        print(f"Matrix {key}:")
        start_time = datetime.now()
        pld, env, _ = power_eng(all_data[key], len(all_data[key]))
        end_time = datetime.now()
        lambda1, v = la.eig(all_data[key])
        print(f"Running Time: {(end_time - start_time).microseconds / 1000}ms")
        print(f"my_lambda1={pld}")
        print(f"my_v={[x[0] for x in env]}")
        Ax = dot(all_data[key], env)
        lambda_x = multiply_with_number(env, pld)
        differ_matrix = [
            [
                0 for _ in range(len(Ax))
            ] for _ in range(len(Ax))
        ]
        for i in range(len(Ax)):
            for j in range(1):
                differ_matrix[i][j] = Ax[i][j] - lambda_x[i][j]
        square_sum = 0
        for i in range(len(differ_matrix)):
            for j in range(len(differ_matrix)):
                square_sum += differ_matrix[i][j] * differ_matrix[i][j]
        differ = sqrt(square_sum)
        print(f"ground_truth_lambda1={max(lambda1)}")
        print(f"differ={differ}")
        print()

    # Inv Power Method
    print("---------------------")
    print(f"Inv Power Method:")
    for key in all_data.keys():
        if "E" == key:
            continue
        # TODO 线性方程组求解
        print()

    # Jacobi Method
    print("---------------------")
    print(f"Jacobi Method:")
    for key in all_data.keys():
        if "E" == key:
            continue
        print(f"Matrix {key}:")
        start_time = datetime.now()
        ev, _ = jacobi_eng(all_data[key], len(all_data[key]))
        end_time = datetime.now()

        print(f"Running Time: {(end_time - start_time).microseconds / 1000}ms")
        print(f"lambdas={ev}")
        ground_truth_lambdas, _ = la.eig(all_data[key])
        print(f"ground_truth_lambdas={ground_truth_lambdas}")
        print(f"ground_truth_differ={abs(max(ev) - max(ground_truth_lambdas))}")
        print()

    # QR Method
    print("---------------------")
    print(f"QR Method:")
    for key in all_data.keys():
        # if "E" == key:
        #     continue
        print(f"Matrix {key}:")
        start_time = datetime.now()
        # TODO QR算法
        ev, _ = jacobi_eng(all_data[key], len(all_data[key]))
        end_time = datetime.now()

        print(f"Running Time: {(end_time - start_time).microseconds / 1000}ms")
        print(f"lambdas={ev}")
        ground_truth_lambdas, _ = la.eig(all_data[key])
        print(f"ground_truth_lambdas={ground_truth_lambdas}")
        print(f"ground_truth_differ={abs(max(ev) - max(ground_truth_lambdas))}")
        print()
