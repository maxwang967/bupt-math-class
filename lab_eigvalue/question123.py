# -*- coding: utf-8 -*-
# @Time    : 2020/5/22 8:32 下午
# @Author  : morningstarwang
# @FileName: question123.py
# @Blog: wangchenxing.com

from copy import deepcopy
from datetime import datetime
from math import sqrt, sin, pi, atan, cos
import global_for_eng as glb
import numpy as np

# TODO this 'numpy' is to remove after release, and all codes using this package can be removed at same time
from numpy import linalg as la


def multiply_with_number(a, number):
    """
    multiply_with_number
    :param a: matrix
    :param number: a constant
    :return:
    """
    aa = deepcopy(a)
    aa = [
        [
            y * number for y in x
        ] for x in aa
    ]
    return aa


def divide_by_number(a, number):
    """
    a divided by number
    :param a: matrix
    :param number: a constant
    :return:
    """
    aa = deepcopy(a)
    aa = [
        [
            y / number for y in x
        ] for x in aa
    ]
    return aa


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


def create_identify_matrix(n):
    """
    create identify matrix (n x n)
    :param n: n dim
    :return: shape=n,n
    """
    return [
        [
            1 if i == j else 0 for j in range(n)
        ] for i in range(n)
    ]


def norm(v):
    sum_num = 0
    for vv in v:
        sum_num += vv * vv
    return sum_num ** 0.5


def lu_decomposition(aa):
    """
    lu decomposition
    :param a: coefficient matrix
    :return: matrix l and matrix u
    """
    a = deepcopy(aa)
    u = a
    l = create_identify_matrix(len(a))
    for j in range(len(a)):
        for i in range(j + 1, len(a)):
            mult = a[i][j] / a[j][j]
            for k in range(len(a[i])):
                u[i][k] = a[i][k] - a[j][k] * mult
            l[i][j] = mult
    return l, u


def compute_lu_result(a, bb, is_initial=False):
    """
    compute lu results for equations
    :param a: coefficient matrix
    :param b: right side matrix
    :return: equation results
    """
    bb = deepcopy(bb)
    if is_initial:
        b = [x[0] for x in bb]
    else:
        b = bb
    l, u = lu_decomposition(a)
    u_multi_x = [0 for _ in range(len(a))]
    for i in range(len(a)):
        right = b[i]
        for j in range(i):
            right -= l[i][j] * u_multi_x[j]
        u_multi_x[i] = right / l[i][i]
    results = [0 for _ in range(len(a[0]))]
    for i in range(len(u) - 1, -1, -1):
        right = u_multi_x[i]
        for j in range(len(u) - 1, i, -1):
            right -= u[i][j] * results[j]
        results[i] = right / u[i][i]
    return results


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


def inv_power_eng(a, n):
    """
    反幂法求矩阵最小特征值
    :param a: 待求矩阵
    :param n: 矩阵大小
    :return: 最小特征值，特征向量
    """
    iter_num = 1000
    # 全1矩阵，shape=(m, 1)
    v0 = [
        [1] for _ in range(n)
    ]
    # env = dot(a, v0)
    env = compute_lu_result(a, v0, is_initial=True)
    k = 0
    pld = 0
    # while abs(max(v0)[0] - max(v1)[0]) > delta:
    while k <= iter_num:
        v0 = env
        # v0_max = max(v0)[0]
        env = compute_lu_result(a, v0)
        # env = dot(a, v0)
        pld = max(map(abs, env))
        env = [[x] for x in env]
        env = divide_by_number(env, pld)
        env = [x[0] for x in env]
        # v1 = dot(a, v0)
        k += 1
        # print(f"delta={abs(max(v0)[0] - max(v1)[0])}")
    # eig = max(v1[0])
    # v = divide_by_number(v1, eig)
    return 1 / pld, env, True


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


def gauss_hessen(aa, n):
    s = [0.0 for _ in range(n)]
    u = [0.0 for _ in range(n)]
    c = 0
    I = [[0.0 for _ in range(n)] for _ in range(n)]
    for r in range(n):
        I[r][r] = 1.0
    H = [[0.0 for _ in range(n)] for _ in range(n)]
    Q = I
    B = I
    for r in range(n - 2):
        s = []
        for i in range(n):
            s.append(aa[i][r])
        e = [0.0 for _ in range(n)]
        for t in range(r + 1):
            s[t] = 0
        all_zero = 1
        for t in range(n):
            if s[t] != 0:
                all_zero = 0
                break
        if all_zero:
            for i in range(n):
                for j in range(n):
                    H[i][j] = I[i][j]
        else:
            s_m = 0
            for i in range(n):
                s_m = s_m + s[i] ** 2
            if aa[r + 1][r] == 0:
                c = sqrt(s_m)
            else:
                if aa[r + 1][r] >= 0:
                    c = -sqrt(s_m)
                else:
                    c = sqrt(s_m)
                e[r + 1] = 1
                for i in range(n):
                    u[i] = s[i] - c * e[i]
                uu = 0
                for i in range(n):
                    uu = uu + u[i] ** 2
                for i in range(n):
                    for j in range(n):
                        H[i][j] = I[i][j] - (2 * u[i] * u[j]) / uu
        B = dot(dot(H, aa), H)
        aa = B
        Q = dot(Q, H)
    return B


def get_q_r_householder(i, j):
    """
    HouseHolder方法获取Q、R矩阵
    :param i
    :param j
    :return: Q、R矩阵
    """
    aa = [
        [0 for _ in range(j - i)] for _ in range(j - i)
    ]
    r = len(aa)
    c = len(aa[0])
    for idx1 in range(j - i):
        for idx2 in range(j - i):
            aa[idx1][idx2] = glb.a[i + idx1][i + idx2]
    Q = create_identify_matrix(len(aa))
    R = deepcopy(aa)
    for cnt in range(r - 1):
        x = [R[x][cnt] for x in range(cnt, len(R))]
        e = [0.0 for x in range(len(x))]
        e[0] = norm(x)
        u = [x[i_xe] - e[i_xe] for i_xe in range(len(x))]
        norm_u = norm(u)
        # if norm_u == 0:
        #     continue
        v = [u[u_idx] / norm_u for u_idx in range(len(u))]
        Q_cnt = create_identify_matrix(r)
        outer = [
           [
               v[v_idx1] * v[v_idx2] for v_idx2 in range(len(v))
           ] for v_idx1 in range(len(v))
        ]
        for idx1 in range(cnt, len(Q_cnt)):
            for idx2 in range(cnt, len(Q_cnt[0])):
                Q_cnt[idx1][idx2] -= 2.0 * outer[idx1 - cnt][idx2 - cnt]
        R = dot(Q_cnt, R)
        Q = dot(Q, Q_cnt)
    return Q, R


def get_q_r(i, j):
    """
    Gram-schmidt正交化方法获取Q、R矩阵
    :param i
    :param j
    :return: Q、R矩阵
    """
    aa = [
        [0 for _ in range(j - i)] for _ in range(j - i)
    ]

    for idx1 in range(j - i):
        for idx2 in range(j - i):
            aa[idx1][idx2] = glb.a[i + idx1][i + idx2]
            # aa[idx1][idx2] = test[i + idx1][i + idx2]
    q = [
        [
            0.0 for _ in range(len(aa[0]))
        ] for _ in range(len(aa))
    ]
    r = [
        [
            0.0 for _ in range(len(aa[0]))
        ] for _ in range(len(aa))
    ]
    cnt = 0
    a_t = transpose(aa)
    for a in a_t:
        u = deepcopy(a)  # u=q
        for i in range(cnt):
            u_minus = 0
            for j in range(len(a)):
                u_minus = u_minus + q[j][i] * a[j]
            # u_minus = dot([[x for x in q[:][i]]], [[x] for x in a])
            r[i][cnt] = u_minus
            for j in range(len(u)):
                u[j] = u[j] - q[j][i] * u_minus
        q_norm = 0
        for j in range(len(u)):
            q_norm = q_norm + u[j] ** 2
        q_norm = q_norm ** 0.5
        r[cnt][cnt] = q_norm
        for j in range(len(u)):
            q[j][cnt] = u[j] / q_norm
        cnt += 1
    return q, r


def qr_aux(i, j):
    # 非二维或一维矩阵
    if j - i > 2:
        has_zero = False
        for idx1 in range(i + 1, j):
            if abs(glb.a[idx1][idx1 - 1]) < 0.0000000000001:
                has_zero = True
                qr_aux(i, idx1)
                qr_aux(idx1, j)
                break
        if not has_zero:
            # QR分解
            for _ in range(10000):
                q, r = get_q_r(i, j)
                rq = dot(r, q)
                for idx1 in range(i, j):
                    for idx2 in range(i, j):
                        glb.a[idx1][idx2] = rq[idx1 - i][idx2 - i]
                has_zero = False
                for idx1 in range(i + 1, j):
                    if abs(glb.a[idx1][idx1 - 1]) < 0.00000000001:
                        has_zero = True
                        if idx1 == i:
                            qr_aux(i, i + 1)
                            qr_aux(i + 1, j)
                        else:
                            qr_aux(i, idx1)
                            qr_aux(idx1, j)
                        break
                if has_zero:
                    count = 0
                    for i1 in range(j):
                        for i2 in range(j):
                            if i1 > i2:
                                count = count + abs(glb.a[i1][i2])
                    if count > 0:
                        if j - i < glb.n:
                            return True
                        else:
                            for i1 in range(glb.n):
                                if i1 not in eng_got:
                                    glb.en.append(glb.a[i1][i1])
                            return True
    # 是二维矩阵
    elif j - i == 2:
        # 手动算特征值
        glb.en.append((glb.a[i][i] + glb.a[i + 1][i + 1] + (
                (glb.a[i][i] - glb.a[i + 1][i + 1]) ** 2 + 4 * glb.a[i + 1][i] * glb.a[i][i + 1]) ** 0.5) / 2)
        glb.en.append((glb.a[i][i] + glb.a[i + 1][i + 1] - (
                (glb.a[i][i] - glb.a[i + 1][i + 1]) ** 2 + 4 * glb.a[i + 1][i] * glb.a[i][i + 1]) ** 0.5) / 2)
        eng_got.append(i)
        eng_got.append(i + 1)
        return True
    # 是一维矩阵
    elif j - i == 1:
        # 手动算特征值
        glb.en.append(glb.a[i][i])
        eng_got.append(i)
        return True
    else:
        return True


def qr_eng(h, m):
    glb.a = deepcopy(h)
    # print(glb.a)
    glb.n = m
    return qr_aux(0, glb.n)


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
        # if "E" == key:
        #     continue
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
        print(f"Matrix {key}:")
        start_time = datetime.now()
        pld, env, _ = inv_power_eng(all_data[key], len(all_data[key]))
        end_time = datetime.now()
        lambda1, v = la.eig(all_data[key])
        print(f"Running Time: {(end_time - start_time).microseconds / 1000}ms")
        print(f"my_lambda1={pld}")
        print(f"my_v={[x for x in env]}")
        env = [[x] for x in env]
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
        print(f"ground_truth_lambda1={lambda1}")
        print(f"differ={differ}")
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
    eng_got = []
    for key in all_data.keys():
        print(f"Matrix {key}:")
        start_time = datetime.now()
        h = gauss_hessen(all_data[key], len(all_data[key]))
        h_test = qr_eng(h, len(h))
        glb.en.sort(key=abs, reverse=True)
        end_time = datetime.now()
        print(f"Running Time: {(end_time - start_time).microseconds / 1000}ms")
        print(f"lambdas={glb.en}")
        ground_truth_lambdas, _ = la.eig(all_data[key])
        print(f"ground_truth_lambdas={ground_truth_lambdas}")
        # print(f"ground_truth_differ={abs(max(glb.en) - max(ground_truth_lambdas))}")
        print()
        glb.en.clear()
        eng_got.clear()

