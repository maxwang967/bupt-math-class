# -*- coding: utf-8 -*-
# @Time    : 2020/5/24 11:32 上午
# @Author  : morningstarwang
# @FileName: test.py
# @Blog: wangchenxing.com
from copy import deepcopy

import numpy as np
import global_for_eng as glb

from eng_value import gauss_hessen, transpose, dot, norm


def gram_schmidt(A):
    """Gram-schmidt正交化"""
    Q = np.zeros_like(A)
    cnt = 0
    for a in A.T:
        u = np.copy(a)
        for i in range(0, cnt):
            u -= np.dot(np.dot(Q[:, i].T, a), Q[:, i])  # 减去待求向量在以求向量上的投影
        if np.linalg.norm(u) != 0:
            e = u / np.linalg.norm(u)
        else:
            e = u
        Q[:, cnt] = e
        cnt += 1
    R = np.dot(Q.T, A)
    return (Q, R)


def givens_rotation(A):
    """Givens变换"""
    (r, c) = np.shape(A)
    Q = np.identity(r)
    R = np.copy(A)
    (rows, cols) = np.tril_indices(r, -1, c)
    print("----")
    print(rows, cols)
    print("----")
    for (row, col) in zip(rows, cols):
        if R[row, col] != 0:  # R[row, col]=0则c=1,s=0,R、Q不变
            r_ = np.hypot(R[col, col], R[row, col])  # d
            c = R[col, col] / r_
            s = -R[row, col] / r_
            G = np.identity(r)
            G[[col, row], [col, row]] = c
            G[row, col] = s
            G[col, row] = -s
            R = np.dot(G, R)  # R=G(n-1,n)*...*G(2n)*...*G(23,1n)*...*G(12)*A
            Q = np.dot(Q, G.T)  # Q=G(n-1,n).T*...*G(2n).T*...*G(23,1n).T*...*G(12).T
    return (Q, R)


def householder_reflection(A):
    """Householder变换"""
    (r, c) = np.shape(A)
    Q = np.identity(r)
    R = np.copy(A)
    for cnt in range(r - 1):
        x = R[cnt:, cnt]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        Q_cnt = np.identity(r)
        Q_cnt[cnt:, cnt:] -= 2.0 * np.outer(v, v)
        R = np.dot(Q_cnt, R)  # R=H(n-1)*...*H(2)*H(1)*A
        Q = np.dot(Q, Q_cnt)  # Q=H(n-1)*...*H(2)*H(1)  H为自逆矩阵
    return (Q, R)


# np.set_printoptions(precision=4, suppress=True)
# A = np.array([[6, 5, 0], [5, -1, 4], [5, 1, -14], [0, 4, 3]], dtype=float)

A = gauss_hessen(
    [
        [611, 196, -192, 407, -8, -52, -49, 29],
        [196, 899, 113, -192, -71, -43, -8, -44],
        [-192, 113, 899, 196, 61, 49, 8, 52],
        [407, -192, 196, 611, 8, 44, 59, -23],
        [-8, -71, 61, 8, 411, -599, 208, 208],
        [-52, -43, 49, 44, -599, 411, 208, 208],
        [-49, -8, 8, 59, 208, 208, 99, -911],
        [29, -44, 52, -23, 208, 208, -911, 99]
    ], 8
)
A = np.array(A)


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
    q = [
        [
            0 for _ in range(len(aa[0]))
        ] for _ in range(len(aa))
    ]
    cnt = 0
    a_t = transpose(aa)
    for a in a_t:
        u = deepcopy(a)
        for i in range(cnt):
            u_minus = dot(
                dot([
                    [x for x in q[:][i]]
                ],
                    [[x] for x in a]),
                [[x for x in q[:][i]]]
            )
            u = [x - u_minus[0][idx] for idx, x in enumerate(u)]
        e = list(map(lambda x: x / norm(u), u))
        for idx in range(len(q)):
            q[idx][cnt] = e[idx]
        cnt += 1
    r = dot(transpose(q), aa)
    return q, r


import copy
import math


class QR(object):

    def __init__(self, data):
        self.M = data
        self.degree = len(data)

    def get_row(self, index):
        res = []
        for i in range(self.degree):
            res.append(self.M[i][index])
        return res

    def get_col(self, index):
        res = []
        for i in range(self.degree):
            res.append(self.M[i][index])
        return res

    @staticmethod
    def dot(m1, m2):
        res = 0
        for i in range(len(m1)):
            res += m1[i] * m2[i]
        return res

    @staticmethod
    def list_multi(k, lt):
        res = []
        for i in range(len(lt)):
            res.append(k * lt[i])
        return res

    @staticmethod
    def one_item(x, yArr):
        res = [0 for i in range(len(x))]
        temp_y_arr = []

        n = len(yArr)
        if n == 0:
            res = x
        else:
            for item in yArr:
                k = QR.dot(x, item) / QR.dot(item, item)
                temp_y_arr.append(QR.list_multi(-k, item))
            temp_y_arr.append(x)

            for item in temp_y_arr:
                for i in range(len(item)):
                    res[i] += item[i]
        return res

    @staticmethod
    def normal(matrix):
        yArr = []
        yArr.append(matrix[0])

        for i in range(1, len(matrix)):
            yArr.append(QR.one_item(matrix[i], yArr))
        return yArr

    @staticmethod
    def normalized(lt):
        res = []
        sm = 0
        for item in lt:
            sm += math.pow(item, 2)
        sm = math.sqrt(sm)
        for item in lt:
            res.append(item / sm)
        return res

    @staticmethod
    def matrix_T(matrix):
        mat = copy.deepcopy(matrix)
        m = len(mat[0])
        n = len(mat)
        for i in range(m):
            for j in range(n):
                if i < j:
                    temp = mat[i][j]
                    mat[i][j] = mat[j][i]
                    mat[j][i] = temp
        return mat

    @staticmethod
    def matrix_multi(mat1, mat2):
        res = []
        rows = len(mat1[0])
        cols = len(mat1)
        for i in range(rows):
            temp = [0 for i in range(cols)]
            res.append(temp)

        for i in range(rows):
            for j in range(cols):
                sm = 0
                for k in range(cols):
                    sm += (mat1[k][i] * mat2[j][k])
                res[j][i] = sm
        return res

    def execute(self):
        xArr = []
        for i in range(self.degree):
            xArr.append(self.get_col(i))
        yArr = QR.normal(xArr)
        self.Q = []
        for item in yArr:
            self.Q.append(QR.normalized(item))

        self.R = QR.matrix_multi(QR.matrix_T(self.Q), xArr)
        return (self.Q, self.R)


# A = [
#     [1, 0, -1, 2, 1],
#     [3, 2, -3, 5, -3],
#     [2, 2, 1, 4, -2],
#     [0, 4, 3, 3, 1],
#     [1, 0, 8, -11, 4]
# ]
# A = [
#     [1, 2, 2],
#     [2, 1, 2],
#     [2, 2, 1]
# ]
a = [[0 for _ in range(11)] for _ in range(11)]
for i in range(1, 11):
    a[i][i - 1] = 1
for i in range(11):
    a[i][10] = -1

# temp = copy.deepcopy(gauss_hessen(A, 8))
# val = []  # 特征值
# times = 20  # 迭代次数
# for i in range(times):
#     qr = QR(temp)
#     (q, r) = qr.execute()
#     temp = QR.matrix_multi(r, q)
#     temp = QR.matrix_T(temp)
#
# for i in range(len(temp)):
#     for j in range(len(temp[0])):
#         if i == j:
#             val.append(temp[i][j])
# # 特征值
# print(val)

# glb.a = A
# (Q, R) = gram_schmidt(glb.a[2:6][2:6])
# print(Q)
# print(R)
# # print(np.dot(Q, R))
# print("---------")
#
# glb.a = A
# q, r = get_q_r(2, 6)
# print(np.array(q))
# print(np.array(r))
#
# (Q, R) = givens_rotation(A)
# print(Q)
# print(R)
# print(np.dot(Q, R))
#
(Q, R) = householder_reflection(a)
print(Q)
print(R)
print(np.dot(Q, R))
