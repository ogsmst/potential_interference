import numpy as np
import matplotlib as mp
from scipy import integrate
from functools import partial


# ここから各設定値
# 電流値
I = 100000

# 各電極のXY座標
er = np.array([[30, 30], [30, 35], [30, 40], [30, 45], [30, 50],
               [35, 30], [35, 35], [35, 40], [35, 45], [35, 50],
               [40, 30], [40, 35], [40, 40], [40, 45], [40, 50],
               [45, 30], [45, 35], [45, 40], [45, 45], [45, 50],
               [50, 30], [50, 35], [50, 40], [50, 45], [50, 50]])

# 大地抵抗率
ro = 200

# 電極長さ
# 構造体電極の長さ
length = np.full(er.shape[0], 5)

# 電極半径
# 各電極の半径を入力
radius = np.full(er.shape[0], 0.5)

# 埋設深さ
# 各電極の埋設深さを入力
depth = np.zeros(er.shape[0])

print(depth)

# 電極間距離のから行列を作成
S = np.zeros([np.shape(er)[0], np.shape(er)[0]])

# インピーダンス行列のから行列を作成
Z = np.zeros([np.shape(er)[0], np.shape(er)[0]])


# インピーダンス導出式
def impedance(z, i, j):
    return (1 / length[i]) * (ro / (4 * np.pi * length[j])) * (
            np.log((length[j] + depth[j] - z + np.sqrt(S[i, j] ** 2 + (length[j] + depth[j] - z) ** 2)) /
                   (depth[j] - z + np.sqrt(S[i, j] ** 2 + (depth[j] - z) ** 2))) +
            np.log((length[j] + depth[j] + z + np.sqrt(S[i, j] ** 2 + (length[j] + depth[j] + z) ** 2)) /
                   (depth[j] + z + np.sqrt(S[i, j] ** 2 + (depth[j] + z) ** 2)))
            )


# ここから処理

# 各電極距離の計算
def er_distance(i, j):
    for ii in i:
        for jj in j:
            if ii == jj:  # ii=jjの時電極距離は電極半径となる
                S[ii, jj] = radius[ii]
            else:
                S[ii, jj] = np.sqrt((er[ii, 0] - er[jj, 0]) ** 2 + (er[ii, 1] - er[jj, 1]) ** 2)


# インピーダンス行列を作成
def create_inpedance_matrix(i, j):
    for ii in i:
        for jj in j:
            matrix_impedance = partial(impedance, i=ii, j=jj)
            range0 = depth[ii]
            range1 = length[ii] + depth[ii]
            Z[ii, jj] = integrate.quad(matrix_impedance, range0, range1)[0]


# 構造体地下杭の数×構造体地下杭の数の行列をそれぞれ作成
# 地下杭同士の距離行列
er_distance(np.arange(np.shape(er)[0]), np.arange(np.shape(er)[0]))

# 地下杭のインピーダンス行列
create_inpedance_matrix(np.arange(np.shape(er)[0]), np.arange(np.shape(er)[0]))
print(Z)


# (x座標, y座標, 電極長さ, 電極半径, 埋設深さ)
info_b = [37, 27.5, 3, 0.005, 0.75]
info_d = [39, 27.5, 3, 0.005, 0.75]

# B種とD種の情報を入力、反映
er = np.append(er, np.array([info_d[0: 2], info_b[0: 2]]), axis=0)  # B, D種の(x, y)座標を挿入
print(er)
length = np.append(length, [info_d[2], info_b[2]])  # B, D種の電極長さを挿入
print(length)
radius = np.append(radius, [info_d[3], info_b[3]])  # B, D種の半径を挿入
print(radius)
depth = np.append(depth,  [info_d[4], info_b[4]])  # B, D種の半径を挿入


# インピーダンス行列の逆行列
Z_inv = np.linalg.inv(Z)

# 4つの小行列に分割
Q_11 = Z_inv[0: np.shape(er)[0] - 1, 0: np.shape(er)[0] - 1]
Q_12 = Z_inv[0: np.shape(er)[0] - 1, np.shape(er)[0] - 1]
Q_21 = Z_inv[np.shape(er)[0] - 1, 0: np.shape(er)[0] - 1]
Q_22 = Z_inv[np.shape(er)[0] - 1, np.shape(er)[0] - 1]

# それぞれの全要素の和を求める
Q_11_sum = np.sum(Q_11)
Q_12_sum = np.sum(Q_12)
Q_21_sum = np.sum(Q_21)
Q_22_sum = np.sum(Q_22)

U_A = Q_22_sum / (Q_22_sum * Q_11_sum - Q_12_sum * Q_21_sum) * I
U_B = - Q_21_sum / (Q_22_sum * Q_11_sum - Q_12_sum * Q_21_sum) * I

print(Z)
print(U_A, '[V]')
print(U_B, '[V]')
print('D種接地極と構造体杭との電位差は', U_A - U_B)
