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

# 電極間距離のから行列を作成（+2はB,D種の二つ分を追加）
S = np.zeros([np.shape(er)[0] + 2, np.shape(er)[0] + 2])

# インピーダンス行列のから行列を作成（+2はB,D種の二つ分を追加）
Z = np.zeros([np.shape(er)[0] + 2, np.shape(er)[0] + 2])


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
for ii in range(np.shape(er)[0]):
    for jj in range(np.shape(er)[0]):
        if ii == jj:  # ii=jjの時電極距離は電極半径となる
            S[ii, jj] = radius[ii]
        else:
            S[ii, jj] = np.sqrt((er[ii, 0] - er[jj, 0]) ** 2 + (er[ii, 1] - er[jj, 1]) ** 2)

# インピーダンス行列を作成
for ii in range(np.shape(er)[0]):
    for jj in range(np.shape(er)[0]):
        matrix_impedance = partial(impedance, i=ii, j=jj)
        range0 = depth[ii]
        range1 = length[ii] + depth[ii]
        Z[ii, jj] = integrate.quad(matrix_impedance, range0, range1)[0]


print(Z)

# D種接地極の座標のパターン
position_D = np.array([39, 27.5])
# B種接地極が動きうる座標のパターン
position_B_interval = 0.1  # B種接地極それぞれの間隔
position_B_interval_number = 201  # B種接地極が何パターンあるか（1軸のみで）0から
iteration_position_B = np.array([[30 + i * position_B_interval, 30 + j * position_B_interval]
                                for i in range(position_B_interval_number) for j in range(position_B_interval_number)])
# 結果の電圧を収める2次元のから行列を作成
V = np.zeros([position_B_interval_number, position_B_interval_number])

for ii in range(position_B_interval_number):
    for jj in range(position_B_interval_number):
        # (x座標, y座標, 電極長さ, 電極半径, 埋設深さ)を設定
        info_d = [position_D[0], position_D[1], 3, 0.005, 0.75]
        info_b = [iteration_position_B[position_B_interval_number * ii, 0], iteration_position_B[position_B_interval_number * ii + jj, 1], 3, 0.005, 0.75]

        # B種とD種の情報を反映
        er = np.append(er, np.array([info_d[0: 2], info_b[0: 2]]), axis=0)  # B, D種の(x, y)座標を挿入
        print(er)
        length = np.append(length, [info_d[2], info_b[2]])  # B, D種の電極長さを挿入
        print(length)
        radius = np.append(radius, [info_d[3], info_b[3]])  # B, D種の半径を挿入
        print(radius)
        depth = np.append(depth,  [info_d[4], info_b[4]])  # B, D種の半径を挿入

        # 各電極距離の計算（B,D種部分のみ追加計算を行う）
        for ii in range(np.shape(er)[0]):
            for jj in range(np.shape(er)[0]):
                if ii < np.shape(er)[0] - 2 and jj < np.shape(er)[0] - 2:  # 既に計算してある部分は省略
                    continue
                if ii == jj:  # ii=jjの時電極距離は電極半径となる
                    S[ii, jj] = radius[ii]
                else:
                    S[ii, jj] = np.sqrt((er[ii, 0] - er[jj, 0]) ** 2 + (er[ii, 1] - er[jj, 1]) ** 2)

        # 構造体地下杭とB種の座標が重なるときは計算を行わない処理
        # 各電極間距離で地下杭の半径以下をFalseそれ以上をTrueとする行列を作成
        S_condition = S > 1
        # B種と地下杭の関係を示す一番右端の列を探査し、全部True時のみ次の処理に進む
        if np.all(S_condition[0 : np.shape(S_condition)[0] - 2 , np.shape(S_condition)[1]]):
            continue

        # インピーダンス行列を作成（B,D種部分のみ追加計算を行う）
        for ii in range(np.shape(er)[0]):
            for jj in range(np.shape(er)[0]):
                if ii < np.shape(er)[0] - 2 and jj < np.shape(er)[0] - 2:  # 既に計算してある部分は省略
                    continue
                matrix_impedance = partial(impedance, i=ii, j=jj)
                range0 = depth[ii]
                range1 = length[ii] + depth[ii]
                Z[ii, jj] = integrate.quad(matrix_impedance, range0, range1)[0]

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
