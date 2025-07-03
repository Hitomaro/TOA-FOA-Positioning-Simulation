import spiceypy as spice
import numpy as np
import matplotlib
import os
import math

# 大域変数定義ここから
localapp = os.environ['LOCALAPPDATA']
# 各種定数を宣言
radius_moon_km = 1737.4                     # 月の半径(km)
radius_moon_m = radius_moon_km * 1000       # 月の半径(m)
c = 3e8                                     # 光速(m/s)
f_carrier_MHz = 2483.5                      # 搬送波周波数(MHz)
f_carrier = f_carrier_MHz * 10**6           # 搬送波周波数(Hz)
clock_drift = 1e-6                          # クロックドリフト
base_clock_bias = 0                         # 基準クロックバイアス
moon_flatting = 0                           # 月を起伏のない球面と仮定する
# 測位関連の技術的仕様を宣言
pos_interval = 10                           # 測位時間間隔(sec)
threshold = 1e-4                            # 十分に測位結果が収束したと判断する距離誤差(m)
max_iter = 10                               # 最小二乗法の最高反復処理回数
elev_mask_angle = 10                        # 最低測位可能仰角(仰角マスク角)
sat_num = 6                                 # 測位衛星数
sat_name = 'sat'                            # 衛星名
# ユーザー位置設定
lat = math.radians(-90)                     # 緯度
lon = math.radians(0)                       # 経度
alt = 0.0                                   # 高度
user_pos = spice.georec(lon, lat, alt, radius_moon_km, moon_flatting)
# 大域変数定義ここまで

def load_spices():
    base_dir = '.\constant_files'
    spice.furnsh(os.path.join(base_dir, 'naif0012.tls'))
    spice.furnsh(os.path.join(base_dir, 'moon_pa_de421_1900-2050.bpc'))
    spice.furnsh(os.path.join(base_dir, 'pck0010.tpc'))
    spice.furnsh(os.path.join(base_dir, 'moon_080317.tf'))
    spice.furnsh(os.path.join(base_dir, 'de421.bsp'))

def load_sat_ids(base_dir, name, num):
    sat_id_list = []
    for i in range(num):
        bsp_file = os.path.join(base_dir, name, num)
        spice.furnsh(bsp_file)
        sat_id_list(i) = spice.spkobj(bsp_file, 1)
    return sat_id_list

def get_sat_pos_vel(et, sat_ids, abcorr):
    '''
    引数
        et: エフェメリス秒
        sat_ids: 衛星ID配列
        abcorr: 補正方法 (例: 'NONE', 'LT+S')
    戻り値
        pos: etにおける衛星位置ベクトル [[x], [y], [z]] (m)
        vel: etにおける衛星速度ベクトル [[vx], [vy], [vz]] (m/s)
    '''

    sat_id = str(sat_ids)
    state_j2000, _ = spice.spkezr(sat_id, et, 'J2000', abcorr, 'MOON')

    rot_mat = spice.sxform('J2000', 'IAU_MOON', et)
    state = rot_mat * state_j2000
    
    pos = state[:3]
    vel = state[3:]
    
    return pos, vel

def elevation_masking(sat_pos, user_pos):
    '''
    引数
        sat_pos: [3x1] 衛星位置（x,y,z 月中心座標系）
        user_pos: [3x1] ユーザー位置 （x,y,z 月中心座標系）
    戻り値
        el: ユーザー視点の衛星仰角(°)
    '''
    los = sat_pos - user_pos
    n_user = user_pos / np.linalg.norm(user_pos)
    r_los = los / np.linalg.norm(los)
    cos_el = np.dot(n_user, r_los)
    el = np.degrees(np.arcsin(cos_el))
    return el

def calc_user_vel(user_pos):
    '''
    引数
        user_pos: [3x1] ユーザー位置(x, y, z: 月中心座標系)
    戻り値
        user_vel: [3x1] ユーザー速度ベクトル(vx, vy, vz)
    '''
    
    omega_moon = 2 * np.pi() / (27.32 * 24 * 60 * 60)
    omega_vec = [[0], [0], [omega_moon]]
    user_vel = np.cross(omega_vec, user_pos)
    return user_vel

def calc_TOA(user_data, sat_pos_in_gravity, sat_pos_ephemeris, clock_bias, noise):
    '''
    引数
        user_data: [1x5] 推定ユーザーデータ（[1~3]ユーザー位置, [4]クロックバイアス, [5]クロックドリフト量（TOAでは使わない）
        sat_pos_in_gravity: [3x1] 真の衛星位置(x, y, z: 月中心座標系)
        sat_pos_ephemeris: [3x1] 軌道データ上の衛星位置(x, y, z: 月中心座標系)
        clock_bias: 真のクロックバイアス
        noise: 距離誤差ノイズ(m)
    戻り値
        residual: 疑似距離 - 観測距離残差 (m)
        A: [1x4] ヤコビ行列
    '''
    # 推定情報の処理
    estimate_user_pos = user_data[:3]                                           # 推定ユーザー位置[x_user, y_user, z_user]
    estimate_clock_bias = user_data[4]                                          # 推定クロックバイアス
    estimate_range = np.linalg.norm(sat_pos_ephemeris - estimate_user_pos)      # 推定位置を用いた衛星-ユーザー間疑似距離
    estimate_pseudo_range = estimate_range + c * estimate_clock_bias            # 推定衛星-ユーザー間距離
    
    # 実情報の処理
    true_range = np.linalg.norm(sat_pos_in_gravity - user_pos)
    true_pseudo_range = true_range + c * clock_bias
    obs_range = true_pseudo_range + noise
    residual = estimate_pseudo_range - obs_range

    # ヤコビ行列Aの作成
    A = []
    A[:3] = np.transpose(estimate_user_pos - sat_pos_in_gravity) / estimate_range
    A[4] = c
    
    return residual, A

def calc_FOA(user_data, sat_pos_in_gravity, sat_vel_in_gravity, sat_pos_ephemeris, sat_vel_ephemeris, user_pos, user_vel):
    '''
    引数
        user_data: [1x5] 推定ユーザーデータ([1~3]ユーザー位置(x, y, z: 月中心座標系), [4]クロックバイアス(FOAでは使わない), [5]クロックドリフト)
        sat_pos_in_gravity: [3x1] 真の衛星位置(x, y, z: 月中心座標系)
        sat_vel_in_gravity: [3x1] 真の衛星速度ベクトル(vx, vy, vz)
        sat_pos_ephemeris: [3x1] 軌道データ上の衛星位置(x, y, z: 月中心座標系)
        sat_vel_ephemeris: [3x1] 軌道データ上の衛星速度ベクトル(vx, vy, vz)
        user_pos: [3x1] 真のユーザー位置(x, y, z: 月中心座標系)
        user_vel: [3x1] 真のユーザー速度ベクトル(vx, vy, vz)
    戻り値
        f_residual: 推定ドップラーシフト量と観測ドップラーシフト量の残差(Hz)
        A: [1x5] ヤコビ行列
    '''

    # 推定ユーザー情報処理
    estimate_user_pos = user_data[:3]                       # 推定ユーザー位置[x_user, y_user, z_user]
    estimate_user_vel = calc_user_vel(estimate_user_pos)    # 推定ユーザー速度[vx_user, vy_user, vz_user]
    estimate_user_drift = user_data[5]                      # 推定ユーザークロックドリフト
    
    # 推定ユーザー観測ドップラーシフト量
    r_est = sat_pos_ephemeris - estimate_user_pos           # 推定衛星位置 - 推定ユーザー位置
    v_est = sat_vel_ephemeris - estimate_user_vel           # 推定衛星速度 - 推定ユーザー速度
    denom_est = np.sqrt(np.sum(r_est)**2, 1)
    estimate_los_unit = ((r_est) / denom_est)
    f_d_est = (f_carrier * np.dot(v_est, estimate_los_unit) / c) + f_carrier * estimate_user_drift
    
    # ユーザー実観測ドップラーシフト量
    r_obs = sat_pos_in_gravity - user_pos
    v_obs = sat_vel_in_gravity - user_vel
    denom_obs = np.sqrt(np.sum(r_obs)**2, 1)
    obs_los_unit = ((r_obs) / denom_obs)
    f_d_obs = (f_carrier * np.dot(v_obs, obs_los_unit) / c) + f_carrier * clock_drift

    # 周波数測定残差
    f_residual = f_d_est - f_d_obs
    
    # ヤコビ行列A
    A = []
    norm_r_est = np.linalg.norm(r_est)
    term1 = v_est / norm_r_est
    term2 = np.dot(v_est, r_est) * r_est / (norm_r_est**3)
    A[:3] = - (f_carrier / c) * (term1 - term2)
    A[4] = 0
    A[5] = f_carrier

    return f_residual, A

def main():
    load_spices()
    
    # 「重力分布影響あり」の衛星軌道情報をGMAT出力から読み取り
    base_dir_in_gravity = os.path.join(localapp, 'GMAT', 'output', 'in_gravity')
    sat_ids_in_gravity = load_sat_ids(base_dir_in_gravity, sat_name, sat_num)
    # 「重力分布影響なし」の衛星軌道情報をGMAT出力から読み取り
    base_dir_ephemeris = os.path.join(localapp, 'GMAT', 'output', 'ephemeris')
    sat_ids_ephemeris = load_sat_ids(base_dir_ephemeris, sat_name, sat_num)

    # 時刻関連情報設定
    sample_bsp_file = os.path.join(base_dir_in_gravity, print(f'{sat_name}1.bsp'))
    cover = spice.spkcov(sample_bsp_file, sat_ids_in_gravity(0), 1)
    start_time_et = cover(0)
    end_time_et = start_time_et + 24*3600
    simulation_time_et = np.arange(start_time_et, end_time_et + pos_interval, pos_interval)
    simulation_time_sec = simulation_time_et - simulation_time_et[0]
    simulation_time_hr = (simulation_time_sec - simulation_time_sec[0]) / 3600

    # 時間ごとのクロックバイアス配列を生成
    clock_bias = base_clock_bias + simulation_time_sec * clock_drift

    # 標準偏差10mの各測位間隔あたりのガウスノイズ(距離誤差)を各衛星ごとに生成
    gausian_noise = np.random.randn(sat_num, len(simulation_time_hr)) * 10

    # 各種結果行列宣言
    output_navigation_result = []
    output_f_doppler = []
    output_elevation = []
    output_gdop = []

    # 測位シミュレーションここから
    for l in range(len(simulation_time_hr)):
        et = simulation_time_et(l)
        # 初回の測位もしくは全基観測不能状態から復帰した初回の観測では初期位置を推定位置とする
        if l == 0 or any(math.isnan(estimate_usr_data)):
            estimate_usr_data = [[0], [0], [-radius_moon_m], [0], [0]]
        else:
            estimate_usr_data[:3] = output_navigation_result[:, l-1]

        sat_timing_position = []
        counter = 0
        timing_clock_bias = clock_bias(l)

        # 時刻ごとのユーザー速度
        [state, _] = spice.spkcpo('USER', et, 'IAU_MOON', 'OBSERVER', 'NONE', user_pos, 'MOON', 'IAU_MOON')
        user_vel = state[4:6]
        
        while True:
            counter += 1
            residual = np.full((2*sat_num, 1), np.nan)
            A = np.zeros((2*sat_num, 5))

            for n in range(sat_name):
                # 真の衛星の軌道データ
                sat_id_in_gravity = sat_ids_in_gravity(n)
                timing_sat_pos_in_gravity, timing_sat_vel_in_gravity = get_sat_pos_vel(et, sat_id_in_gravity, 'NONE')
                sat_timing_position[:,n] = timing_sat_pos_in_gravity

                # 時間ごとのバイアス込みの時刻情報
                estimate_bias = estimate_usr_data(4)
                estimate_et = et + clock_bias - estimate_bias

                # エフェメリス上の軌道データ
                sat_id_ephemeris = sat_ids_ephemeris(n)
                timing_sat_pos_ephemeris, timing_sat_vel_ephemeris = get_sat_pos_vel(et, sat_id_ephemeris, 'NONE')

                # 距離誤差ノイズ
                noise = gausian_noise(n, l)
                
                # 衛星仰角計算
                sat_el = elevation_masking(timing_sat_pos_in_gravity, user_pos)
                output_elevation(n, l) = sat_el

                # 衛星可視仰角判定
                if sat_el > elev_mask_angle:
                    # TOA測位計算処理
                    residual(n), A[n, :4] = calc_TOA(estimate_usr_data, timing_sat_pos_in_gravity, timing_sat_pos_ephemeris, timing_clock_bias, noise)
                    # FOA測位計算処理
                    residual(n+sat_num), A[n+sat_num, :] = calc_FOA(estimate_usr_data, timing_sat_pos_in_gravity, timing_sat_pos_ephemeris, timing_sat_vel_in_gravity, timing_sat_vel_ephemeris, user_vel)
                else:
                    residual(n) = np.nan()
                    residual(n+sat_num) = np.nan()
                    A[n, :] = np.nan()

if __name__ == '__main__':
    main()