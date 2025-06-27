import spiceypy as spice
import numpy as np
import matplotlib
import os
import math

localapp = os.environ['LOCALAPPDATA']

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

def main():
    load_spices()
    
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
        
        while True:
            counter += 1
            risidual = np.full((2*sat_num, 1), np.nan)

            for n in range(sat_name):
                # 衛星ごとのTOA計算処理
                sat_id_in_gravity = sat_ids_in_gravity(n)
                timing_sat_pos_in_gravity, timing_sat_vel_in_gravity = get_sat_pos_vel(et, sat_id_in_gravity)

if __name__ == '__main__':
    main()