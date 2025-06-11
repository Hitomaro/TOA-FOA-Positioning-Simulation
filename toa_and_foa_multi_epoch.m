warning('off','all');

% SPICEカーネルを呼び出し
function load_spice_kernels()
    base_dir = '.\constant_files';
    addpath('C:\Program Files\mice\lib\');
    addpath('C:\Program Files\mice\src\mice');
    cspice_furnsh(fullfile(base_dir, 'naif0012.tls'));                                  % 時刻カーネル
    cspice_furnsh(fullfile(base_dir, 'moon_pa_de421_1900-2050.bpc'));                   % 月の姿勢
    cspice_furnsh(fullfile(base_dir, 'pck00010.tpc'));                                  % 惑星定数（テキストPCK、IAU系含む）
    cspice_furnsh(fullfile(base_dir, 'moon_080317.tf'));                                % フレーム定義（IAU_MOONなど）
    cspice_furnsh(fullfile(base_dir, 'de421.bsp'));                                     % 月の軌道情報（必要に応じて）
end

% 各衛星のID番号を取得
function sat_ids = load_sat_ids(base_dir, sat_name, sat_num)
    sat_ids = zeros(1, sat_num);
    for i = 1:sat_num
        bsp_file = fullfile(base_dir, sprintf('%s%d.bsp', sat_name, i));
        cspice_furnsh(bsp_file);
        sat_ids(1, i) = cspice_spkobj(bsp_file, 1);
    end
end

% ユーザー速度計算関数
function user_vel = calc_user_vel(user_pos)
    % 入力
    %   r_moon_m: 月半径 (m)
    %   user_pos: [3x1] ユーザー位置 (x, y, z: 月中心座標系)
    % 出力
    %   user_vel: [3x1] ユーザー速度 (vx, vy, vz)

    omega_moon = 2 * pi / (27.32 * 24 * 60 * 60);                                       % 月の自転角速度（単位: rad/s）
    omega_vec = [0; 0; omega_moon];                                                     % 回転軸ベクトル z方向
    user_vel = cross(omega_vec, user_pos);                                              % 速度ベクトル = 自転による位置ベクトルとの外積
end


% 衛星 - ユーザー間仰角計算関数
function el = elevation_mask(sat_pos, user_pos)
    % 入力
    %   sat_pos: [3x1] 衛星位置 (x, y, z: 月中心座標系)
    %   user_pos: [3x1] ユーザー位置 (x, y, z: 月中心座標系)
    % 出力
    %   el: ユーザーから見た衛星仰角 (°)
    los = sat_pos - user_pos;                                                           % 視線ベクトル
    n_user = user_pos / norm(user_pos);                                                 % ユーザー天頂方向
    r_los = los / norm(los);                                                            % 単位視線ベクトル

    cos_el = dot(n_user, r_los);                                                        % 仰角のcos
    el = asind(cos_el);                                                                 % 仰角(°)
end


function [pos, vel] = get_sat_pos_vel(et, sat_name, abcorr)
    % 入力
    %   et: SPICE時刻（ET秒）
    %   sat_name: 衛星名文字列 (例: 'MOONSAT1')
    %   abcorr: 補正方法 (例: 'NONE'、'LT+S')
    % 出力
    %   pos: SPICE時刻における衛星位置 (m)
    %   vel: SPICE時刻における衛星の速度ベクトル (m/s)

    sat_id = int2str(sat_name);
    [state_j2000, ~] = cspice_spkezr(sat_id, et, 'J2000', abcorr, 'MOON');              % CSPICEによるアルマナックからの座標(J2000)取得処理

    rot_mat = cspice_sxform('J2000', 'IAU_MOON', et);                                   % 座標変換行列取得（J2000 → IAU_MOON)
    state = rot_mat * state_j2000;                                                      % 座標変換

    % stateは6要素のベクトル [x;y;z;vx;vy;vz] (km, km/s)
    pos = state(1:3) * 1000;                                                            % 位置(m)
    vel = state(4:6) * 1000;                                                            % 速度ベクトル(m/s)
end


% TOAヤコビ行列A計算関数
function [residual, A] = calc_TOA(estimate_user_data, timing_sat_pos, user_pos, noise, c, clock_bias)
    % 入力
    %   estimate_user_data: [1x5] 推定ユーザーデータ ([1~3]ユーザー位置、[4]クロックバイアス、[5]クロックドリフト量(この関数では使わない))
    %   timing_sat_pos: [3x1] 衛星位置 (x, y, z: 月中心座標系)
    %   user_pos: [3x1] 真のユーザー位置 (x, y, z: 月中心座標系)
    %   noise: 距離誤差ノイズ (m)
    %   c: 光速 (3.0 x 10^8 m/s)
    %   clock_bias: 真のクロックバイアス (秒)
    % 出力
    %   residual: 疑似距離 - 観測距離残差 (m)
    %   A: [1x4] ヤコビ行列

    % 推定情報処理
    estimate_user_pos = estimate_user_data(1:3);                                        % 推定ユーザー位置
    estimate_clock_bias = estimate_user_data(4);                                        % 推定クロックバイアス
    estimate_range = norm(timing_sat_pos - estimate_user_pos);                          % 推測位置を用いた衛星-ユーザー間疑似距離
    estimate_pseudo_range = estimate_range + c * estimate_clock_bias;                   % 推定衛星-ユーザー間距離

    % 実観測情報処理                                                                        
    true_range = norm(timing_sat_pos - user_pos);                                       % 真の衛星-ユーザー間距離
    pseudo_range = true_range + c * clock_bias;                                         % ユーザーが観測する衛星-ユーザー間距離
    obs_range = pseudo_range + noise;                                                   % 距離ノイズを含んだユーザー観測衛星 - ユーザー間距離
    residual = estimate_pseudo_range - obs_range;                                       % 測定距離残差

    % ヤコビ行列Aを追加
    A(1:3) = (estimate_user_pos - timing_sat_pos)' / estimate_range;
    A(4) = c;
end


% FOAヤコビ行列A計算関数
function [f_residual, A, f_d_obs] = calc_FOA(f_carrier, estimate_user_data, timing_sat_pos, timing_sat_vel, user_pos, user_vel, drift, c)
    % 入力
    %   f_carrier: 搬送波周波数 (Hz)
    %   estimate_user_data: [1x5] 推定ユーザーデータ ([1~3]ユーザー位置、[4]クロックバイアス(この関数では使わない)、[5]クロックドリフト量)
    %   timing_sat_pos: [3x1] 衛星位置 (x, y, z: 月中心座標系)
    %   timing_sat_vel: [3x1] 衛星速度ベクトル (vx, vy, vz)
    %   user_pos: [3x1] 真のユーザー位置 (x, y, z: 月中心座標系)
    %   user_vel: [3x1] 真のユーザー速度 (vx, vy, vz)
    %   drift: 真のクロックドリフト量
    %   c: 光速
    % 出力
    %   f_residual: 推定ドップラーシフト量と観測ドップラーシフト量の残差 (Hz)
    %   A: [1x5] ヤコビ行列

    % 推定ユーザー位置
    estimate_user_pos = estimate_user_data(1:3);
    estimate_user_vel = calc_user_vel(estimate_user_pos);
    estimate_user_drift = estimate_user_data(5);

    % 推定ユーザー観測ドップラーシフト量
    r_est = timing_sat_pos - estimate_user_pos;
    v_est = timing_sat_vel - estimate_user_vel;
    denom_est = sqrt(sum((r_est).^2, 1));
    estimate_los_unit = ((r_est) ./ denom_est);
    f_d_est = (f_carrier * dot(v_est, estimate_los_unit) / c) + f_carrier * estimate_user_drift;

    % ユーザー観測ドップラーシフト量
    r_obs = timing_sat_pos - user_pos;
    v_obs = timing_sat_vel - user_vel;
    denom_obs = sqrt(sum((r_obs).^2, 1));
    obs_los_unit = ((r_obs) ./ denom_obs);
    f_d_obs = (f_carrier * dot(v_obs, obs_los_unit) / c) + f_carrier * drift;
    
    % 周波数測定残差
    f_residual = f_d_est - f_d_obs;
    
    % ヤコビ行列Aを追加
    norm_r_est = norm(r_est);
    term1 = v_est / norm_r_est;
    term2 = dot(v_est, r_est) * r_est / (norm_r_est^3);
    A(1:3) = - (f_carrier / c) * (term1 - term2);
    A(4) = 0;
    A(5) = f_carrier;
end


% ------------------------ ここからメイン処理 ------------------------
clear;

load_spice_kernels();

% 各種定数
r_moon_km = 1737.4;                                                                     % 月半径(km)
r_moon_m = r_moon_km * 1000;                                                            % 月半径(m)
c = 3.0 * 10^8;                                                                         % 光速(m/s)
f_carrier = 2483.5;                                                                     % LNSS搬送波周波数(MHz)
f_carrier = f_carrier * 10^6;                                                           % LNSS搬送波周波数(Hz)
drift = 1.0e-6;                                                                         % クロックドリフト
base_clock_bias = 0;                                                                    % 基準クロックバイアス

positioning_interval = 10;                                                              % 測位時間間隔(秒)
number_of_positioning = 3;                                                              % 1回の測位あたりの観測回数
threshold = 1e-4;                                                                       % 十分収束したと判断する前回測位結果との位置誤差の基準値
max_iter = 10;                                                                          % 最小二乗法の最高反復計算処理回数
elev_mask_angle = 10;                                                                   % 仰角マスク角度(度)
sat_num = 6;                                                                            % 衛星数

% 衛星位置をGMAT出力から読み込み
sat_name = 'sat';                                                                       % エフェメリスファイルの衛星名
base_dir = '.\orbit_data\ELFO_true_sats';                                               % エフェメリスファイルがあるフォルダ名
sat_ids = load_sat_ids(base_dir, sat_name, sat_num);                                    % 各衛星のIDを取得

% 衛星搭載クロック時刻
bsp_file = fullfile(base_dir, sprintf('%s1.bsp', sat_name));
cover = cspice_spkcov(bsp_file, sat_ids(1), 1);
start_time_et = cover(1)+100;
end_time_et = 24*3600;%cover(end);
simulation_time_et = start_time_et:positioning_interval:end_time_et;
simulation_time_sec = (simulation_time_et - simulation_time_et(1));
simulation_time_hours = ((simulation_time_sec - simulation_time_sec(1)) / 3600);

all_clock_bias = base_clock_bias + simulation_time_sec .* drift;                        % 時間ごとのクロックバイアス値生成

% 標準偏差10mの時間当たりのガウスノイズ(距離誤差)を各衛星ごとに生成
rng(53);
gausian_noise = randn(sat_num, length(simulation_time_hours)) .* 10;

% ユーザー位置
x = 0;
y = 0;
z = -sqrt(r_moon_m^2 - x^2 - y^2);
user_pos = [x; y; z];           % [1x3] 真のユーザー位置 (m)

% ユーザー速度
user_vel = calc_user_vel(user_pos);


% 各種結果配列を宣言
output_time_estimate_pos = nan(3,length(simulation_time_hours));                        % 時間ごとの測位結果
output_f_doppler = nan(sat_num, length(simulation_time_hours));                         % 時間ごとの各衛星から発射された電波のドップラーシフト量
usr_sat_elevation = zeros(sat_num, length(simulation_time_hours));                      % ユーザーから見た衛星の仰角情報
% gdop = zeros(1, length(simulation_time_hours));                                       % 時間ごとのGDOP値


% 衛星測位プロセス ここから
for l=1:length(simulation_time_hours)
    et = simulation_time_et(l);
    % 初回の測位もしくは測位不能状態から復帰した初回の観測では初期位置を推定位置とする
    if l == 1 || any(isnan(estimate_user_data))
        estimate_user_data = [10;10;-r_moon_m+10;0;0];
        %estimate_user_data=[0;0;-r_moon_m;0;0];
    else
        % それ以外では前回の測位結果を推定位置とする
        estimate_user_data(1:3) = output_time_estimate_pos(:, l-1);
    end
    sat_timing_position = zeros(3, sat_num);
    iter = 0;
    while true
        iter = iter + 1;
        residual = nan((2*sat_num)*3, 1);
        A = zeros((2*sat_num)*number_of_positioning, 5);
        et_dash = et - (positioning_interval*number_of_positioning);
        for m=1:number_of_positioning
            A_dash = zeros(2*sat_num, 5);
            residual_dash = nan((2*sat_num), 1);
            for n=1:sat_num
                sat_id = sat_ids(1, n);
                [timing_sat_pos, timing_sat_vel] = get_sat_pos_vel(et_dash, sat_id, 'NONE');
                sat_timing_position(:, n) = timing_sat_pos;

                noise = gausian_noise(n,l);                                                 % 距離誤差ノイズ
                
                % ユーザー - 衛星間仰角計算
                sat_el = elevation_mask(timing_sat_pos, user_pos);
                usr_sat_elevation(n,l) = sat_el;

                % 衛星可視仰角条件判定
                if sat_el > elev_mask_angle
                    % 衛星ごとのTOA計算処理
                    clock_bias = all_clock_bias(l);
                    [residual_dash(n), A_dash(n, 1:4)] = calc_TOA(estimate_user_data, timing_sat_pos, user_pos, noise, c, clock_bias);
                    % 衛星ごとのFOA計算処理
                    [residual_dash(n+sat_num), A_dash(n+sat_num, :), output_f_doppler(n, l)] = calc_FOA(f_carrier, estimate_user_data, timing_sat_pos, timing_sat_vel, user_pos, user_vel, drift, c);
                else
                    % 仰角が基準以下なら観測結果を不可視(NaN)に上書き
                    residual_dash(n) = NaN;
                    residual_dash(n+sat_num) = NaN;
                    A_dash(n, :) = [NaN, NaN, NaN, NaN, NaN];
                    A_dash(n+sat_num, :) = [NaN, NaN, NaN, NaN, NaN];
                end
            end
            array_start = ((sat_num*2)*m) - ((sat_num*2)-1);
            array_end = (sat_num*2)*m;
            A(array_start:array_end, :) = A_dash;
            residual(array_start:array_end, :) = residual_dash;
            et_dash = et_dash + positioning_interval;
        end
        
        vaild_idx = ~any(isnan(A), 2);                                                  % A行列にNaN(ユーザーから観測不可能)が入っていないものを抽出
        A_vaild = A(vaild_idx, :);                                                      % 観測可能な衛星のみに絞ってA行列を再構成
        residual_vaild = residual(vaild_idx);                                           % 観測可能な衛星のみに絞って測定残差行列を再構成

        total_rank = rank(A_vaild);                                                     % 観測行列(A行列)の階数を計算
        
        % 可観測衛星が3基ない場合、もしくはA行列の階数が5未満の場合測位不能とみなす
        if length(residual_vaild) < 3 || total_rank < 5
            estimate_user_data = [NaN, NaN, NaN, NaN, NaN];
            break;
        else
            X = (A_vaild'*A_vaild) \ (A_vaild' * residual_vaild);                       % 最小二乗法計算(疑似距離からのx, y, zの誤差を出す)
            estimate_user_data = estimate_user_data - X;                                % 推定位置に最小二乗法で出した各方向の誤差を入れて推定位置を修正

            % 最小二乗法で出した推定位置の誤差が「基準値以下」もしくは繰り返し回数が「基準以下」なら反復計算処理を終了
            if norm(X) < threshold || iter >= max_iter
                break;
            end
        end
    end
    output_time_estimate_pos(:, l) = estimate_user_data(1:3);
end
% ------------------------ ここまで ------------------------

% 実際の位置と測位位置との距離 (m)
est_err = vecnorm(user_pos - output_time_estimate_pos);
% N%信頼区間を計算
sorted_est_err = sort(est_err);
CI_95 = prctile(sorted_est_err, 95);
fprintf("95%%信頼区間: %.2f(m)\n", CI_95);


figure;
scatter3(user_pos(1,:), user_pos(2,:), user_pos(3,:), 200, 'blue', 'filled', 'd');
hold on;
scatter3(output_time_estimate_pos(1,:), output_time_estimate_pos(2,:), output_time_estimate_pos(3,:), 20, simulation_time_hours, 'filled');
colorbar;
colormap(jet);
hold off;
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
axis equal;
grid on;

figure;
plot(simulation_time_hours, est_err);
ylabel("estimate posinoning error(m)");
xlabel("time (hr)");
grid on;