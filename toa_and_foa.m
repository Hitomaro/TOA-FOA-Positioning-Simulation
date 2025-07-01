warning('off','all');

% SPICEカーネルを呼び出し(読み飛ばしてOK)
function load_spice_kernels()
    base_dir = '.\constant_files';
    addpath(fullfile(getenv('LOCALAPPDATA'), 'mice', 'lib'));
    addpath(fullfile(getenv('LOCALAPPDATA'), 'mice', 'src', 'mice'));
    cspice_furnsh(fullfile(base_dir, 'naif0012.tls'));                  % 時刻カーネル
    cspice_furnsh(fullfile(base_dir, 'moon_pa_de421_1900-2050.bpc'));   % 月の姿勢
    cspice_furnsh(fullfile(base_dir, 'pck00010.tpc'));                  % 惑星定数（テキストPCK、IAU系含む）
    cspice_furnsh(fullfile(base_dir, 'moon_080317.tf'));                % フレーム定義（IAU_MOONなど）
    cspice_furnsh(fullfile(base_dir, 'de421.bsp'));                     % 月の軌道情報（必要に応じて）
end

% 各衛星のID番号を取得(読み飛ばしてOK)
function sat_ids = load_sat_ids(base_dir, sat_name, sat_num)
    sat_ids = zeros(1, sat_num);
    for i = 1:sat_num
        bsp_file = fullfile(base_dir, sprintf('%s%d.bsp', sat_name, i));
        cspice_furnsh(bsp_file);
        sat_ids(1, i) = cspice_spkobj(bsp_file, 1);
    end
end

% ユーザー速度計算関数(多分読み飛ばしてOK)
function user_vel = calc_user_vel(user_pos)
    % 入力
    %   r_moon_m: 月半径 (m)
    %   user_pos: [3x1] ユーザー位置 (x, y, z: 月中心座標系)
    % 出力
    %   user_vel: [3x1] ユーザー速度 (vx, vy, vz)

    omega_moon = 2 * pi / (27.32 * 24 * 60 * 60);    % 月の自転角速度（単位: rad/s）
    omega_vec = [0; 0; omega_moon];                  % 回転軸ベクトル z方向
    user_vel = cross(omega_vec, user_pos);           % 速度ベクトル = 自転による位置ベクトルとの外積
end


% 衛星 - ユーザー間仰角計算関数(多分読み飛ばしてOK)
function el = elevation_mask(sat_pos, user_pos)
    % 入力
    %   sat_pos: [3x1] 衛星位置 (x, y, z: 月中心座標系)
    %   user_pos: [3x1] ユーザー位置 (x, y, z: 月中心座標系)
    % 出力
    %   el: ユーザーから見た衛星仰角 (°)

    los = sat_pos - user_pos;               % 視線ベクトル
    n_user = user_pos / norm(user_pos);     % ユーザー天頂方向
    r_los = los / norm(los);                % 単位視線ベクトル
    cos_el = dot(n_user, r_los);            % 仰角のcos
    el = asind(cos_el);                     % 仰角(°)
end


% SPKファイル(軌道情報ファイル)から時刻ごとの衛星位置と衛星速度ベクトルを取り出す関数(読み飛ばしてOK)
function [pos, vel] = get_sat_pos_vel(et, sat_name, abcorr)
    % 入力
    %   et: SPICE時刻（ET秒）
    %   sat_name: 衛星名文字列 (例: 'MOONSAT1')
    %   abcorr: 補正方法 (例: 'NONE'、'LT+S')
    % 出力
    %   pos: SPICE時刻における衛星位置 (m)
    %   vel: SPICE時刻における衛星の速度ベクトル (m/s)

    sat_id = int2str(sat_name);
    [state_j2000, ~] = cspice_spkezr(sat_id, et, 'J2000', abcorr, 'MOON');      % SPICEによるアルマナックからの座標(J2000)取得処理

    rot_mat = cspice_sxform('J2000', 'IAU_MOON', et);                           % 座標変換行列取得（J2000 → IAU_MOON)
    state = rot_mat * state_j2000;                                              % 座標変換

    % stateは6要素のベクトル [x;y;z;vx;vy;vz] (km, km/s)
    pos = state(1:3) * 1000;                                                    % 位置(m)
    vel = state(4:6) * 1000;                                                    % 速度ベクトル(m/s)
end


% TOAヤコビ行列A計算関数
function [residual, A] = calc_TOA(estimate_user_data, timing_sat_pos, estimate_timing_sat_pos, user_pos, noise, c, clock_bias)
    % 入力
    %   estimate_user_data: [1x5] 推定ユーザーデータ ([1~3]ユーザー位置、[4]クロックバイアス、[5]クロックドリフト量(この関数では使わない))
    %   timing_sat_pos: [3x1] 真の衛星位置 (x, y, z: 月中心座標系)
    %   estimate_timing_sat_pos: [3x1] 推定衛星位置 (x, y, z: 月中心座標系)
    %   user_pos: [3x1] 真のユーザー位置 (x, y, z: 月中心座標系)
    %   noise: 距離誤差ノイズ (m)
    %   c: 光速 (3.0 x 10^8 m/s)
    %   clock_bias: 真のクロックバイアス (秒)
    % 出力
    %   residual: 疑似距離 - 観測距離残差 (m)
    %   A: [1x4] ヤコビ行列

    % 推定情報処理
    estimate_user_pos = estimate_user_data(1:3);                            % 推定ユーザー位置[x_user, y_user, z_user]
    estimate_clock_bias = estimate_user_data(4);                            % 推定クロックバイアス
    estimate_range = norm(estimate_timing_sat_pos - estimate_user_pos);     % 推測位置を用いた衛星-ユーザー間疑似距離
    estimate_pseudo_range = estimate_range + c * estimate_clock_bias;       % 推定衛星-ユーザー間距離

    % 実観測情報処理                                                                        
    true_range = norm(timing_sat_pos - user_pos);                           % 真の衛星-ユーザー間距離
    pseudo_range = true_range + c * clock_bias;                             % ユーザーが観測する衛星-ユーザー間距離
    obs_range = pseudo_range + noise;                                       % 距離ノイズを含んだユーザー観測衛星 - ユーザー間距離
    residual = estimate_pseudo_range - obs_range;                           % 測定距離残差

    % ヤコビ行列Aを追加
    A(1:3) = (estimate_user_pos - timing_sat_pos)' / estimate_range;
    A(4) = c;
end


% FOAヤコビ行列A計算関数
function [f_residual, A, f_d_obs] = calc_FOA(f_carrier, estimate_user_data, timing_sat_pos, timing_sat_vel, estimate_timing_sat_pos, estimate_timing_sat_vel, user_pos, user_vel, drift, c)
    % 入力
    %   f_carrier: 搬送波周波数 (Hz)
    %   estimate_user_data: [1x5] 推定ユーザーデータ ([1~3]ユーザー位置、[4]クロックバイアス(この関数では使わない)、[5]クロックドリフト量)
    %   timing_sat_pos: [3x1] 真の衛星位置 (x, y, z: 月中心座標系)
    %   timing_sat_vel: [3x1] 真の衛星速度ベクトル (vx, vy, vz)
    %   estimate_timing_sat_pos: [3x1] 推定衛星位置 (x, y, z: 月中心座標系)
    %   estimate_timing_sat_vel: [3x1] 真の衛星速度ベクトル (vx, vy, vz)
    %   user_pos: [3x1] 真のユーザー位置 (x, y, z: 月中心座標系)
    %   user_vel: [3x1] 真のユーザー速度 (vx, vy, vz)
    %   drift: 真のクロックドリフト量
    %   c: 光速
    % 出力
    %   f_residual: 推定ドップラーシフト量と観測ドップラーシフト量の残差 (Hz)
    %   A: [1x5] ヤコビ行列

    % 推定ユーザー情報を取得
    estimate_user_pos = estimate_user_data(1:3);            % 推定ユーザー位置[x_user, y_user, z_user]
    estimate_user_vel = calc_user_vel(estimate_user_pos);   % 推定ユーザー速度[vx_user, vy_user, vz_user]
    estimate_user_drift = estimate_user_data(5);            % 推定ユーザークロックドリフト

    % 推定ユーザー観測ドップラーシフト量
    r_est = estimate_timing_sat_pos - estimate_user_pos;    % 推定衛星 - ユーザー間距離[x_sat-x_user, y_sat-y_user, z_sat-z_user]
    v_est = estimate_timing_sat_vel - estimate_user_vel;    % 推定衛星 - ユーザー間相対速度[vx_sat-vx_user, vy_sat-vy_user, vz_sat-vz_user]
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

% DOP計算関数(多分読み飛ばしてOK)
function GDOP = calc_DOP(sat_num, sat_position, user_pos, elev_mask_angle)
    G = zeros(sat_num, 4);
    for n=1:sat_num
        sat_pos = sat_position(:, n);
        if elevation_mask(sat_pos, user_pos) > elev_mask_angle
            r = norm(sat_pos - user_pos);
            G(n, 1:3) = (sat_pos - user_pos) / r;
            G(n, 4) = 1;
        else
            G(n, :) = NaN;
        end
    end
    G = G(~any(isnan(G), 2), :);
    Q = pinv(G' * G);
    GDOP = sqrt(trace(Q));
end


% ------------------------ ここからメイン処理 ------------------------
clear;

load_spice_kernels();

% 各種定数
r_moon_km = 1737.4;              % 月半径(km)
r_moon_m = r_moon_km * 1000;     % 月半径(m)
c = 3.0 * 10^8;                  % 光速(m/s)
f_carrier = 2483.5;              % LNSS搬送波周波数(MHz)
f_carrier = f_carrier * 10^6;    % LNSS搬送波周波数(Hz)
drift = 1.0e-6;                  % クロックドリフト
base_clock_bias = 0;             % 基準クロックバイアス

positioning_interval = 10;       % 測位時間間隔(秒)
threshold = 1e-4;                % 十分収束したと判断する前回測位結果との位置誤差の基準値
max_iter = 10;                   % 最小二乗法の最高反復計算処理回数
elev_mask_angle = 15;            % 仰角マスク角度(度)
sat_num = 6;                     % 衛星数


% 重力分布の影響を受けている衛星位置をGMAT出力から読み込み
sat_name_true = 'sat';                                                                      % ファイルの衛星名
base_dir_true = fullfile(getenv("LOCALAPPDATA"), 'GMAT', 'output', 'in_gravity');           % 軌道ファイルがあるフォルダ名
sat_ids_true = load_sat_ids(base_dir_true, sat_name_true, sat_num);                         % 各衛星のIDを取得

% 重力分布の影響を受けない理論上の衛星位置をGMAT出力から読み込み
sat_name_ephemeris = 'sat';                                                                 % ファイルの衛星名
base_dir_ephemeris = fullfile(getenv("LOCALAPPDATA"), 'GMAT', 'output', 'in_gravity');       % 軌道ファイルがあるフォルダ名
sat_ids_ephemeris = load_sat_ids(base_dir_ephemeris, sat_name_ephemeris, sat_num);          % 各衛星のIDを取得

% 衛星搭載クロック時刻 (本質ではないので読み飛ばしてOK)
% エフェメリス上の開始時刻と終了時刻を取得して秒や時間に変換しているだけ
bsp_file = fullfile(base_dir_true, sprintf('%s1.bsp', sat_name_true));
cover = cspice_spkcov(bsp_file, sat_ids_true(1), 1);
start_time_et = cover(1)+100;
end_time_et = start_time_et + 24*3600;%cover(end);
simulation_time_et = start_time_et:positioning_interval:end_time_et;
simulation_time_sec = (simulation_time_et - simulation_time_et(1));
simulation_time_hours = ((simulation_time_sec - simulation_time_sec(1)) / 3600);

% 時間ごとのクロックバイアス値生成
all_clock_bias = base_clock_bias + simulation_time_sec .* drift;

% 標準偏差10mの時間当たりのガウスノイズ(距離誤差)を各衛星ごとに生成
% 下のrng(i);のiに当たる値を固定にすると常に同じ乱数が出てくる(再現性を持たせるため)
rng(91);
gausian_noise = randn(sat_num, length(simulation_time_hours)) .* 10;

% ユーザー位置
x = 0;
y = 0;
z = -sqrt(r_moon_m^2 - x^2 - y^2);
user_pos = [x; y; z];           % [1x3] 真のユーザー位置 (m)

% ユーザー速度
user_vel = calc_user_vel(user_pos);


% 各種結果配列を宣言
output_time_estimate_pos = zeros(3,length(simulation_time_hours));      % 時間ごとの測位結果
output_f_doppler = nan(sat_num, length(simulation_time_hours));         % 時間ごとの各衛星から発射された電波のドップラーシフト量
usr_sat_elevation = zeros(sat_num, length(simulation_time_hours));      % ユーザーから見た衛星の仰角情報
gdop = zeros(1, length(simulation_time_hours));                         % 時間ごとのGDOP値


% 衛星測位プロセス ここから
for l=1:length(simulation_time_hours)
    et = simulation_time_et(l);
    % 初回の測位もしくは測位不能状態から復帰した初回の観測では初期位置を推定位置とする
    if l == 1 || any(isnan(estimate_user_data))
        %estimate_user_data = [10;10;-r_moon_m+10;0;0];
        estimate_user_data=[0;0;-r_moon_m;0;0];
    else
        % それ以外では前回の測位結果を推定位置とする
        estimate_user_data(1:3) = output_time_estimate_pos(:, l-1);
    end
    sat_timing_position = zeros(3, sat_num);    % グラフで使う用の各時刻における衛星位置を入れるための配列(読み飛ばしてOK)
    iter = 0;                                   % 最小二乗法の繰り返し回数カウント用変数
    clock_bias = all_clock_bias(l);
    while true
        iter = iter + 1;
        residual = nan(2*sat_num, 1);
        A = zeros(2*sat_num, 5);
        for n=1:sat_num
            % 衛星ごとのTOA計算処理
            sat_id_true = sat_ids_true(1, n);
            [timing_sat_pos_true, timing_sat_vel_true] = get_sat_pos_vel(et, sat_id_true, 'NONE');
            sat_timing_position(:, n) = timing_sat_pos_true;
            
            estimate_bias = estimate_user_data(4);
            estimate_et = et + clock_bias - estimate_bias;

            % エフェメリス上の軌道データ
            sat_id_ephemeris = sat_ids_ephemeris(1,n);
            [estimate_timing_sat_pos, estimate_timing_sat_vel] = get_sat_pos_vel(estimate_et, sat_id_ephemeris, 'NONE');

            % 距離誤差ノイズ
            noise = gausian_noise(n,l);

            % ユーザー - 衛星間仰角計算
            sat_el = elevation_mask(timing_sat_pos_true, user_pos);
            usr_sat_elevation(n,l) = sat_el;

            % 衛星可視仰角条件判定
            if sat_el > elev_mask_angle
                [residual(n), A(n, 1:4)] = calc_TOA(estimate_user_data, timing_sat_pos_true, estimate_timing_sat_pos,  user_pos, noise, c, clock_bias);
                % 衛星ごとのFOA計算処理
                [residual(n+sat_num), A(n+sat_num, :), output_f_doppler(n, l)] = calc_FOA(f_carrier, estimate_user_data, timing_sat_pos_true, timing_sat_vel_true, estimate_timing_sat_pos, estimate_timing_sat_vel, user_pos, user_vel, drift, c);
            else
                % 仰角が基準以下なら観測結果を不可視(NaN)に上書き
                residual(n) = NaN;
                residual(n+sat_num) = NaN;
                A(n, :) = [NaN, NaN, NaN, NaN, NaN];
                A(n+sat_num, :) = [NaN, NaN, NaN, NaN, NaN];
            end
        end
        
        % 行列の再構成
        vaild_idx = ~any(isnan(A), 2);              % A行列にNaN(ユーザーから観測不可能)が入っていないものを抽出
        A_vaild = A(vaild_idx, :);                  % 観測可能な衛星のみに絞ってA行列を再構成
        residual_vaild = residual(vaild_idx);       % 観測可能な衛星のみに絞って測定残差行列を再構成
        
        % 観測行列(A行列)の階数を計算
        total_rank = rank(A_vaild);
        
        % 可観測衛星が3基ない場合、もしくはA行列の階数が5未満の場合測位不能とみなす
        if length(residual_vaild) < 3 || total_rank < 5
            estimate_user_data = [NaN, NaN, NaN, NaN, NaN];
            break;
        else
            X = (A_vaild'*A_vaild) \ (A_vaild' * residual_vaild);       % 最小二乗法計算(疑似距離からのx, y, zの誤差を出す)
            estimate_user_data = estimate_user_data - X;                % 推定位置に最小二乗法で出した各方向の誤差を入れて推定位置を修正

            % 最小二乗法で出した推定位置の誤差が「基準値以下」もしくは繰り返し回数が「基準以下」なら反復計算処理を終了
            if norm(X) < threshold || iter >= max_iter
                break;
            end
        end
    end
    if length(residual_vaild) < 3 || total_rank < 5
        gdop(1, l) = NaN;
    else
        gdop(1, l) = calc_DOP(sat_num, sat_timing_position, user_pos, elev_mask_angle);
    end
    output_time_estimate_pos(:,l) = estimate_user_data(1:3);
end
% ここまで

% 実際の位置と測位位置との距離 (m)
est_err = vecnorm(user_pos - output_time_estimate_pos);

sorted_est_err = sort(est_err);
CI_95 = prctile(sorted_est_err, 95);
fprintf("95%%信頼区間: %.2f(m)\n", CI_95);

usr_sat_elevation(usr_sat_elevation < elev_mask_angle) = NaN;

figure;
scatter3(user_pos(1,:), user_pos(2,:), user_pos(3,:), 200, 'blue', 'filled', 'd');
hold on;
scatter3(output_time_estimate_pos(1,:), output_time_estimate_pos(2,:), output_time_estimate_pos(3,:), 20, simulation_time_hours, 'filled');
colorbar;
colormap(jet);
hold off;
xlabel('x(m)');
ylabel('y(m)')
zlabel('z(m)');
axis equal;
grid on;

% figure;
% plot(simulation_time_hours, est_err);
% xlabel("time (hour)");
% ylabel("estimate position error (m)");
% grid on;

%figure;
%yyaxis left
%plot(simulation_time_hours, usr_sat_elevation);
%ylabel("user - sat elevation (degree)");
%yyaxis right
%plot(simulation_time_hours, output_f_doppler);
%ylabel("doppler shift amount (Hz)")
%xlabel("time (hour)");
%grid on;

% 時間当たりの仰角及び測位誤差変化グラフ
figure;
yyaxis left
plot(simulation_time_hours, usr_sat_elevation);
ylabel("user - sat elevation (degree)");
yyaxis right
plot(simulation_time_hours, est_err);
ylabel("user position error (m)")
xlabel("time (hr)");
grid on;

% 時間当たりのドップラーシフト量及び測位誤差変化グラフ
figure;
yyaxis right
plot(simulation_time_hours, output_f_doppler);
hold off;
ylabel('doppler shift amount (Hz)');
yyaxis left
plot(simulation_time_hours, est_err);
ylabel("estimate position error (m)");
grid on;
xlabel('time (hr)');
grid on;

% 時間当たりのGDOP及び測位誤差変化グラフ
figure;
yyaxis left
plot(simulation_time_hours, gdop);
ylabel("GDOP");
yyaxis right
plot(simulation_time_hours, est_err);
ylabel("user position error (m)")
xlabel("time (hr)");
grid on;

% 時間当たりのGDOP変化グラフ
figure;
plot(simulation_time_hours, gdop);
ylabel("GDOP")
xlabel("time (hr)");
grid on;

data1 = [simulation_time_hours; est_err];
writematrix(data1, 'output.xlsx');