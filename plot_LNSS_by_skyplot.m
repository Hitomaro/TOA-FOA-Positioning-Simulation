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
function sat_ids = load_sats(base_dir, sat_name, sat_num)
    sat_ids = zeros(1, sat_num);
    for i = 1:sat_num
        bsp_file = fullfile(base_dir, sprintf('%s%d.bsp', sat_name, i));
        cspice_furnsh(bsp_file);
        sat_ids(1, i) = cspice_spkobj(bsp_file, 1);
    end
end


% 各種定数
r_moon_km = 1737.4;              % 月半径(km)
r_moon_m = r_moon_km * 1000;     % 月半径(m)
r_moon_flatting = 0.0;

base_dir = fullfile(getenv("LOCALAPPDATA"), 'GMAT', 'output', 'in_gravity');
sat_name = "sat";
sat_num = 6;

load_spice_kernels()
sat_ids = load_sats(base_dir, sat_name, sat_num);

% ユーザー位置設定
lat = deg2rad(-90);     % 緯度
lon = deg2rad(0);       % 経度
alt = 0;                % 高度
user_pos = cspice_georec(lon, lat, alt, 'MOON', r_moon_km, r_moon_flatting);

% 時刻設定
cover = cspice_spkcov(bsp_file, sat_ids(1), 1);
start_time_et = cover(1) + 100;
input_time = str2double(input('Input satellite position time(hr) > ', 's'));
show_time_et = start_time_et + (input_time * 3600);

% 位置描画
figure;
hold on;
colors = lines(sat_num);

for s = 1:sat_num
    [state, ~] = spice_spkpos(sat_ids(1, s), show_time_et, 'MOON', 'LT+S', user_pos);
    los = state;

    % 座標変換行列R

end